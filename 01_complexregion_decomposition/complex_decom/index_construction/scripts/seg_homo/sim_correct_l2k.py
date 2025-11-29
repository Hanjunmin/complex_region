import itertools
import pandas as pd
import pandas as pd
import re
from Bio.Seq import Seq
import pickle
from kmerhash import kmerhasher, hashes2seq
import json
import multiprocessing as mp
import pysam


import argparse
parser = argparse.ArgumentParser(description='SD-Sim')
parser.add_argument('--o', help='output file')
parser.add_argument('--r', help='reference path')
parser.add_argument('--p', help='processes of running')
parser.add_argument('--d1', help='inter-seg kmer database')
parser.add_argument('--d2', help='segment database')
args = parser.parse_args()


def group_columns_by_values(dfa):
    groups = {}
    for col in dfa.columns:
        col_values = tuple(dfa[col])
        if col_values in groups:
            groups[col_values].append(col)
        else:
            groups[col_values] = [col]
    return groups

def kmer_cal(seq,size):
    kmers = [str(seq[i:i+size]) for i in range(len(seq) - size + 1)]
    kmers = list(filter(lambda kmer: 'N' not in kmer, kmers))
    rkmers = [str(seq[i:i+size].reverse_complement()) for i in range(len(seq) - size + 1)]
    rkmers = list(filter(lambda kmer: 'N' not in kmer, rkmers))
    result = [min(kmer1, kmer2, key=lambda x: kmerhasher(x, size)) for kmer1, kmer2 in zip(kmers, rkmers)]
    result= {x.lower() for x in result}
    return(result)


def refine_conti(replacements):
    savekey=replacements.keys()
    sal={}
    delal=[]
    for item in savekey:
        #print(item)
        if item in delal:
            continue
        havecor="FASLE"
        valueitem=re.findall(r'([><]\d+)',item)
        for key in savekey:
            if key==item or key in delal:
                continue
            valu1eitem=re.findall(r'([><]\d+)', key)
            if valu1eitem[0]==valueitem[-1]:
                sal[item+"|"+key]=item+"".join(valu1eitem[1:])
                #print(sal)
                delal.extend([key,item])
                havecor="TRUE"
                break
            if valu1eitem[-1]==valueitem[0]:
                sal[key+"|"+item]=key+"".join(valueitem[1:])
                #print(sal)
                delal.extend([key,item])
                havecor="TRUE"
                break
        #print(havecor)
        if havecor=="FASLE":
            sal[item]="".join(valueitem)
    lastsal=sal
    while delal!=[]:
        delal=[]
        lastsal=sal
        sal={}
        for item,itemvalue in lastsal.items():
            #print(item)
            #print(delal)
            havecor="FASLE"
            if item in delal:
                continue
            valueitem=re.findall(r'([><]\d+)',itemvalue)
            for key,keyvalue in lastsal.items():
                if key==item or key in delal:
                    continue
                valu1eitem=re.findall(r'([><]\d+)', keyvalue)
                if valu1eitem[0]==valueitem[-1]:
                    sal[item+"|"+key]=itemvalue+"".join(valu1eitem[1:])
                    delal.extend([key,item])
                    havecor="TRUE"
                    break
                if valu1eitem[-1]==valueitem[0]:
                    sal[key+"|"+item]=keyvalue+"".join(valueitem[1:])
                    delal.extend([key,item])
                    havecor="TRUE"
                    break
            #print(havecor)
            if havecor=="FASLE":
                sal[item]=lastsal[item]
                #print(sal)
    return sal

import lmdb
import msgpack

df = pd.read_csv("vcf.extr",sep="\t",header=None,keep_default_na=False,na_values=[],low_memory=False)
df.iloc[:,9:549] = df.iloc[:,9:549].applymap(str) ###modify过了不然有的是数字有的是字符

header = pd.read_csv("headern", sep="\t", keep_default_na=False, na_values=[], nrows=1)
df.columns=header.columns

with open("new.win.long", 'r') as f:
    lines = f.readlines()

win = pd.read_csv("new.win.long", sep="\t", keep_default_na=False, na_values=[],header=None)
win['name']=win[0]+":"+win[1].astype("str")+"-"+win[2].astype("str")
winn= win.groupby(5, as_index=False).first()
lines=winn.values.tolist()

print("OK")
extracted_ids = df['ID'].str.extract(r'(\d+).*?(\d+)')
### 4.T2T position file
#T2Tpos = pd.read_csv("/dssg/home/acct-clsmyf/clsmyf-user1/data/chr1.posCHM.csv",sep="\t")
T2Tpos = pd.read_csv("node.extr",sep="\t")
T2Tpos['node']=T2Tpos['node'].astype(str)
print(T2Tpos)
startdict=dict(zip(T2Tpos['node'], T2Tpos['start_position']))
enddict=dict(zip(T2Tpos['node'], T2Tpos['end_position']))
df['st1']=list(map(startdict.get, extracted_ids[0]))
df['en1']=list(map(startdict.get, extracted_ids[0]))
df['st2']=list(map(startdict.get, extracted_ids[1]))
df['en2']=list(map(startdict.get, extracted_ids[1]))

df['st']=df[['st1','en1','st2','en2']].min(axis=1)
df['en']=df[['st1','en1','st2','en2']].max(axis=1)
df = df.drop(['st1', 'en1','st2','en2'], axis=1)  # axis=1 表示操作列

import re
def process_line(line):
    env = lmdb.open(args.d2, map_size=20 * 1024 * 1024 * 1024,create=False,readonly=True, lock=False,max_readers=1000)  # 
    allsim={}
    cor=[]
    print(line)
    wstart=int(line[2])
    wend=int(line[3])
    CHR=line[1]
    dfsel=df[(df['st']<wend) & (df['en']>wstart)]
    dfall=dfsel
    if len(dfall)!=0:
        snw=min(dfsel['st'])
        enw=max(dfsel['en'])
    else:
        snw=wstart
        enw=wend
    snw1=min(wstart,snw)
    enw1=max(wend,enw)
    T2Tpath=T2Tpos[(T2Tpos['start_position'] <=enw ) & (T2Tpos['end_position'] >=snw1)]
    T2Tpath['merged'] =  T2Tpath['ori'] + T2Tpath['node'].astype(str)
    fa_series = T2Tpath.loc[T2Tpath['ori'] == '<', 'fa']
    fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
    T2Tpath.loc[T2Tpath['ori'] == '<', 'fa'] = fa_series
    T2Tseq=''.join(T2Tpath['fa'])
    allT2Tpath=''.join(T2Tpath['merged'])
    name=CHR+":"+str(snw1)+"-"+str(enw1)
    namewin=CHR+":"+str(wstart)+"-"+str(wend)
    dfseltry=df[(df['st']<enw1) & (df['en']>snw1)]
    if len(dfsel)!=len(dfseltry):
        print("NO!!!!!!!!!!!!!!!!!!!!!!!",line)
    if dfall.shape[0]==0:
        FA=args.r
        fa1 = pysam.FastaFile(FA)
        samplesim={}
        samplesim[('T2TCHM13')]=str(fa1.fetch(region=CHR + ":" + str(wstart) + "-" + str(wend)))
        allsim[namewin]=samplesim
        cor.append([namewin,namewin])
        return allsim,cor
    dfall.iloc[:,9:549] = dfall.iloc[:,9:549].applymap(str) ###motify过了不然有的是数字有的是字符
    dfa=dfall.iloc[:,9:549]

    groups=group_columns_by_values(dfa)
    dfv=dfall.iloc[:,:8]
    dfall['AT_content'] = dfall['INFO'].apply(lambda x: re.search(r'AT=([^;]+);', x).group(1) if re.search(r'AT=([^;]+);', x) else None)
    dfall = dfall.reset_index(drop=True)
    AT=dfall['AT_content']  ##新建一个字典去存储不同的位置的segment
    var_segment={}
    for varindex, at_content in dfall['AT_content'].items():
        pathdic={}
        for gtindex, path in enumerate(at_content.split(",")):
            pathdic[gtindex]=path # symbols = re.findall(r'[<>]',path) # numbers = re.findall(r'\d+', path) [symbols,numbers]
        pathdic['.']=">00000"
        var_segment[varindex]=pathdic
    samplesim={}
    for genotype in groups.keys():
        if genotype.count('.')/(len(genotype))>0.5:
            samplesim[tuple(groups[genotype])] = "NNNNNNNNNNNNNNNNNNNNN"
            continue
        if all(x == 0 for x in genotype):
            samplesim[tuple(groups[genotype])] = T2Tseq
        else:
            replacements={}
            replacementsdot={}
            for varindex, value in enumerate(genotype):
                if value=='0' or genotype==0 :
                    continue
                if value!=0 and value!=".":
                    replacements[var_segment[varindex][0]]=var_segment[varindex][int(value)]
                if value==".":
                    replacementsdot[var_segment[varindex][0]]=var_segment[varindex][value]
            if len(replacementsdot)!=0:
                all_numbers = [num for key in replacementsdot.keys() for num in re.findall(r'\d+', key)]
                if len(all_numbers)!=len(set(all_numbers)):
                    refi=refine_conti(replacementsdot)
                    replacementsdot={key: ">00000" for key in refi.values()}
            text=allT2Tpath
            for old_pattern, new_pattern in replacements.items():
                text = text.replace(old_pattern, new_pattern)
            for old_pattern, new_pattern in replacementsdot.items():
                text = text.replace(old_pattern, new_pattern)
            symbols = re.findall(r'[<>]',text)
            numbers = re.findall(r'\d+', text)
            end=pd.DataFrame()
            end.index=numbers
            with env.begin() as txn:
                end['fa'] = ["NNNNNNNNNNNNNNNNNNNNN" if key == "00000" else msgpack.unpackb(txn.get(key.encode('utf-8'))) for key in numbers] 
            #end['fa']=list(map(fasta_dict.get, numbers))  #end.isnull().values.any()
            end['ori']=symbols
            fa_series = end.loc[end['ori'] == '<', 'fa']
            fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
            end.loc[end['ori'] == '<', 'fa'] = fa_series
            samseq=''.join(end['fa'])
            samplesim[tuple(groups[genotype])] = samseq
    samplesim[('T2TCHM13')]=T2Tseq
    allsim[namewin]=samplesim
    cor.append([namewin,name])
    return allsim,cor



with mp.Pool(processes=int(args.p)) as pool:
    res=pool.map(process_line, lines)



with open("res.pkl", "wb") as f:  # 'wb' 表示二进制写入
    pickle.dump(res, f)


del df
#del fasta_dict
del T2Tpos
### 5.kmer database
#with open('/dssg/home/acct-clsmyf/clsmyf-user1/data/SD_lowerall829.pkl', 'rb') as f:
with open(args.d1+"inter_seg_kmer.pkl", 'rb') as f:
    loaded_set = pickle.load(f)



def process_item(item):
    keys=next(iter(item[0]))
    id=winn[winn['name']==keys][5].values[0]
    subwin=win[win[5]==id]
    # print(keys)
    simdf=[]
    for sample in item[0][keys]:
        samseq=item[0][keys][sample]
        samkmer=kmer_cal(Seq(samseq),31)
        if len(samkmer)!=len(set(samkmer)):
            print(key,"sam...NOHNONONONONO-----------------")
        samkmer=set(samkmer)
        for T2Twinn in list(subwin['name']):
            FA=args.r
            fa1 = pysam.FastaFile(FA)
            T2Tseq=str(fa1.fetch(region=T2Twinn))
            T2Tseqkmer=kmer_cal(Seq(T2Tseq),31)
            if len(T2Tseqkmer)!=len(set(T2Tseqkmer)):
                print(T2Twinn,key,"OHNONONONONO-----------------")
            T2Tseqkmer=set(T2Tseqkmer)
            if len(samkmer)!=0:
                T2TSDkmer=T2Tseqkmer & loaded_set ###对应纯粹window下的SDkmer
                sampleSDkmer=T2TSDkmer & samkmer ###sam与T2TCHM13一致的SDkmer
                samplecorkmer=T2Tseqkmer & samkmer ###sam与T2TCHM13一致的kmer
                samplenocorkmer=samkmer-T2Tseqkmer ###sam与T2TCHM13不一致的kmer
                samplenocorSDkmer=samplenocorkmer & loaded_set ###sam与T2TCHM13不一致的SDkmer
                if (len(T2Tseqkmer)==0):
                    sim1="."
                else:
                    sim1 = len(sampleSDkmer) / len(T2Tseqkmer)
                if (len(samplenocorkmer)==0):
                    sim2="."
                else:
                    sim2 = len(samplenocorSDkmer) / len(samplenocorkmer)
                sim3= len(samkmer & loaded_set) / len(samkmer)
            else:
                sim="."
                sim1="."
                sim2="."
                sim3="."
            if sample=="T2TCHM13":
                simdf.append(["T2TCHM13",sim1,sim2,sim3,item[1][0][1],T2Twinn])
            else:
                for sam in sample:
                    simdf.append([sam,sim1,sim2,sim3,item[1][0][1],T2Twinn])
    return simdf

res1=[]
for ress in res:
    if ress[1][0][0] in list(winn[winn[4]-winn[3]<=20000]['name']):
        res1.append(ress)

with mp.Pool(processes=int(args.p)) as pool:
    results1 = pool.map(process_item, res1)


with open("results1.pkl", "wb") as f:  # 'wb' 表示二进制写入
    pickle.dump(results1, f)


res2=[]
for ress in res:
    if ress[1][0][0] in list(winn[winn[4]-winn[3]>20000]['name']):
        res2.append(ress)

with mp.Pool(processes=int(args.p)) as pool:
    results2 = pool.map(process_item, res2)

with open("results2.pkl", "wb") as f:  # 'wb' 表示二进制写入
    pickle.dump(results2, f)


# for ress in res:
#     if  ress[1][0][0]=="chr1:674001-675000":
#         break



##sim1:样本中一致的SDkmer/(窗口对应的kmer数量)
##sim2:不一致的kmer中的SDkmer/(所有不一致的kmer)
##sim3:样本SDkmer/所有的SDkmer


allis=[]
for xx in results1:
    allis.extend(xx)


for xx in results2:
    allis.extend(xx)


data=pd.DataFrame(allis)
data.to_csv(args.o,sep='\t')
