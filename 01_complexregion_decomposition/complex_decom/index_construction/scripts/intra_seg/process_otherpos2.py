import pandas as pd
from collections import defaultdict
import lmdb
import re
import msgpack
import warnings# 忽略所有警告
import argparse
import json
import multiprocessing as mp
import pysam
from Bio.Seq import Seq
import edlib
import pandas as pd
import subprocess
import sys
from pyfaidx import Fasta
import warnings
warnings.filterwarnings('ignore')
# sys.path.append('/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/')
# from vnre_fun_moti import *

parser = argparse.ArgumentParser(description='Pan-VNTR-coordi')
parser.add_argument('--s', help='sample')
parser.add_argument('--o', help='output')
parser.add_argument('--w', help='the fold of W path')
parser.add_argument('--json', help='the file of json')
parser.add_argument('--script', help='script of VNTR calculation')
parser.add_argument('--fapath', help='Fasta path')
args = parser.parse_args()

sys.path.append(args.script)  ##'/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/'
from vnre_fun_moti import *

# INPUT=args.i
# INPUT="INSOTHER.al"
# data=pd.read_csv(INPUT,sep="\t",header=None)
ALSAM=pd.read_csv("INSOTHER.al.sam", sep='\t')
ALSAM=list(set(ALSAM['7']))
findpos_cleaned=pd.read_csv("INSOTHER.al.clean", sep='\t')


# env = lmdb.open('/share/home/zhanglab/user/maoyafei/project/pansd/pandb/CHM13-APGp1-HPRCp1-HGSVCp3_MC_db_msgpack', create=False,readonly=True, lock=False)  # 10 GB
with open(args.json, 'rt', encoding='utf-8') as gz_file:
    fasta_dict = json.load(gz_file)




all_positions = []
for pos_str in findpos_cleaned['13']:
    positions = [int(x) for x in pos_str.split(',')]
    all_positions.extend(positions)

allpos=list(set(all_positions)) ###相关的segment
allpos=[str(pos) for pos in allpos]

sam=args.s


# sam="C001-CHA-E01-Mat"

samdict={}
allWfold=args.w
Wfile=allWfold+sam+".W" #"/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/Wpath/"
with open(Wfile, 'r') as f:
    lines = f.readlines()

local_dict = defaultdict(list)  # 使用局部字典减少锁竞争
for line in lines:
    line_list = line.strip().split('\t')
    samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
    local_dict[samname] = line_list[6]
# subdf=data[data['sam']==sam]

keys = ALSAM
local_dict = {k: local_dict[k] for k in keys if k in local_dict}
allist=[]
dic_pos_s={}  ####这里生成一个字典 key为不同的segment value为对应的pos
dic_pos_e={}  ####这里生成一个字典 key为不同的segment value为对应的pos

dic_path={}  ####这里生成一个字典 key为不同的segment value为对应的path name
for key in local_dict.keys():
    print(key)
    value=local_dict[key]
    symbols = re.findall(r'[<>]',value)
    numbers = re.findall(r'\d+', value)
    end = pd.DataFrame()
    end['ori']=symbols
    end['node']=numbers
    # with env.begin(write=False) as txn:
    #     cursor = txn.cursor()
    #     values = cursor.getmulti([key.encode('utf-8') for key in numbers])
    #     end['fa'] = [msgpack.unpackb(value[1]) for value in values if value[1] is not None]
    end['fa']  = pd.DataFrame(list(map(fasta_dict.get, numbers)))
    end['base_count'] = end['fa'].apply(lambda x: len(x))
    end['pos']=end['base_count'].cumsum()
    end['pos_s']=end['pos']-end['base_count']
    filtered_end = end[end['node'].isin(allpos)]
    filtered_end['path']=key
    subdic=filtered_end.set_index('node')['path'].to_dict() ###每个样本都生成相关的segment对应的位置信息
    dic_path.update(subdic)
    allist.append(filtered_end)
    subdic=filtered_end.set_index('node')['pos_s'].to_dict() ###每个样本都生成相关的segment对应的位置信息
    dic_pos_s[key]=subdic
    subdic=filtered_end.set_index('node')['pos'].to_dict() ###每个样本都生成相关的segment对应的位置信息
    dic_pos_e[key]=subdic

df=pd.concat(allist)
segment_positions=list(df['node'])
def check_positions_in_segments(pos_str):
    positions = set(str(x) for x in pos_str.split(','))
    corres_path=list(set([dic_path.get(pos) for pos in positions]))
    issub=positions.issubset(segment_positions)
    if len(corres_path)==1:
        return [issub,corres_path[0]]
    else:
        return [issub, None]


# 应用到DataFrame
res_series=findpos_cleaned['13'].apply(check_positions_in_segments) ##判断要检查的VNTR对应的segment是否存在于样本中
res_df = pd.DataFrame.from_records(res_series.tolist(), columns=['in_segments', 'corres_path'])
findpos_cleaned=pd.concat([findpos_cleaned, res_df], axis=1)

examin=findpos_cleaned[findpos_cleaned['in_segments']==True]
examin['path_part1'] = examin['corres_path'].str.split('@').str[-3]
examin['path_part2'] = examin['corres_path'].str.split('@').str[-2]
examin['chr']=examin['path_part1']+"@"+examin['path_part2']


###开始迭代计算拷贝数



import pandas as pd
import pickle
from collections import defaultdict
import glob
import argparse
import pysam
import sys
from collections import Counter
from Bio.Seq import Seq
import warnings
warnings.filterwarnings('ignore')


def generate_motif_pos(iterid,trans,molis):
    motiflen=len(row['5'])
    st=stal-iterid*motiflen-20
    en=enal+iterid*motiflen+20
    if st==0:
        st=st+1
    if st<0:
        st=1
    if en>maxlength:
        en=maxlength
    T2TSeq = fa.fetch(region=row['chr']+":"+str(st)+"-"+str(en))
    ref=T2TSeq.upper()
    if trans==2: ##chan2
        ref=str(Seq(T2TSeq).reverse_complement())
    storelistn=[]
    motif=cnmotif
    resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
    expand_motif(resdict,motif, ref, interval,0.2,alllenlist,molis)
    if alllenlist!=[]:
        cn=cal_cn1(alllenlist,len(motif),0.3)
        cn=list(cn)
        motiflist=[]
        for mo in cn[3]:
            motiflist.append(ref[mo[0]:(mo[1]+1)])
        cn.append(motiflist)
        storelistn.append(cn)   
    return storelistn,st,en-st ##chan1



# # initcn="/share/home/zhanglab/user/maoyafei/project/ENDVNTR/STRVNTR130/chr1.STR_VNTR.out" ##获取对应的motif
# initcn=args.init
# file_name=args.i
# #file_name="chr13_72othersampos.pkl"
# #initcn="chr13.other.all"
# T2Tpos = pd.read_csv(initcn,sep="\t",header=None)
# T2Tpos['inde'] = T2Tpos[0] + "_" + T2Tpos[1].astype(str) + "_" + T2Tpos[2].astype(str)
# motifdic=T2Tpos.set_index('inde')[9].str.split(",").to_dict()

###前面还是要加一步先call CNV 得到很多种motif 最后再按照相似性计算


alldflis=[]
print(sam)
FA=args.fapath+sam+".fa" #"/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/FA/"
fa = pysam.FastaFile(FA)
for index,row in examin.iterrows():
    maxlength=fa.get_reference_length(row['chr'])
    corpath=row['corres_path']
    motiflen=len(row['5'])
    seglis=row['13'].split(",") + [str(row['12'])]
    a=min(pd.Series(seglis).map(dic_pos_s[corpath]))
    b=max(pd.Series(seglis).map(dic_pos_e[corpath]))
    stal=a
    if stal==0:
        stal=1
    enal=b
    ####确定初始motif和对应的拷贝数
    T2TSeq = fa.fetch(region=row['chr']+":"+str(stal)+"-"+str(enal))
    T2TSeq=T2TSeq.upper()
    ini1=0       
    cnn1="a"   
    for cnmotif in [row['5']]: ###相当于找到最concensus的motif序列然后进行注释
        # print(cnmotif)
        initcn=0  
        outmax="a"
        nextcn=initcn+1
        initid=1
        lastcn=0
        conti=0
        iternum=0
        while nextcn>=initcn & conti<3:
            iternum=iternum+1
            if iternum>10:
                break
            [storelist,st,lena]=generate_motif_pos(initid,1,[row['5']]) ##chan3
            if len(storelist)==0:
                initid=initid+2
                continue
            out1=process_vamos(storelist)
            out1.append(lena)
            out1.append(st)
            # print(out1)
            nextcn=out1[0][0]
            initid=initid+2
            if nextcn==lastcn  or nextcn<lastcn:
                conti=conti+1
            if nextcn>lastcn:
                cncopy=out1
            if conti>=3:
                outmax=cncopy
                outcn=outmax[0][0]
                if outcn>ini1:
                    ini1=lastcn
                    cnn1=outmax
                break
            lastcn=nextcn
    # print(ini1)
    ini2=0  
    cnn2="a"        
    for cnmotif in [row['5']]: ###相当于找到最concensus的motif序列然后进行注释
        if ini1>10:
            break
        # print(cnmotif)
        initcn=0  
        outmax="a"
        nextcn=initcn+1
        initid=1
        lastcn=0
        conti=0
        iternum=0
        while nextcn>=initcn & conti<3:
            iternum=iternum+1
            if iternum>10:
                break
            [storelist,st,lena]=generate_motif_pos(initid,2,[row['5']]) ##chan3
            if len(storelist)==0:
                initid=initid+2
                continue
            out1=process_vamos(storelist)
            out1.append(lena)
            out1.append(st)
            nextcn=out1[0][0]
            initid=initid+2
            if nextcn==lastcn  or nextcn<lastcn:
                conti=conti+1
            if nextcn>lastcn:
                cncopy=out1
            if conti>=3:
                outmax=cncopy
                outcn=outmax[0][0]
                if outcn>ini2:
                    ini2=lastcn
                    cnn2=outmax
                break
            lastcn=nextcn
    if ini2>ini1:
        outmax=cnn2
        trans=2
    else:
        outmax=cnn1
        trans=1
    if ini2==0 and ini1==0:
        poslis=[row[0],0,sam,row['chr'],a,b,row['5'],-1,row['1'],row['2'],row['3'],row['5']]
    else:
        molis=[outmax[0][4]] if not isinstance(outmax[0][4], list) else set(outmax[0][4])
        if outmax[0][2]-outmax[0][1]!=0:
            positions = []
            for start, end in outmax[0][3]:
                positions.extend(range(start, end + 1))  # 展开区间为单个位置
            position_counts = Counter(positions)
            count_one_positions = [pos for pos, count in position_counts.items() if count == 1]
            simindex=len(count_one_positions)/((outmax[0][2]-outmax[0][1])+1)
            #simindex=(sum(end - start for (start, end) in outmax[0][3]))/(outmax[0][2]-outmax[0][1])
        else:
            simindex=-1
        if trans==2:
            poslis=[row[0],outmax[0][0],sam,row['chr'],outmax[-1]+outmax[-2]-outmax[0][2],outmax[-1]+outmax[-2]-outmax[0][1],",".join(molis),simindex,row['1'],row['2'],row['3'],row['5']]
        else:
            poslis=[row[0],outmax[0][0],sam,row['chr'],outmax[-1]+outmax[0][1],outmax[-1]+outmax[0][2],",".join(molis),simindex,row['1'],row['2'],row['3'],row['5']]
    alldflis.append(poslis)


alldflisdf=pd.DataFrame(alldflis)
alldflisdf.to_csv(args.o, index=True, sep="\t")

