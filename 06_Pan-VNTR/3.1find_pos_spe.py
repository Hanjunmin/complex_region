import pysam
import pandas as pd
import pickle
from collections import defaultdict
import glob
import argparse
import pysam
import sys
sys.path.append('./project/ENDVNTR')
from vntr_fun import *
from Bio.Seq import Seq
import warnings
warnings.filterwarnings('ignore')
from collections import Counter

parser = argparse.ArgumentParser(description='OTHER-SDcoordi')
parser.add_argument('--sam', help='sample')
args = parser.parse_args()


def generate_motif_pos(iterid):
    motiflen=len(cnmotif)
    st=stal-iterid*motiflen-20
    en=enal+iterid*motiflen+20
    if st==0:
        st=st+1
    if st<0:
        st=1
    if en>maxlength:
        en=maxlength
    T2TSeq = fa.fetch(region=subchr+":"+str(st)+"-"+str(en))
    ref=T2TSeq.upper()
    storelistn=[]
    motif=cnmotif
    resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
    expand_motif(resdict,motif, ref, interval,0.2,alllenlist)
    if alllenlist!=[]:
        cn=cal_cn1(alllenlist,len(motif),0.3)
        cn=list(cn)
        motiflist=[]
        for mo in cn[3]:
            motiflist.append(ref[mo[0]:(mo[1]+1)])
        cn.append(motiflist)
        storelistn.append(cn)   
    return storelistn,st

def refinecor(num):
    if num<1:
        num=1
    if num>maxlength:
        num=maxlength
    return num
    



# seq="GGGTTGGTCCTGCCCCTCAGAGGCAGAAAGGGCAGAAGGGGCAGAAGAG*******GTGGAGAAGGTGGAAGCGGGGCAGAAGAG*******GTGGAGAAGGTGGAAGGGAGGCAGGAGAGGCAGGCACATTCAGTGTAAC"
# seq=seq.replace("*", "")
# motif="GGGGCAGAAGAGGTGGAGAAGGTGGAAG"

# sam="C001-CHA-E01-Mat"        
# subchr="C001-CHA-E01_Mat_chr10"
# s=124354061     
# e=124354077



sam=args.sam   #"NA19650_hap2"
print(sam)
#sam="GRCh38"
FA="./allfa/"+sam+".fa"
fa = pysam.FastaFile(FA)
with open(sam+".out", 'r') as f:
    lines = f.readlines()

id=0
x=[]
alldflis=[]
for line in lines[1:]:
    line_list = line.strip().split('\t')
    subchr=line_list[13]
    subchrs=int(line_list[14])
    maxlength=fa.get_reference_length(subchr)
    posa=[]
    for subl in line_list[12].split("@"):
        posa.append((int(subl.split("_")[0]),int(subl.split("_")[1])))
    have="FALSE"
    for possam in posa:
        s=possam[0]+subchrs
        e=possam[1]+subchrs
        seq=line_list[7].replace("N****", "")
        seq=seq.replace("*****N", "")
        seq=seq.replace("*", "")
        motif=line_list[4]
        first=refinecor(s-len(seq))
        second=refinecor(e+len(seq))
        T2TSeq = fa.fetch(region=subchr+":"+str(first)+"-"+str(second))
        T2TSeq=T2TSeq.upper()
        if T2TSeq.find(seq) != -1:
            seq=seq
            have="TRUE"
            break
        # elif T2TSeq.find(str(Seq(seq).reverse_complement()))!=-1:
        #     seq=str(Seq(seq).reverse_complement())
    if have=="FALSE":
        id=id+1
        x.append(line)
        print(sam)
        print(line)
        print("NO")
        break
    stal=first+T2TSeq.find(seq)
    enal=stal+len(seq)
    T2TSeq = fa.fetch(region=subchr+":"+str(stal)+"-"+str(enal))
    cnmax=0
    cnmotif="a"
    ref=T2TSeq
    resdict,interval, motiflist, alllenlist1 = defaultdict(list),(0, len(ref)), [motif], []
    expand_motif(resdict,motif, ref, interval,0.2,alllenlist1)
    samseqrev=str(Seq(T2TSeq).reverse_complement())
    ref=samseqrev
    resdict,interval, motiflist, alllenlist2 = defaultdict(list),(0, len(ref)), [motif], []
    expand_motif(resdict,motif, ref, interval,0.2,alllenlist2)
    if alllenlist1!=[]:
        enall=cal_cn1(alllenlist1,len(motif),0.3)
        cn1=enall[0]
    if alllenlist2!=[]:
        enall=cal_cn1(alllenlist2,len(motif),0.3)
        cn2=enall[0]
    if alllenlist1!=[]:
        if cn1>cnmax:
            cnmax=cn1
            cnmotif=motif
    if alllenlist2!=[]:
        if cn2>cnmax:
            cnmax=cn2
            cnmotif=str(Seq(motif).reverse_complement())
    if cnmax==0:
        poslis=[line_list[0],line_list[1],line_list[2],0,motif]
        alldflis.append(poslis)
        continue
    else:
        initcn=cnmax  
        nextcn=initcn+1
        initid=1
    while nextcn>initcn:
        [storelist,st]=generate_motif_pos(initid)
        out1=process_vamos(storelist)
        nextcn=out1[0][0]
        initid=initid+2
        if nextcn>initcn:
            initcn=nextcn
    molis=[out1[0][4]] if not isinstance(out1[0][4], list) else set(out1[0][4])

    positions = []
    for start, end in out1[0][3]:
        positions.extend(range(start, end + 1))  # 展开区间为单个位置
    position_counts = Counter(positions)
    count_one_positions = [pos for pos, count in position_counts.items() if count == 1]
    simindex=len(count_one_positions)/((out1[0][2]-out1[0][1])+1)
    poslis=[line_list[0],line_list[1],line_list[2],subchr,nextcn,sam,len(list(molis)[0]),st+out1[0][1],st+out1[0][2],",".join(molis),simindex]
    alldflis.append(poslis)

        



####str(Seq(motif).reverse_complement())
pd.DataFrame(alldflis).to_csv(sam+".END", index=False,sep="\t")  # 不保存索引列


