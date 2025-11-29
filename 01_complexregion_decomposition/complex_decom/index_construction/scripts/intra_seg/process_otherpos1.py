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
# sys.path.append('/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/')
# from vntr_fun import *


warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Pan-VNTR-coordi')
parser.add_argument('--i', help='input file')
parser.add_argument('--script', help='script of VNTR calculation')
parser.add_argument('--r', help='reference')
parser.add_argument('--json', help='the file of json')

args = parser.parse_args()
sys.path.append(args.script)  ##'/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/'
from vntr_fun import *


INPUT=args.i
INPUT="INSOTHER.al"
data=pd.read_csv(INPUT,sep="\t",header=None)
split_columns = data[7].str.split('@', expand=True)
data['sam']=split_columns[0]
data['chr']=split_columns[2]+"@"+split_columns[3]


# env = lmdb.open('/share/home/zhanglab/user/maoyafei/project/pansd/pandb/CHM13-APGp1-HPRCp1-HGSVCp3_MC_db_msgpack', create=False,readonly=True, lock=False)  # 10 GB
with open(args.json, 'rt', encoding='utf-8') as gz_file:
    fasta_dict = json.load(gz_file)



###先看一下T2TCHM13对应是否有CNV CN>2
##左右扩ref，判断是否有额外的拷贝
fafile=args.r
fa = pysam.FastaFile(fafile)
fasta = Fasta(fafile)
chromosome_max_positions = {chrom: len(sequence) for chrom, sequence in fasta.items()}

alldf=[]
i=1
with open("INSOTHER.al", 'r') as file:
    for line in file:
        i=i+1
        print(i)
        lis=line.strip().split('\t')
        motif=lis[5]
        max_value=chromosome_max_positions[lis[1]]
        st=(int(float(lis[2]))-2*len(motif))
        if lis[3]=="INS":
            lis[3]=lis[2]
        en=(int(float(lis[3]))+2*len(motif))
        name =lis[1] + ":" + str(st) + "-" +  str(en)
        if st<=0:
            name =lis[1] + ":" + str(1) + "-" + str(en)
        if en>max_value:
            name =lis[1] + ":" + str(st) + "-" + str(max_value)
        Seq1lef = (fa.fetch(region=name)).upper() ##refseqall=lis[7].split("*******")[0]+lis[7].split("*******")[2]
        refseqallrev=str(Seq(Seq1lef).reverse_complement())  ##refseqallrev=str(Seq(refseqall).reverse_complement())
        ref=refseqallrev
        resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
        expand_motif(resdict,motif, ref, interval,0.2,alllenlist)
        if alllenlist!=[]:
            motiflist1=[]
            enall1=cal_cn1(alllenlist,len(motif),0.3)
        else:
            enall1=[0]
        ref=Seq1lef
        resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
        expand_motif(resdict,motif, ref, interval,0.2,alllenlist)
        if alllenlist!=[]:
            motiflist2=[]
            enall2=cal_cn1(alllenlist,len(motif),0.3)
        else:
            enall2=[0]
        cnnow=max(enall1[0],enall2[0])
        if cnnow>2:
            continue
        else:
            alldf.append(line.strip().split('\t'))


END=pd.DataFrame(alldf)
findpos=END.iloc[:,[1,2,3,5,12,13]]
findpos_cleaned = findpos.drop_duplicates()
END[7].to_csv("INSOTHER.al.sam", sep='\t', index=False)
findpos_cleaned.to_csv("INSOTHER.al.clean", sep='\t', index=False)

samlis = list(set(data['sam']))
df = pd.DataFrame({'sample': samlis})
df.to_csv("samlis", sep='\t', index=False)
