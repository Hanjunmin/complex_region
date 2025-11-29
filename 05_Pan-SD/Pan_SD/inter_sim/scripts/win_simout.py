import pandas as pd
import numpy as np
import glob
import os
import pandas as pd
import re
from Bio.Seq import Seq
import pickle
from kmerhash import kmerhasher, hashes2seq
import json
import multiprocessing as mp
import pysam
import argparse
import lmdb
import msgpack


parser = argparse.ArgumentParser(description='SD-Sim')
parser.add_argument('--p', help='running path')
parser.add_argument('--r', help='reference file')
parser.add_argument('--d', help='database of inter_seg homo kmer')
args = parser.parse_args()




combined_df = pd.DataFrame()
#base_dir="/home/jmhan/project/APG/github_final/PanSD/output/"
base_dir=args.p
df = pd.read_csv(base_dir+"01_intersim/sim.out", sep='\t',index_col=0)
win=df['win']
df=df.apply(pd.to_numeric, errors='coerce')
dfmax=pd.DataFrame(df.iloc[:, :-1].max(axis=1, skipna=True))
dfmax['wind']=win
combined_df = pd.concat([combined_df, dfmax], axis=0)
combined_df.to_csv('simmax.al', sep='\t', index=True)




alw=pd.read_csv(base_dir+"00_preproc/WGS.window.bed",sep="\t",header=None)
bef=pd.read_csv("simmax.al",sep="\t",header=0)
alw['name']=alw[0]+":"+alw[1].astype("str")+"-"+alw[2].astype("str")
noww=alw[~alw['name'].isin(bef.iloc[:, 0])]
def kmer_cal(seq,size):
    kmers = [str(seq[i:i+size]) for i in range(len(seq) - size + 1)]
    kmers = list(filter(lambda kmer: 'N' not in kmer, kmers))
    rkmers = [str(seq[i:i+size].reverse_complement()) for i in range(len(seq) - size + 1)]
    rkmers = list(filter(lambda kmer: 'N' not in kmer, rkmers))
    result = [min(kmer1, kmer2, key=lambda x: kmerhasher(x, size)) for kmer1, kmer2 in zip(kmers, rkmers)]
    result= {x.lower() for x in result}
    return(set(result))


with open(args.d, 'rb') as f:
        loaded_set = pickle.load(f)

print("OK")
#fasta_dict['00000']="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
### 2.file title
def process_line(line):
    print(line)
    allsim={}
    cor=[]
    wstart=int(line[1])
    wend=int(line[2])
    FA=args.r
    fa1 = pysam.FastaFile(FA)
    samplesim={}
    namewin=line[0] + ":" + str(wstart) + "-" + str(wend)
    samseq=str(fa1.fetch(region=line[0] + ":" + str(wstart) + "-" + str(wend)))
    samkmer=kmer_cal(Seq(samseq),31)
    sim = len(samkmer & loaded_set) / len(samkmer)
    allsim[namewin]=sim
    return allsim

nowlis=noww.values.tolist()
with mp.Pool(processes=50) as pool:
    res=pool.map(process_line, nowlis)

merged_dict = {}
for d in res:
    merged_dict.update(d)

simdf = pd.DataFrame({
    'win': merged_dict.keys(),
    'T2TCHM13init': merged_dict.values()
})


en1=bef.iloc[:,1:3]
en1.columns=['maxsim','win']
simdf.columns=['win','maxsim']

end=pd.concat([simdf,en1],axis=0)
end.to_csv('simmax.al.add', index=False, float_format='%.10f',sep="\t")