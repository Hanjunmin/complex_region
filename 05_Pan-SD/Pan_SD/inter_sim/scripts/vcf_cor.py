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
parser.add_argument('--seg', help='chromsome segment position')
parser.add_argument('--vcf', help='chromsome vcf file')
parser.add_argument('--h', help='vcf file header')
parser.add_argument('--o', help='output file')
args = parser.parse_args()


print("OK")
### 2.file title
header = pd.read_csv(args.h, sep="\t", keep_default_na=False, na_values=[], nrows=1)
print("OK")
### 3.vcf file
#df = pd.read_csv("chr1.vcf",sep="\t",header=None,keep_default_na=False,na_values=[],low_memory=False)
df = pd.read_csv(args.vcf,sep="\t",header=None,keep_default_na=False,na_values=[],low_memory=False)
print("OK")
df.columns=header.columns
extracted_ids = df['ID'].str.extract(r'(\d+).*?(\d+)')
### 4.T2T position file
#T2Tpos = pd.read_csv("/dssg/home/acct-clsmyf/clsmyf-user1/data/pos/chrY.posCHM.csv",sep="\t")
T2Tpos = pd.read_csv(args.seg,sep="\t")
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
df[['#CHROM','st','en']].to_csv(args.o, sep='\t', index=False, header=False)
#simdf_transposed.to_csv("chr1.sim",sep='\t')



