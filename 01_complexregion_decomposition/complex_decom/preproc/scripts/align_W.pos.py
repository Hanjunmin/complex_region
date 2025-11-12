import lmdb
import msgpack
from collections import defaultdict
import gzip
import json
import os
import pandas as pd
import  pysam
import edlib
from Bio.Seq import Seq
import re
import numpy as np
import multiprocessing as mp
import math
import argparse
import concurrent.futures
from itertools import chain
import subprocess
import psutil
from datetime import datetime
import pickle


from multiprocessing import Pool

def process_chromosome(chromsome):
    print(f"Processing chromosome: {chromsome}")
    subchr = pair[pair['REF_CHR'] == chromsome]
    T2Tchmpos = "/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.posCHM.txt"
    T2Tpos = pd.read_csv(T2Tchmpos, sep="\t", header=None)
    T2Tpos.columns = ["CHR", "start_position", "end_position", "node", "orient"]
    T2Tpos['node'] = T2Tpos['node'].astype(str)
    node_to_start_dict = T2Tpos.set_index('node')['start_position'].to_dict()
    results = []
    for index, row in subchr.iterrows():
        refkey = row['REF_path']
        querykey = row['SAM_path']
        quennode = queryaldic[querykey]
        refnode = refaldic[refkey]
        internode = set(quennode) & set(refnode)
        interpos = list(map(lambda x: node_to_start_dict.get(x), internode))
        Ss = min(interpos)
        Es = max(interpos)
        Es=list(T2Tpos[T2Tpos['start_position'] == Es]['end_position'])[0]
        results.append((index, Ss, Es))
    return results




parser = argparse.ArgumentParser(description='Pan-SD-idencordinate')
parser.add_argument('--sam', help='processid')
parser.add_argument('--p', help='W path')
parser.add_argument('--r', help='reference')

args = parser.parse_args()

sam=args.sam   #"HG00512_hap1"
allWfold=args.p
queryW=allWfold+sam+".W"
ref=args.r
refW=allWfold+ref+".W"



'''querynode path'''
with open(queryW, 'r') as f:
    lines = f.readlines()

queryaldic = {}
for line in lines:
	line_list = line.strip().split('\t')
	samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
	nodelist = (line_list[6]).replace('<', ' ').replace('>', ' ').split()
	queryaldic[samname] = nodelist

'''refnode path'''
with open(refW, 'r') as f:
    lines = f.readlines()

refaldic = {}
for line in lines:
	line_list = line.strip().split('\t')
	samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
	nodelist = (line_list[6]).replace('<', ' ').replace('>', ' ').split()
	refaldic[samname] = nodelist






pair=pd.read_csv(sam+"_REFcor.list",sep="\t",header=None)
pair.columns=['REF_path','SAM_path']
pair[['REFSAM', 'x1', 'REF_CHR', 'REF_s', 'REF_e']] = pair['REF_path'].str.split('@', expand=True)
pair[['SAM', 'x2', 'SAM_CHR', 'SAM_s', 'SAM_e']] = pair['SAM_path'].str.split('@', expand=True)
pair['QUE_REFS']=0
pair['QUE_REFE']=0


chromosomes = set(list(pair['REF_CHR']))
with Pool(processes=4) as pool: 
    all_results = pool.map(process_chromosome, chromosomes)

all_results_flat = [result for sublist in all_results for result in sublist]
for index, Ss, Es in all_results_flat:
    pair.loc[index, 'QUE_REFS'] = Ss
    pair.loc[index, 'QUE_REFE'] = Es


END=pair[['REFSAM','REF_CHR', 'REF_s', 'REF_e','SAM','SAM_CHR', 'SAM_s', 'SAM_e','QUE_REFS','QUE_REFE','REF_path','SAM_path']]

END.to_csv(sam+"_REFcor.pos", sep="\t", index=False, header=True)
