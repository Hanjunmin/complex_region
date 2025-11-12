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
parser = argparse.ArgumentParser(description='Pan-SD-idencordinate')
parser.add_argument('--pr', help='path region')
args = parser.parse_args()



subpath=args.pr
# subpath="/home/jmhan/COMPLEXFEATURE/09_bubble_only/01_GFA/chr1/chr1:1-2000000.path"
sdpa=pd.read_csv("Sregional.node",sep="\t",header=None)


aldicmulti = defaultdict(list)
aldicmultipath=defaultdict(list)
with open(args.pr, 'r') as f:
    lines = f.readlines()


for line in lines:
    line_list = line.strip().split('\t')
    samname=line_list[1]+"|"+line_list[2]+"|"+line_list[3]+"|"+line_list[4]+"|"+line_list[5]
    # print(samname)
    nodelist = (line_list[6]).replace('>', ' ').replace(',', ' ').replace('<', ' ').split()
    nodelistdict = {value: index for index, value in enumerate(nodelist)}
    aldicmulti[samname] = nodelistdict
    aldicmultipath[samname] = line_list[6]


chrpathsam=aldicmulti.keys()

def process_item(item):
    print(item)
    alist=[]
    t2tmotiflist=[]
    i=0
    vntrid=item
    i=i+1
    storep=[]
    # print(vntrid)
    vntridcopy_number={}
    result_dict = {}
    allnodesetsub=list(sdpa[sdpa[2]==item][1])
    allnodesetsub= set([str(i) for i in allnodesetsub])
    numbers = allnodesetsub
    allnodesetsub= [str(i) for i in allnodesetsub]
    for sam in chrpathsam:
        # print(sam)
        pa=aldicmultipath[sam]
        samno=aldicmulti[sam]
        end=pd.DataFrame(numbers)
        end.index=numbers
        end['node']=end.index
        end['sampos'] = end['node'].apply(lambda x: samno.get(x))
        no=end.dropna(subset=['sampos'])
        if len(no)==0:
            continue
        fro=min(no['sampos'])
        to=max(no['sampos'])
        ndori=re.findall(r'[><]',pa)
        nd = re.findall(r'\d+', pa)
        sanoden=nd[int(fro):(int(to)+1)]
        nd_ori = pd.DataFrame({'nd': sanoden, 'ndori': ndori[int(fro):(int(to)+1)]})
        nd_ori['ndori'] = nd_ori['ndori'] .replace({'<': '-,', '>': '+,'})
        nd_ori['mer']=nd_ori['nd']+nd_ori['ndori']
        path=''.join(nd_ori['mer'].astype(str))
        if len(set(sanoden)-set(allnodesetsub))==0:
            result_dict[vntrid+"@"+sam]=path
        else:
            with open("errlog", "a") as log_file:
                log_file.write(vntrid+"@"+sam+"PATH ERROR!" + "\n")
                # print("PATH ERROR!")  
    return result_dict

with mp.Pool(processes=10) as pool:
    print("E")
    results = pool.map(process_item, set(sdpa[2]))


storep=[]
for res in results:
    for key,value in res.items():
        storep.append(['P',key,value,"*"])

pathdf=pd.DataFrame(storep)
pathdf.to_csv("path.txt", sep="\t", index=False, header=False)