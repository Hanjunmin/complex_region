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
parser.add_argument('--sam', help='processid')
parser.add_argument('--p', help='Wpath')
parser.add_argument('--r', help='reference name')

args = parser.parse_args()


samin=args.sam
# samin="HG00512_hap1"
# ref="CHM13v2"
ref=args.r
# allWfold="./project/APG/github_final/01_1p36/DATA/graph/Wpath/"
allWfold=args.p
queryW=allWfold+samin+".W"
refW=allWfold+ref+".W"

'''querynode path'''
with open(queryW, 'r') as f:
    lines = f.readlines()

queryaldic = {}
for line in lines:
	line_list = line.strip().split('\t')
	samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
	nodelist = (line_list[6]).replace('<', ' ').replace('>', ' ').split()
	queryaldic[samname] = set(nodelist)

'''refnode path'''
with open(refW, 'r') as f:
    lines = f.readlines()

refaldic = {}
for line in lines:
	line_list = line.strip().split('\t')
	samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
	nodelist = (line_list[6]).replace('<', ' ').replace('>', ' ').split()
	refaldic[samname] = set(nodelist)

Wpair=[]
for key,value in refaldic.items():
	print(key)
	for key1,value1 in list(queryaldic.items()):
		if len(value & value1)!=0:
			Wpair.append([key,key1])
			del queryaldic[key1]  # 删除 key1
			# print("OK",key,key1)

Wpairdf=pd.DataFrame(Wpair)
weno=queryaldic.keys()
Wpairdf.to_csv(samin+"_REFcor.list", sep="\t", index=False, header=False)
pd.DataFrame(weno).to_csv(samin+"_NOalign.list", sep="\t", index=False, header=False)
