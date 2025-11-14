import lmdb
import msgpack
import os
import re
import multiprocessing
import argparse
import csv
import re
import numpy as np
import pandas as pd
import json
import gzip
import hdbscan
from sklearn.cluster import DBSCAN
import lmdb
import msgpack
from collections import defaultdict
from io import StringIO
from multiprocessing import Pool
from collections import defaultdict
import subprocess

###################################

def process_single_line(line):
	line_list = line.strip().split('\t')
	base_samname = f"{line_list[1]}_{line_list[2]}_{line_list[3]}_{line_list[4]}_{line_list[5]}"
	letters_to_find = [part for part in re.split(r'[<>]', line_list[6]) if part]
	number_positions = {number: index + 1 for index, number in enumerate(letters_to_find)}
	endfa = pd.DataFrame()
	orilis=re.findall(r'[<>]', line_list[6])
	endfa['node']=letters_to_find 
	endfa['ori']=orilis
	with env.begin() as txn:
		endfa['fa'] = [msgpack.unpackb(txn.get(key.encode('utf-8'))) for key in letters_to_find] #endfa['fa']=list(map(fasta_dict.get, letters_to_find))  #end.isnull().values.any()
	endfa['start_position'] = endfa['fa'].apply(len).cumsum() - endfa['fa'].apply(len) + 1
	endfa['end_position'] = endfa['start_position'] + endfa['fa'].str.len()-1
	endfa['chr_position'] = base_samname +"_"+ endfa['start_position'].astype(str)
	endfa['chr']=base_samname
	chunk_size=50000
	returnlis=[]
	for start_idx in range(0, len(endfa), chunk_size):
		end_idx = min(start_idx + chunk_size, len(endfa))
		chunk = endfa.iloc[start_idx:end_idx]
		ls=min(chunk['start_position'])
		lr=max(chunk['end_position'])
		chunk_samname = f"{base_samname}_rows_{ls}_to_{lr}"
		df_dict = defaultdict(list)
		if len(set(chunk['node']))==len(chunk):
			for node, start, end, chr_pos, chr_val,chr_ori in zip(chunk['node'], chunk['start_position'], chunk['end_position'], chunk['chr_position'], chunk['chr'],chunk['ori']):
				df_dict[node]=[[start, end, chr_pos, chr_val,chr_ori]]
		else:
			for node, start, end, chr_pos, chr_val,chr_ori in zip(chunk['node'], chunk['start_position'], chunk['end_position'], chunk['chr_position'], chunk['chr'],chunk['ori']):
				df_dict[node].append([start, end, chr_pos, chr_val,chr_ori])
		chunk = chunk.copy()  # Make an explicit copy
		chunk.loc[:, 'nindex'] = (chunk.index.astype('int') + 1)
		chunk.index=chunk['node']
		key_value_dict = dict(chunk[['node', 'nindex']].values)
		letters_to_findsub=list(chunk['node'])
		orilischunk=list(chunk['ori'])
		returnlis.append([chunk_samname.encode('utf-8'), msgpack.packb(df_dict),msgpack.packb(orilischunk),msgpack.packb(key_value_dict),chunk_samname,letters_to_findsub])
	return base_samname,returnlis
# return samname.encode('utf-8'), msgpack.packb(df_dict),msgpack.packb(orilis),msgpack.packb(number_positions),samname,letters_to_find

def process_lines(lines,envfinaldic,envpos,envori,batch_size=20, num_workers=20):
	cachedfdict = {}
	cacheori = {}
	cachepos = {}
	cacheseg = {}
	i = 0
	pool = multiprocessing.Pool(num_workers)
	for samnamb, lis in pool.imap(process_single_line, lines):
		print(samnamb)
		for samnameenc,packed_dfdict,packed_ori,packed_pos,samename,seg in lis:
			print(samnameenc)
			cachedfdict[samnameenc] = packed_dfdict
			cacheori[samnameenc] = packed_ori
			cachepos[samnameenc] = packed_pos
			cacheseg[samename]=seg
			i += 1
			if (i % batch_size == 0):
				with envfinaldic.begin(write=True) as txn:
					for k, v in cachedfdict.items():
						txn.put(k, v)
						cachedfdict = {}
				with envpos.begin(write=True) as txn:
					for k, v in cachepos.items():
						txn.put(k, v)
						cachepos = {}
				with envori.begin(write=True) as txn:
					for k, v in cacheori.items():
						txn.put(k, v)
						cacheori = {}
	if cachedfdict:
		with envfinaldic.begin(write=True) as txn:
			for k, v in cachedfdict.items():
				txn.put(k, v)
	if cachepos:
		with envpos.begin(write=True) as txn:
			for k, v in cachepos.items():
				txn.put(k, v)
	if cacheori:
		with envori.begin(write=True) as txn:
			for k, v in cacheori.items():
				txn.put(k, v)
	pool.close()
	pool.join()
	return cacheseg



# with envpos.begin() as txn:
#     cursor = txn.cursor()  # 获取数据库的游标
#     all_data = {}  # 用来存储所有解包后的字典数据
#     for key, value in cursor:
#         unpacked_data = msgpack.unpackb(value)
#         all_data[key] = unpacked_data  # 将解包后的字典存储在 all_data 中，使用键作为字典的键


# xxxx=[]
# for val in all_data.values():
# 	xxxx.extend(val.keys())




parser = argparse.ArgumentParser(description='Pan-SD-idencordinate')
parser.add_argument('--W', help='graph path')
parser.add_argument('--S', help='corresponding sample name')
parser.add_argument('--p', help='work folder absolute path')
args = parser.parse_args()


sam=args.S
win=args.W
workdir=args.p

# python pansd.py --W ./uniqW/"$line".uniq.W --S "$line" --i ${inputallres} --o ./output/"$line"can_sd.corregion.bed --p /home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts
# sam="C001-CHA-E01-Mat"
# win="/home/jmhan/pangenome/WGSfea/allW/APG3allWnew/C001-CHA-E01-Mat.W"
# workdir="/home/jmhan/pan-SD/alldb/testsplitWdb"

envfinaldic = lmdb.open(workdir+'/db/'+sam+'/'+'finaldb', map_size=20 * 1024 * 1024 * 1024)  # 请根据需求修改路径和大小
envpos = lmdb.open(workdir+'/db/'+sam+'/'+'segpos', map_size=20 * 1024 * 1024 * 1024)  # 10 GB
envori = lmdb.open(workdir+'/db/'+sam+'/'+'segori', map_size=20 * 1024 * 1024 * 1024)  # 10 GB
with open(win, 'r') as f:
    lines = f.readlines()

env = lmdb.open(workdir+'/seg_db', map_size=20 * 1024 * 1024 * 1024,create=False,readonly=True, lock=False,max_readers=1000)  # 10 GB
nodedic=process_lines(lines, envfinaldic,envpos,envori)
print("SUCCESS!!!!!!!!")






