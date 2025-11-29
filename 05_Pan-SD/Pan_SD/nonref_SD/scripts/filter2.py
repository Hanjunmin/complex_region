##从这里开始
import pandas as pd
import networkx as nx
import re
import argparse
import csv
import re
import numpy as np
import pandas as pd
import json
import os
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser(description='filter the path')
parser.add_argument('--s', help='samlis path')
parser.add_argument('--b', help='base path')
args = parser.parse_args()
base_dir=args.b
samlis=args.s


import pandas as pd
samdf = pd.read_csv(samlis,sep="\t",header=None)
samname=list(samdf[samdf[1]!="reference"][0])

# with open('/share/home/zhanglab/user/maoyafei/project/panVNTR/data/CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
# 	fasta_dict = json.load(gz_file)

with open(base_dir+"/00_preproc/reg.json", 'rt', encoding='utf-8') as gz_file:
	fasta_dict = json.load(gz_file)



df_pair=pd.read_csv("df.pair", sep="\t",header=None)
	

###去掉自身冗余的的pair
filtal=[]
mapping = {}
for sam in samname:
	pa=base_dir+"/03_nonref/02_ITER/2_SD_GENE/"+sam+"/processend/y.temp"
	samdfpa=base_dir+"/03_nonref/02_ITER/2_SD_GENE/"+sam+"/processend/filt_split.path"
	if os.path.exists(samdfpa) and os.path.getsize(samdfpa) > 0:
		samdf=pd.read_csv(base_dir+"/03_nonref/02_ITER/2_SD_GENE/"+sam+"/processend/filt_split.path", sep="\t",header=None)
	if os.path.exists(pa) and os.path.getsize(pa) > 0:
		ytemp=pd.read_csv(pa, sep="\t",header=None)
		G = nx.Graph()
		for index,ytempr in ytemp.iterrows():
			G.add_edge(ytempr[0]+":"+str(ytempr[1])+"-"+str(ytempr[2]), ytempr[3]+":"+str(ytempr[4])+"-"+str(ytempr[5]))
		connected_components = list(nx.connected_components(G))
		for com in connected_components:
			print(com)
			nownam = next(iter(com))
			for diname in list(com):
				mapping[diname] = nownam
		first_elements = [next(iter(s)) for s in connected_components if s]
		dellist= set.union(*connected_components)
		deplpath=dellist-set(first_elements)
		filtered_df = samdf[~samdf[1].isin(deplpath)]
		filtal.append(filtered_df)
	else:
		filtal.append(samdf)

df_pair[0] = df_pair[1].map(mapping).fillna(df_pair[0]) 
df_pair[0] = df_pair[0].map(mapping).fillna(df_pair[0])

result = pd.concat(filtal, ignore_index=True)
###去掉之间冗余的pair
import random
aldic={}
alldict10={}
for sam in list(set(result[0])):
	print(sam)
	samdic={}
	samdic10={}
	subsa=result[result[0]==sam]
	for index,row in subsa.iterrows():
		letters_to_find =[part for part in re.split(r'[<>]', row[4]) if part]
		samdic[index]=letters_to_find
		if len(letters_to_find) >= 30:
			samdic10[index]=set(random.sample(letters_to_find, 30))
		else:
			samdic10[index]=set(letters_to_find)
	aldic[sam]=samdic
	alldict10[sam]=samdic10


from itertools import combinations
alcompos=list(combinations(list(set(result[0])), 2)) ###所有的样本的组合
i=1
for comp in alcompos:
	i=i+1
	print(i)
	sam1=comp[0]
	sam2=comp[1]
	A=alldict10[sam1]
	B=alldict10[sam2]
	keys_to_delete=[]
	for key, value in A.items():
		intersam=[sample for sample, path in B.items() if len(set(value) & set(path))]
		# print(intersam)
		if intersam:
			for internx in intersam:
				nowcho=set(aldic[sam1][key])
				nextcho=set(aldic[sam2][internx])
				inter=nowcho & nextcho
				nowchofa=list(map(fasta_dict.get, nowcho))  #end.isnull().values.any()
				nextchofa=list(map(fasta_dict.get, nextcho))  #end.isnull().values.any()
				interfa=list(map(fasta_dict.get, inter))  #end.isnull().values.any()
				lena=sum([len(sublist) for sublist in nowchofa])
				lenb=sum([len(sublist) for sublist in nextchofa])
				leninter=sum([len(sublist) for sublist in interfa])
				ninterratio=leninter/max(lena,lenb)
				if ninterratio>=0.95:
					print("yes")
					print(inter)
					print(key)
					deln=result.loc[key,1]
					chan=result.loc[internx,1]
					df_pair.loc[df_pair[0]==deln,0]=chan
					keys_to_delete.append(key)
	for key in set(keys_to_delete):
		del aldic[sam1][key]
		del alldict10[sam1][key]



saveal=[]
for sam in list(set(result[0])):
	saveal.extend(aldic[sam])

resultal=result.loc[saveal]
#resultal=resultal.iloc[:,1:5]
resultal.to_csv("allpair.sort.END", sep="\t", index=False, header=False)
df_pair.to_csv("df.pair.END", sep="\t", index=False, header=False)



