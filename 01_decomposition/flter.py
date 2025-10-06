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

with open('/share/home/zhanglab/user/maoyafei/project/panVNTR/data/CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
	fasta_dict = json.load(gz_file)


alpaircut=pd.read_csv("allpair.sort",sep="\t",header=None)
alpaircut['note']=0
alpaircut['node']=alpaircut[4].str.extract(r'(\d+)[><]')
alpaircut['group_id'] = (alpaircut['node'] != alpaircut['node'].shift()).cumsum() - 1
grouped = alpaircut.groupby('group_id')
alpaircut['len']=alpaircut[3]-alpaircut[2]







duplicate_pairs = []
dedup_dfs = []
for group_id, group_df in grouped:
	print(group_id)
	len_to_kept_row = {}
	for idx, row in group_df.iterrows():
		length = row['len']
		if length not in len_to_kept_row:
			len_to_kept_row[length] = row[1]
			duplicate_pairs.append((len_to_kept_row[length], row[1]))
		else:
			duplicate_pairs.append((len_to_kept_row[length], row[1]))
	dedup_df = group_df.drop_duplicates(subset=['len'], keep='first')
	dedup_dfs.append(dedup_df)  # 添加到列表

final_df = pd.concat(dedup_dfs, ignore_index=True)
df_pair = pd.DataFrame(duplicate_pairs) ##所有的pair都有了，之后要更新左边的区域，因为还有冗余


aldic={}
for index,row in final_df.iterrows():
	print(index)
	letters_to_find =[part for part in re.split(r'[<>]', row[4]) if part]
	aldic[index]=letters_to_find

grouped = final_df.groupby('group_id')



import multiprocessing
import networkx as nx
from itertools import combinations

def process_group(group):
	component=[]
	group_id, group_df = group
	print(group_id)
	new_dict = {key: aldic[key] for key in group_df.index}
	corr=[]
	for key1 in new_dict.keys():
		for key2 in new_dict.keys():
			if key2>key1:
				nowcho=set(aldic[key1])
				nextcho=set(aldic[key2])
				inter=nowcho & nextcho
				nowchofa=list(map(fasta_dict.get, nowcho))  #end.isnull().values.any()
				nextchofa=list(map(fasta_dict.get, nextcho))  #end.isnull().values.any()
				interfa=list(map(fasta_dict.get, inter))  #end.isnull().values.any()
				lena=sum([len(sublist) for sublist in nowchofa])
				lenb=sum([len(sublist) for sublist in nextchofa])
				leninter=sum([len(sublist) for sublist in interfa])
				corr.append([key1,key2,leninter/max(lena,lenb)])
	corre=[roww for roww in corr if roww[2] > 0.95]
	G = nx.Graph()
	for rowq in corre:
		G.add_edge(rowq[0], rowq[1])
	connected_components = list(nx.connected_components(G))
	component.extend(connected_components)
	return component

with multiprocessing.Pool(30) as pool:
    results = pool.map(process_group, grouped)
    

components = []
for group_components in results:
    components.extend(group_components)



dellist= set.union(*components)
allis = set(range(0, len(final_df)))
save1 = allis - dellist
save1df = final_df.loc[list(save1)] ###没有构成一个components的就不管了

first_elements = [next(iter(s)) for s in components if s]
save2df = final_df.loc[first_elements]

# 预先生成所有映射关系
mapping = {}
for com in components:
    print(com)
    nowindex = next(iter(com))
    chanin = final_df.loc[nowindex, 1]
    for idx in list(com):
        mapping[final_df.loc[idx, 1]] = chanin



# 批量更新df_pair
df_pair[0] = df_pair[1].map(mapping).fillna(df_pair[0]) 
df_pair[0] = df_pair[0].map(mapping).fillna(df_pair[0])


end=pd.concat([save1df, save2df])


##上面的第一轮为根据顶点进行过滤
##下面这轮继续去冗余

##按照样本分组
for sam in list(set(end[0])):
	print(sam)
	samdf=end[end[0]==sam]
	samdf['s']=samdf[1].str.split(':|-').str[-2]
	samdf['e']=samdf[1].str.split('-').str[-1]
	samdf['chr']=samdf[1].str.split(':').str[0]
	samdf.to_csv("/share/home/zhanglab/user/maoyafei/project/PANSDEND/END/02_expan_zd/"+sam+"/filt_split.path", sep="\t", index=False, header=False)


df_pair.to_csv("df.pair", sep="\t", index=False, header=False)

