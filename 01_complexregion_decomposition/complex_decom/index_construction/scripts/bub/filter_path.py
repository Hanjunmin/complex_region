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
parser = argparse.ArgumentParser(description='Filter the paths')
parser.add_argument('--i', help='path file')
parser.add_argument('--json', help='json file')
args = parser.parse_args()


with open(args.json, 'rt', encoding='utf-8') as gz_file:  ##/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.json
	fasta_dict = json.load(gz_file)

alpaircut=pd.read_csv(args.i,sep="\t",header=None)


# alpaircut=pd.read_csv("chr1:12333500-12757000.path.dou.bef",sep="\t",header=None)
alpaircut['note']=0
alpaircut['node']=alpaircut[2].str.extract(r'(\d+)[+-]')
alpaircut = alpaircut.sort_values('node')
alpaircut['group_id'] = (alpaircut['node'] != alpaircut['node'].shift()).cumsum() - 1
grouped = alpaircut.groupby('group_id')
alpaircut['len']=alpaircut[0]





aldic={}
for index,row in alpaircut.iterrows():
	print(index)
	letters_to_find =[part for part in re.split(",", row[2].replace("+","").replace("-","")) if part]
	aldic[index]=letters_to_find

grouped = alpaircut.groupby('group_id')

import multiprocessing
import networkx as nx
from itertools import combinations

import numpy as np
from scipy.sparse import csr_matrix
from itertools import combinations

# 构建序列到索引的映射
all_sequences = set()
for sequences in aldic.values():
    all_sequences.update(sequences)

sequence_to_idx = {seq: i for i, seq in enumerate(all_sequences)} ##赋予每个segment一个id




def process_group(group):
	component=[]
	group_id, group_df = group
	print(group_id)
	new_dict = {key: aldic[key] for key in group_df.index}
	keys = list(new_dict.keys())
	n_keys = len(keys)
	n_seqs = len(all_sequences)
	seq_lengths = np.array([len(fasta_dict[seq]) for seq in all_sequences]) ##每个segment赋予一个长度标识
	matrix = np.zeros((n_keys, n_seqs), dtype=bool) ##group对应的元素构建一个稀疏矩阵
	key_lengths = np.zeros(n_keys)
	for i, key in enumerate(keys):
		sequences = aldic[key]
		indices = [sequence_to_idx[seq] for seq in sequences]
		matrix[i, indices] = True
		key_lengths[i] = np.sum(seq_lengths[indices])
	corr = []
	for i, j in combinations(range(n_keys), 2):
		print(i,j)
		intersection = matrix[i] & matrix[j]
		if np.any(intersection):
			len_inter = np.sum(seq_lengths[intersection])
			max_len = max(key_lengths[i], key_lengths[j])
			corr.append([keys[i], keys[j], len_inter / max_len])
	corre=[roww for roww in corr if roww[2] > 0.97]
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
allis = set(alpaircut.index)
save1 = allis - dellist
save1df = alpaircut.loc[list(save1)] ###没有构成一个components的就不管了
first_elements = [next(iter(s)) for s in components if s]
save2df = alpaircut.loc[first_elements]

end=pd.concat([save1df, save2df])

###对应的pair的情况
pairs = []
for s in components:
    lst = list(s)
    if len(lst) > 1:
        a1 = lst[0]
        for a2 in lst:  
            pairs.append((a1, a2))

for ele in save1:
	pairs.append((ele, ele))


pairdf = pd.DataFrame(pairs, columns=['A1', 'A2'])
ID_dict = dict(zip(alpaircut.index, alpaircut[1]))
pairdf['A1_name'] = pairdf['A1'].map(ID_dict)
pairdf['A2_name'] = pairdf['A2'].map(ID_dict)
fin=end.iloc[:,:6]


fin.to_csv(args.i+".filt", sep="\t", index=False, header=False) ### >0.97的路径被过滤掉 之后改为一个输入的参数
pairdf.to_csv(args.i+".pair", sep="\t", index=False, header=False) ## 所有的路径被替换生成的pair
