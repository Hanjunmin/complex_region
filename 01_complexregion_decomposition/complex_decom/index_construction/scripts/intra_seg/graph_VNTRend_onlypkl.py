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
import sys
sys.path.append('/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/')
from vntr_fun import *



parser = argparse.ArgumentParser(description='Pan-VNTR')
# parser.add_argument('--CHR', help='chromsome')
parser.add_argument('--T2Tpos', help='segment of T2T')
parser.add_argument('--vntr', help='VNTR')
# parser.add_argument('--s', help='reruns')
# parser.add_argument('--e', help='rerune')
args = parser.parse_args()


fafile="/home/jmhan/breakpoints/chm13v2.0.fa"
allWfold="/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/Wpath/"
# trfpath=args.vntr
trfpath="can.VNTR.pos"
# CHRn=args.CHR
# T2Tchmpos=args.T2Tpos
T2Tchmpos="/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.posCHM.txt"
# samstart=int(args.s)
# samend=int(args.e)
T2Tpos = pd.read_csv(T2Tchmpos,sep="\t",header=None)
# T2Tcen = pd.read_csv("/share/home/zhanglab/user/maoyafei/data/chm13v2.0.refine.sat.bed",sep="\t",header=None)
# T2Tcen.columns=['chr','ref_start','ref_end']
# cenlines=int(list(T2Tcen[T2Tcen['chr']==CHRn]['ref_start'])[0])
# cenlinee=int(list(T2Tcen[T2Tcen['chr']==CHRn]['ref_end'])[0])
T2Tpos.columns=["CHR","start_position","end_position","node","orient"]
T2Tpos['node'] = T2Tpos['node'].astype(str)
node_to_start_dict = T2Tpos.set_index('node')['start_position'].to_dict() ##存放了T2T的node对应的位置信息
pos_to_start_dict = T2Tpos.set_index('start_position')['node'].to_dict() ##存放了T2T的node对应的位置信息
allnode=T2Tpos['node']
allnodeset=set(allnode)
######### allnodeset=set([str(node) for node in allnodeset])


############下面这个是HPRC的结果
### 下面的处理路径都在/home/jmhan/pangenome/WGSfea/SDsim/chr1/
### less -S chr1.vcf  |cut -f 1-4,8 |less -S |sed 's/AT=/\t/g'|sed 's/;NS/\t/g' |less -S |cut -f 1-4,6 |less -S |grep -v "^#" |less -S >chr1.simply
### less -S chr1.simply  |awk 'OFS="\t"{print $1,$2,$2+length($4)-1,$0}' |less -S >chr1.simply.bed
### bedtools intersect -a chr1.simply.bed -b /home/jmhan/breakpoints/trf/test/vntrpri.bed -wa|less -S >chr1.vntr.vcfcor
### less -S /home/jmhan/pangenome/WGSfea/allW/CHM13.Wallposend.tsv |awk '$1=="chr1"' >CHM13.chr1.Wallposend.tsv



with open('/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.json', 'rt', encoding='utf-8') as gz_file:
    fasta_dict = json.load(gz_file)


# vntrnode = pd.read_csv("vntrinter.node",sep="\t",header=None)
# vntrnodelist=list(vntrnode[0])




directory_path = allWfold
Wfile_list = []
for filename in os.listdir(directory_path):
    if filename.endswith('.W'):  # Ensure it's a file
        Wfile_list.append(directory_path+filename)

###上面这个有点大了，分成五个list运行吧
def split_list(lst, n):
    k, m = divmod(len(lst), n)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


##test一下



def df_fa(daf):
	fa_series = daf.loc[daf['ori'] == '<', 'fa']
	fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
	daf.loc[daf['ori'] == '<', 'fa'] = fa_series
	samseq=''.join(daf['fa'])
	return samseq










def find_pos(node, sampleW): ##这里的sampleW是一个样本的所有染色体组成的字典
    samwhich = []
    pathnow = []
    for path, samname in sampleW.items():
        nodelis=path.replace('<', ' ').replace('>', ' ').split()
        if node in nodelis:  # 修改为使用 `in` 检查
            pathnow.append(path)
            samwhich.append(samname)
            return pathnow, samwhich
    return []




def find_specipos(node,pa):
	plist=[">"+node+">",">"+node+"<","<"+node+"<","<"+node+">"]
	for nop in plist:
		whi=pa.find(nop)
		if whi!=-1:
			return [whi,whi+len(nop)]
	




def find_chr_pa(allnodeset, sampleW): ##这里的sampleW是一个样本的所有染色体组成的字典
    samchr=[]
    for path, samname in sampleW.items():
        nodelis=path.replace('<', ' ').replace('>', ' ').split()
        if len(allnodeset & set(nodelis))!=0:  # 修改为使用 `in` 检查
            samchr.append(samname)
    return samchr



# def test_intersect(lis):
# 	new_lis=list(range(lis[0][0],lis[0][1]+1))
# 	new_tup=[(lis[0])]
# 	for li in range(1,len(lis)):
# 		ne=list(range(lis[li][0],lis[li][1]+1))
# 		if len(set(ne)&set(new_lis))/len(set(ne))>=0.7:
# 			continue
# 		else:
# 			new_lis.extend(ne)
# 			new_tup.append((lis[li]))
# 	return new_tup
	




def process_file(file):
    print(file)
    with open(file, 'r') as f:
        lines = f.readlines()
    aldic = {}
    for line in lines:
        line_list = line.strip().split('\t')
        samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
        aldic[line_list[6]] = samname
    return aldic


def process_file2(file):
    print(file)
    with open(file, 'r') as f:
        lines = f.readlines()
    local_dict = defaultdict(list)  # 使用局部字典减少锁竞争
    for line in lines:
        line_list = line.strip().split('\t')
        samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
        if samname in chrpathsam:
            nodelist = (line_list[6]).replace('<', ' ').replace('>', ' ').split()
            local_dict[samname] = nodelist
    return local_dict

def process_file3(file):
    print(file)
    with open(file, 'r') as f:
        lines = f.readlines()
    local_dict = defaultdict(list)  # 使用局部字典减少锁竞争
    for line in lines:
        line_list = line.strip().split('\t')
        samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
        if samname in chrpathsam:
            local_dict[samname] = line_list[6]
    return local_dict



def process_file4(target_value,lines):
    print(target_value)
    values = np.fromiter(interdict[target_value].values(), dtype=int)
    sorted_values=np.sort(values)
    local_canpos = defaultdict(list)
    for line in lines:
        allduc={}
        split_line = line.strip().split("\t")
        CHR=split_line[0]
        if CHR!=CHRn:
            continue
        start=int(split_line[1])
        end=int(split_line[2])
        name=CHR+"_"+str(start)+"_"+str(end)
        startsam=np.searchsorted(sorted_values, start)
        if startsam==0 or startsam==len(sorted_values): ##所以相当于是忽略了边界因素？还没想清楚
            continue
        endsam=np.searchsorted(sorted_values, end)
        if endsam==0 or endsam==len(sorted_values):
            continue
        spos=sorted_values[startsam-1]
        epos=sorted_values[endsam]
        allduc[target_value]=[pos_to_start_dict[spos],pos_to_start_dict[epos]]
        local_canpos[name].append(allduc)
    return local_canpos

		

def process_sample(sampledic):
    result = find_chr_panew(allnodeset, sampledic)
    print(result)  # 确保缩进正确
    return result

def find_chr_panew(allnodeset, sampleW): ##这里的sampleW是一个样本的所有染色体组成的字典
    samchr=[]
    for path, samname in sampleW.items():
        nodelis=path.replace('<', ' ').replace('>', ' ').split()
        if not set(nodelis).isdisjoint(allnodeset):  # 修改为使用 `in` 检查
            samchr.append(samname)
    return samchr

def df2hapdf(df_combined,df_mark_zerodf):
	df_combined['intefa'] = None
	df_filtered = df_combined[df_combined['mark'] != 0]
	concatenated_values = (
		df_filtered.groupby('mark')['newfa']
		.apply(lambda x: ''.join(x.dropna()))  # Concatenate non-null 'newfa' values
	).reset_index()
	df_combined = df_combined.merge(concatenated_values, on='mark', how='left', suffixes=('', '_concat'))
	df_combined['intefa'] = df_combined['newfa_concat']
	df_combined.drop(columns=['newfa_concat'], inplace=True)
	df_combinedclus = df_combined.merge(df_mark_zerodf[['cluster']], left_index=True, right_index=True, how='left', suffixes=('', '_concat'))
	df_filtered = df_combinedclus[df_combinedclus['mark'] == 0]
	concatenated_values = (
		df_filtered.groupby('cluster')['newfa']
		.apply(lambda x: ''.join(x.dropna()))  # Concatenate non-null 'newfa' values
	).reset_index()
	df_combinedclus = df_combinedclus.merge(concatenated_values, on='cluster', how='left', suffixes=('', '_concat'))
	return df_combinedclus
	

def makefasta(input_file,allseq):
	with open(input_file, 'w') as f:
		f.write('>sequence\n')  # 描述行
		f.write(allseq+ '\n')
	trf_command = f'trf {input_file}  2 7 7 80 10 50 500 -h -l 20 -ngs'
	result = subprocess.run(['bash', '-c', trf_command], check=True, text=True, capture_output=True)
	result=result.stdout.split('\n')
	result=[item for item in result[1:] if item]
	return result

def findt2tloc(seq, T2T, T2Trev, t2ts):
    res1 = edlib.align(seq, T2T, mode='HW', task="path")
    l1="FALSE"
    for loc in res1['locations']:
        if (loc[1] - loc[0] + 1)/len(seq)>=0.95:
            l1="TRUE"
            t2tsn1 = t2ts + loc[0]
            t2ten1 = t2ts + loc[1]+1         
    res2 = edlib.align(seq, T2Trev, mode='HW', task="path")
    l2="FALSE"
    for loc in res2['locations']:
        if (loc[1] - loc[0] + 1)/len(seq)>=0.95:
            l2="TRUE"
            t2tsn2 = t2ts + (len(T2Trev)-loc[1])-1
            t2ten2 = t2ts + (len(T2Trev)-loc[0])
            break
    if l1=="TRUE" and l2=="TRUE":
        if res1['editDistance']>res2['editDistance']:
            return t2tsn2, t2ten2
        else:
            return t2tsn1, t2ten1
    if l1=="TRUE" and l2=="FALSE":
        return t2tsn1, t2ten1
    if l1=="FALSE" and l2=="TRUE":
        return t2tsn2, t2ten2
    if l1=="FALSE" and l2=="FALSE":
        return None
    
    







split_lists = split_list(Wfile_list, 100)
allCON = defaultdict(dict) 

for sams in range(len(split_lists)): # for samlis in split_lists:
	samlisid=sams
	allPA=[]
	graphvntr=[]
	samlis=split_lists[sams]
	Wfile_list=samlis
	alist=[] 
	with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
		futures = {executor.submit(process_file, file): file for file in Wfile_list}
		for future in concurrent.futures.as_completed(futures):
			result = future.result()
			alist.append(result)
	chrpathsam=[]
	with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
		futures = {executor.submit(process_sample, sampledic): sampledic for sampledic in alist}
		for future in concurrent.futures.as_completed(futures):
			res = future.result()
			chrpathsam.extend(res)  # 将结果扩展到chrpathsam
	del alist
	aldicmulti = defaultdict(list)
	with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
		futures = {executor.submit(process_file2, file): file for file in Wfile_list}
		for future in concurrent.futures.as_completed(futures):
			result = future.result()
			for key, value in result.items():
				aldicmulti[key] = value  # 合并结果到主字典
	aldicmultipath=defaultdict(list)
	with concurrent.futures.ProcessPoolExecutor(max_workers=20) as executor:
		futures = {executor.submit(process_file3, file): file for file in Wfile_list}
		for future in concurrent.futures.as_completed(futures):
			result = future.result()
			for key, value in result.items():
				aldicmultipath[key] = value  # 合并结果到主字典
	interdict={} ###存放每个样本与T2Tnode有交集的node以及对应的pos
	for target_value in set(chrpathsam):
		print(target_value)
		interset=allnodeset & set(aldicmulti[target_value]) ### target_value="HG01243_1_JAHEOY010000244.1_3031_869230"
		new_dict = {key: node_to_start_dict.get(key) for key in interset}
		interdict[target_value]=new_dict
	# trfpath="/share/home/zhanglab/user/maoyafei/project/ENDVNTR/STRVNTR130/chrY.STR_VNTR.out"
	with open(trfpath, 'r') as f:
		lines = f.readlines()
	motifdic={}
	for line in lines:
		split_line = line.strip().split("\t")
		CHR=split_line[0]
		start=int(split_line[1])
		end=int(split_line[2])
		name=CHR+"_"+str(start)+"_"+str(end)
		motifdic[name]=split_line[4].split(",")
	canpos = defaultdict(list)
	allcorresnode=[]
	for target_value in chrpathsam:
		print(target_value)
		values = np.fromiter(interdict[target_value].values(), dtype=int)
		sorted_values=np.sort(values)
		for line in lines:
			allduc={}
			split_line = line.strip().split("\t")
			CHR=split_line[0]
			start=int(split_line[1])
			end=int(split_line[2])
			name=CHR+"_"+str(start)+"_"+str(end)
			startsam=np.searchsorted(sorted_values, start)
			if startsam==0 or startsam==len(sorted_values): ##所以相当于是忽略了边界因素？还没想清楚
				continue
			endsam=np.searchsorted(sorted_values, end)
			if endsam==0 or endsam==len(sorted_values):
				continue
			spos=sorted_values[startsam-1]
			epos=sorted_values[endsam]
			if epos-spos>1000000:
				continue
			allduc[target_value]=[pos_to_start_dict[spos],pos_to_start_dict[epos]]
			canpos[name].append(allduc)
			allcorresnode.extend([pos_to_start_dict[spos],pos_to_start_dict[epos]])
	sampathpos={} ###生成所有样本在与VNTR相关区域的样本坐标
	for subsam in aldicmultipath.keys():
		symbols = re.findall(r'[<>]',aldicmultipath[subsam])
		numbers = re.findall(r'\d+', aldicmultipath[subsam])
		end = pd.DataFrame({'fa': list(map(fasta_dict.get, numbers))}, index=numbers)
		end['ori']=symbols
		end['base_count'] = end['fa'].apply(lambda x: len(x))
		end['pos']=end['base_count'].cumsum()
		end['pos_s']=end['pos']-end['base_count']
		filtend=end[end.index.isin(allcorresnode)]
		sampathpos[subsam]=filtend.iloc[:, [4, 3]].apply(tuple, axis=1).to_dict()
	alvntrsampos={} ###生成所有样本在与VNTR相关区域的真正坐标
	for vntrid in canpos.keys():
		alvntrsampossub={}
		result_dict = {}
		samall=canpos[vntrid]
		for single_dict in samall:
			result_dict.update(single_dict)
		for sam,wichnode in result_dict.items():
			sampos1=sampathpos[sam][wichnode[0]][0]
			sampos2=sampathpos[sam][wichnode[1]][1]
			min_val = min(sampathpos[sam][wichnode[0]][0], sampathpos[sam][wichnode[0]][1],
              sampathpos[sam][wichnode[1]][0], sampathpos[sam][wichnode[1]][1])
			max_val = max(sampathpos[sam][wichnode[0]][0], sampathpos[sam][wichnode[0]][1],
			              sampathpos[sam][wichnode[1]][0], sampathpos[sam][wichnode[1]][1])
			alvntrsampossub[sam]=[min_val,max_val]
		alvntrsampos[vntrid]=alvntrsampossub
	with open(str(samlisid)+"pos.pkl", "wb") as f:
		pickle.dump(alvntrsampos, f)
	