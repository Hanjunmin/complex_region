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
sys.path.append('./project/ENDVNTR')
from vntr_fun import *

##./project/ENDVNTR
#####一次运行50个haplotype 200G 20 core

## less -S /share/home/project/zhanglab/APG/Graph/MC/CHM13-APGp1-HPRCp1-HGSVCp3/CHM13-APGp1-HPRCp1-HGSVCp3/CHM13-APGp1-HPRCp1-HGSVCp3_MC.vcf.gz  |cut -f 1-9 >1_9.all.vcf
# bedtools intersect -a <(cat /dssg/home/acct-clsmyf/clsmyf-user1/allW/vntrpri.bed |cut -f 1-3)  -b /dssg/home/acct-clsmyf/clsmyf-user1/allW/CHM13.chr1.Wallposend.tsv -wb|less -S |cut -f7>vntrinter.node
# for file in *.csv; do
#  	new_filename="${file%.csv}.txt"
#   	awk 'OFS="\t"{print $7, $4,$5,$1,$2}' "$file" |tail -n +2 >$new_filename
# done

parser = argparse.ArgumentParser(description='Pan-VNTR')
parser.add_argument('--CHR', help='chromsome')
parser.add_argument('--T2Tpos', help='segment of T2T')
parser.add_argument('--s', help='reruns')
parser.add_argument('--e', help='rerune')
args = parser.parse_args()


fafile="./LSGvar/genome/chm13v2.0.fa"
allWfold="./decompose/splinew/"
CHRn=args.CHR
T2Tchmpos=args.T2Tpos
samstart=int(args.s)
samend=int(args.e)
# trfpath="/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/Pan-VNTR/init/split/integrate_trf/out/chr2.vntr"
# CHRn="chr2"
# samstart=0
# T2Tchmpos="./project/panVNTR/data/chrY.posCHM.txt"
T2Tpos = pd.read_csv(T2Tchmpos,sep="\t",header=None)
T2Tcen = pd.read_csv("./project/panVNTR/data/hm_centroend.tsv",sep="\t")
cenlines=int(list(T2Tcen[T2Tcen['chr']==CHRn]['ref_start'])[0])
cenlinee=int(list(T2Tcen[T2Tcen['chr']==CHRn]['ref_end'])[0])
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



with open('./project/panVNTR/data/CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
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


# ##test一下
# def process_item(item):
# 	copy_number={}
# 	delnode=[]
# 	copy_numberlis={}
# 	alist=[]
# 	t2tmotiflist=[]
# 	i=0
# 	vntrid=item
# 	i=i+1
# 	storep=[]
# 	print(vntrid)
# 	vntridcopy_number={}
# 	result_dict = {}
# 	samall=canpos[vntrid]
# 	for single_dict in samall:
# 		result_dict.update(single_dict)
# 	for sam,wichnode in result_dict.items():
# 		pa=aldicmultipath[sam]
# 		pa=pa+"<"
# 		pos1, pos2 = find_specipos(wichnode[0], pa), find_specipos(wichnode[1], pa)
# 		corvntrpath=pa[min(min(pos1),min(pos2)):max(max(pos1),max(pos2))]
# 		symbols = re.findall(r'[<>]',corvntrpath)
# 		numbers = re.findall(r'\d+', corvntrpath)
# 		delnode.extend(numbers)
# 		end = pd.DataFrame({'fa': list(map(fasta_dict.get, numbers))}, index=numbers)
# 		end['ori']=symbols[0:len(end)]
# 		fa_series = end.loc[end['ori'] == '<', 'fa']
# 		fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
# 		end.loc[end['ori'] == '<', 'fa'] = fa_series
# 		newend=end
# 		newend['ori'] = newend['ori'].replace({'<': '-', '>': '+'})
# 		newend['node']=newend.index
# 		newend['merged'] = newend['node'] + newend['ori']
# 		allpa=','.join(newend['merged'])
# 		storep.append(['P',sam,allpa,"*"])
# 		samseq=''.join(end['fa'])
# 		ref=samseq
# 		if "CHM13v2" in sam:
# 			newend['T2Tpos'] = newend['node'].apply(lambda x: node_to_start_dict.get(x))
# 			newend['T2Tend'] = newend['T2Tpos'] + newend['fa'].str.len()
# 			t2ts=list(newend['T2Tpos'])[0]
# 			t2te=list(newend['T2Tend'])[-1]
# 			if len(ref)>100000:
# 				continue
# 			cnmax=0
# 			for motif in motifdic[vntrid]:
# 				resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
# 				expand_motif(resdict,motif, ref, interval,0.2,alllenlist)
# 				if alllenlist!=[]:
# 					motiflist=[]
# 					enall=cal_cn1(alllenlist,len(motif),0.3)
# 					for mo in enall[3]:
# 					    motiflist.append(ref[mo[0]:(mo[1]+1)])
# 					cn=enall[0]
# 					if cn>cnmax:
# 						cnmax=cn
# 						vntridcopy_number[sam]=cn
# 						enstatis=enall
# 						enmotiflist=motiflist
# 			if 'enstatis' in locals():
# 				t2tmotiflist.append([item,t2ts+enstatis[1]-1,t2ts+enstatis[2],cnmax,",".join(set(enmotiflist))])
# 				vntridcopy_number[sam]=cnmax
# 			if cnmax==0:
# 				vntridcopy_number[sam]="NONONO"
# 		else:
# 			if len(ref)>100000:
# 				continue
# 			cnmax=0
# 			cnmotif="a"
# 			for motif in motifdic[vntrid]:
# 				resdict,interval, motiflist, alllenlist1 = defaultdict(list),(0, len(ref)), [motif], []
# 				expand_motif(resdict,motif, ref, interval,0.2,alllenlist1)
# 				samseqrev=str(Seq(samseq).reverse_complement())
# 				ref=samseqrev
# 				resdict,interval, motiflist, alllenlist2 = defaultdict(list),(0, len(ref)), [motif], []
# 				expand_motif(resdict,motif, ref, interval,0.2,alllenlist2)
# 				if alllenlist1!=[]:
# 					enall=cal_cn1(alllenlist1,len(motif),0.3)
# 					print(enall[0])
# 					cn1=enall[0]
# 				if alllenlist2!=[]:
# 					enall=cal_cn1(alllenlist2,len(motif),0.3)
# 					print(enall[0])
# 					cn2=enall[0]
# 				if alllenlist1!=[] and alllenlist2!=[]:
# 					if max(cn1,cn2)>cnmax:
# 						cnmax=max(cn1,cn2)
# 						vntridcopy_number[sam]=max(cn1,cn2)
# 				if alllenlist1!=[] and alllenlist2==[]:
# 					if cn1>cnmax:
# 						cnmax=cn1
# 						cnmotif=motif
# 						vntridcopy_number[sam]=cn1
# 				if alllenlist1==[] and alllenlist2!=[]:
# 					if cn2>cnmax:
# 						cnmax=cn2
# 						cnmotif=motif
# 						vntridcopy_number[sam]=cn2
# 			if cnmax==0:
# 				vntridcopy_number[sam]="NONONO"
# 		copy_number[vntrid]=vntridcopy_number
# 	t2tmotiflist_unique = [list(item) for item in set(tuple(sublist) for sublist in t2tmotiflist)]
# 	delnode=set(delnode)
# 	return copy_number,storep,delnode,t2tmotiflist_unique

##test一下
def process_itemsub(item):
	delnode=[]
	i=0
	vntrid=item
	i=i+1
	print(vntrid)
	result_dict = {}
	samall=canpos[vntrid]
	for single_dict in samall:
		result_dict.update(single_dict)
	for sam,wichnode in result_dict.items():
		pa=aldicmultipath[sam]
		pa=pa+"<"
		pos1, pos2 = find_specipos(wichnode[0], pa), find_specipos(wichnode[1], pa)
		corvntrpath=pa[min(min(pos1),min(pos2)):max(max(pos1),max(pos2))]
		numbers = re.findall(r'\d+', corvntrpath)
		delnode.extend(numbers)
	delnode=set(delnode)
	return delnode



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

from collections import Counter

def process_other(target_value):
	out=[]
	# idd=idd+1
	# print(f"########process {idd}/{len(chrpathsam)}")
	print(target_value)
	df = pd.DataFrame(list(interdict[target_value].items()), columns=['Key', 'Value'])
	insert_list=aldicmulti[target_value]
	pa=aldicmultipath[target_value]
	symbols = re.findall(r'[<>]',pa)
	insert_df = pd.DataFrame({'Key': insert_list})
	df_combined = pd.merge(insert_df, df, on='Key', how='left')
	df_combined['fa'] = list(map(fasta_dict.get, df_combined['Key']))
	df_combined['newfa']=df_combined['fa']
	df_combined['ori'] = symbols
	df_combined['base_count'] = df_combined['fa'].apply(lambda x: len(x))
	df_combined['pos']=df_combined['base_count'].cumsum()
	df_combined['T2Tpos'] = df_combined['Key'].apply(lambda x: node_to_start_dict.get(x))
	df_combined['T2Tend'] = df_combined['T2Tpos']+df_combined['base_count']
	fa_series = df_combined.loc[df_combined['ori'] == '<', 'fa']
	fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
	df_combined.loc[df_combined['ori'] == '<', 'newfa'] = fa_series
	df_clean = df_combined.dropna()
	mark=1
	df_combined['mark']=0
	temp_list=[]
	for i in range(1,len(df_clean)):
		if df_clean.index[i] == df_clean.index[i-1] + 1:
			temp_list.extend([df_clean.index[i-1],df_clean.index[i]])
		else:
			temp_list.extend([df_clean.index[i-1],df_clean.index[i]])
			df_combined.loc[list(set(temp_list)),'mark']=mark
			temp_list=[df_clean.index[i]]
			mark=mark+1
	if temp_list:
		df_combined.loc[list(set(temp_list)), 'mark'] = mark
	now=0
	df_mark_nonzero=df_combined[df_combined['mark']!=0] ###1.相同T2Tseg下的deletion位置 大于20bp的搞出来
	import pyranges as pr
	x=df_mark_nonzero.iloc[:,[9,7,8]]
	x.columns=['Chromosome','Start','End']
	gr = pr.PyRanges(x)
	merged = gr.merge(slack=20) 
	merged_df = merged.df
	if len(merged_df)!=0:
		mark_counter = Counter(merged_df['Chromosome'])
		marks_gt_2 = [mark for mark, count in mark_counter.items() if count > 2]
		merged_df[merged_df['Chromosome'].isin(marks_gt_2)]
		merged_df["gap"] = merged_df.groupby("Chromosome")["Start"].shift(-1) - merged_df["End"]   ###有的染色体不一样！！
		merged_df=merged_df[merged_df['Chromosome'].isin(marks_gt_2)]
		gap_coords = (
			merged_df.assign(
				gap_chrom=merged_df["Chromosome"],
				gap_start=merged_df["End"],
				gap_end=merged_df.groupby("Chromosome")["Start"].shift(-1)
			)
			.dropna(subset=["gap_end"])  # 删除无 gap 的行
			.loc[:, ["gap_chrom", "gap_start", "gap_end"]]  # 只保留坐标列
		)
		for index,gap in gap_coords.iterrows():
			out.append([target_value,gap["gap_chrom"],gap["gap_start"],gap["gap_end"],0,"DEL"])
	df_mark_zero = df_combined[df_combined['mark'] == 0].copy()
	df_mark_zero['cluster'] = (df_mark_zero.index.to_series().diff() != 1).cumsum()
	df2hap=df2hapdf(df_combined,df_mark_zero)
	
	nan_rows = df2hap[df2hap['Value'].isna()].copy()
	nan_rows.loc[:, 'fa_length'] = nan_rows['newfa_concat'].apply(lambda x: len(x) if pd.notna(x) else 0)
	non_nan_rows = df2hap[~df2hap['Value'].isna()].copy()###*******
	# nanx=nan_rows[nan_rows['fa_length']>5]
	nanx=nan_rows
	fa = pysam.FastaFile(fafile)
	chrlength=fa.get_reference_length(CHRn)
	for nanclus in set(nanx['cluster']):
		labe11=False
		label2=False
		print(nanclus)
		cho=nanx[nanx['cluster']==nanclus]
		keyNod=list(cho['Key'])[0];keyNod1=list(cho['Key'])[-1]
		index1=cho.index[0]
		index2=cho.index[-1]
		now=now+1
		nowins=df2hap.loc[index1]
		insfa=nowins['newfa_concat']
		if index1!=0:
			lefmark=df2hap.loc[index1-1,"Key"]
			lefpos=df2hap.loc[index1-1,"T2Tpos"]
			lefend=df2hap.loc[index1-1,"T2Tend"]
		else:
			lefmark=-1
			lefpos=-1
			lefend=-1
		if index2!=len(df2hap)-1:
			rigmark=df2hap.loc[index2+1,"Key"]
			rigpos=df2hap.loc[index2+1,"T2Tpos"]
			rigend=df2hap.loc[index2+1,"T2Tend"]
		else:	
			rigmark=-2
			rigpos=-2
			rigend=-2
		if len(insfa)>20 or rigpos-lefend>20 or lefpos-rigend>20:
			out.append([target_value,lefmark,rigmark,lefpos,rigpos,insfa])
	return out



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
with open(CHRn+"_time", "a") as file:
	file.write(f"start ##### {datetime.now()} \n")

for sams in range(len(split_lists)): # for samlis in split_lists:
	if sams <samstart: ###samstart=1说明我要运行split_lists[1]
		continue
	if sams >samend: ###samstart=1说明我要运行split_lists[1]
		continue
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
	with mp.Pool(processes=30) as pool:
		results = pool.map(process_other, chrpathsam)
	for co in results:
		for res in co:
			graphvntr.append(res)
	with open(CHRn+"_time", "a") as file:
		g = psutil.Process(os.getpid())
		file.write(f"Job ID: {samlisid}/{len(split_lists)} :Memory usage(graph-p): {g.memory_info().rss / 1024 ** 3} GB##### {datetime.now()} \n")
	del aldicmulti
	del aldicmultipath
	del interdict
	graphvntrdf=(pd.DataFrame(graphvntr))
	allPAdf=(pd.DataFrame(allPA))
	graphvntrdf.to_csv(CHRn+"_"+str(samlisid)+'INSpos.csv', index=True,sep="\t")  # 不保存索引列
	del allPAdf
	del graphvntrdf
	del allPA
	del graphvntr
	with open(CHRn+"_time", "a") as file:
		g = psutil.Process(os.getpid())
		file.write(f"Job ID: {samlisid}/{len(split_lists)} :Memory usage(end): {g.memory_info().rss / 1024 ** 3} GB##### {datetime.now()} \n")
