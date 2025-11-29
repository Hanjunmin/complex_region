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
# sys.path.append('/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/')
# from vntr_fun import *



parser = argparse.ArgumentParser(description='Pan-VNTR')
# parser.add_argument('--CHR', help='chromsome')
parser.add_argument('--script', help='script of VNTR calculation')
parser.add_argument('--T2Tpos', help='segment of T2T')
parser.add_argument('--vntr', help='VNTR')
parser.add_argument('--r', help='reference')
parser.add_argument('--w', help='the fold of W path')
parser.add_argument('--json', help='the file of json')
# parser.add_argument('--s', help='reruns')
# parser.add_argument('--e', help='rerune')
args = parser.parse_args()

sys.path.append(args.script)  ##'/home/jmhan/project/APG/github_final/01_1p36/index_construction/script/intra_seg/'
from vntr_fun import *

fafile=args.r #"/home/jmhan/breakpoints/chm13v2.0.fa"
allWfold=args.w #"/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/Wpath/"
CHRn="chr1"   ##这个需要改！！！
# trfpath=args.vntr
trfpath="can.VNTR.pos"
# CHRn=args.CHR
# T2Tchmpos=args.T2Tpos
T2Tchmpos=args.T2Tpos #"/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.posCHM.txt"
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



with open(args.json, 'rt', encoding='utf-8') as gz_file: #'/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.json'
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
def process_item(item):
	copy_number={}
	delnode=[]
	copy_numberlis={}
	alist=[]
	t2tmotiflist=[]
	i=0
	vntrid=item
	i=i+1
	storep=[]
	print(vntrid)
	vntridcopy_number={}
	result_dict = {}
	samall=canpos[vntrid]
	for single_dict in samall:
		result_dict.update(single_dict)
	for sam,wichnode in result_dict.items():
		pa=aldicmultipath[sam]
		pa=pa+"<"
		pos1, pos2 = find_specipos(wichnode[0], pa), find_specipos(wichnode[1], pa)
		corvntrpath=pa[min(min(pos1),min(pos2)):max(max(pos1),max(pos2))]
		symbols = re.findall(r'[<>]',corvntrpath)
		numbers = re.findall(r'\d+', corvntrpath)
		delnode.extend(numbers)
		end = pd.DataFrame({'fa': list(map(fasta_dict.get, numbers))}, index=numbers)
		end['ori']=symbols[0:len(end)]
		fa_series = end.loc[end['ori'] == '<', 'fa']
		fa_series = fa_series.apply(lambda x: str(Seq(x).reverse_complement()))
		end.loc[end['ori'] == '<', 'fa'] = fa_series
		newend=end
		newend['ori'] = newend['ori'].replace({'<': '-', '>': '+'})
		newend['node']=newend.index
		newend['merged'] = newend['node'] + newend['ori']
		allpa=','.join(newend['merged'])
		storep.append(['P',sam,allpa,"*"])
		samseq=''.join(end['fa'])
		ref=samseq
		if "CHM13v2" in sam:
			newend['T2Tpos'] = newend['node'].apply(lambda x: node_to_start_dict.get(x))
			newend['T2Tend'] = newend['T2Tpos'] + newend['fa'].str.len()
			t2ts=list(newend['T2Tpos'])[0]
			t2te=list(newend['T2Tend'])[-1]
			if len(ref)>100000:
				continue
			cnmax=0
			for motif in motifdic[vntrid]:
				resdict,interval, motiflist, alllenlist = defaultdict(list),(0, len(ref)), [motif], []
				expand_motif(resdict,motif, ref, interval,0.2,alllenlist)
				if alllenlist!=[]:
					motiflist=[]
					enall=cal_cn1(alllenlist,len(motif),0.3)
					for mo in enall[3]:
					    motiflist.append(ref[mo[0]:(mo[1]+1)])
					cn=enall[0]
					if cn>cnmax:
						cnmax=cn
						vntridcopy_number[sam]=cn
						enstatis=enall
						enmotiflist=motiflist
			if 'enstatis' in locals():
				t2tmotiflist.append([item,t2ts+enstatis[1]-1,t2ts+enstatis[2],cnmax,",".join(set(enmotiflist))])
				vntridcopy_number[sam]=cnmax
			if cnmax==0:
				vntridcopy_number[sam]="NONONO"
		else:
			if len(ref)>100000:
				continue
			cnmax=0
			cnmotif="a"
			for motif in motifdic[vntrid]:
				resdict,interval, motiflist, alllenlist1 = defaultdict(list),(0, len(ref)), [motif], []
				expand_motif(resdict,motif, ref, interval,0.2,alllenlist1)
				samseqrev=str(Seq(samseq).reverse_complement())
				ref=samseqrev
				resdict,interval, motiflist, alllenlist2 = defaultdict(list),(0, len(ref)), [motif], []
				expand_motif(resdict,motif, ref, interval,0.2,alllenlist2)
				if alllenlist1!=[]:
					enall=cal_cn1(alllenlist1,len(motif),0.3)
					print(enall[0])
					cn1=enall[0]
				if alllenlist2!=[]:
					enall=cal_cn1(alllenlist2,len(motif),0.3)
					print(enall[0])
					cn2=enall[0]
				if alllenlist1!=[] and alllenlist2!=[]:
					if max(cn1,cn2)>cnmax:
						cnmax=max(cn1,cn2)
						vntridcopy_number[sam]=max(cn1,cn2)
				if alllenlist1!=[] and alllenlist2==[]:
					if cn1>cnmax:
						cnmax=cn1
						cnmotif=motif
						vntridcopy_number[sam]=cn1
				if alllenlist1==[] and alllenlist2!=[]:
					if cn2>cnmax:
						cnmax=cn2
						cnmotif=motif
						vntridcopy_number[sam]=cn2
			if cnmax==0:
				vntridcopy_number[sam]="NONONO"
		copy_number[vntrid]=vntridcopy_number
	t2tmotiflist_unique = [list(item) for item in set(tuple(sublist) for sublist in t2tmotiflist)]
	delnode=set(delnode)
	return copy_number,storep,delnode,t2tmotiflist_unique

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
        # if CHR!=CHRn:
        #     continue
        CHRn=CHR
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

def process_other(target_value):
	# idd=idd+1
	# print(f"########process {idd}/{len(chrpathsam)}")
	othervntr=[]
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
	df_mark_zero = df_combined[df_combined['mark'] == 0].copy()
	df_mark_zero['cluster'] = (df_mark_zero.index.to_series().diff() != 1).cumsum()
	df2hap=df2hapdf(df_combined,df_mark_zero)
	nan_rows = df2hap[df2hap['Value'].isna()].copy()
	nan_rows.loc[:, 'fa_length'] = nan_rows['newfa_concat'].apply(lambda x: len(x) if pd.notna(x) else 0)
	non_nan_rows = df2hap[~df2hap['Value'].isna()].copy()###*******
	nanx=nan_rows[nan_rows['fa_length']>5]
	# nanx=nanx[~nanx['Key'].isin(del_region)]
	fa = pysam.FastaFile(fafile)
	for nanclus in set(nanx['cluster']):
		labe11=False
		label2=False
		print(nanclus)
		cho=nanx[nanx['cluster']==nanclus]
		keyNod=list(cho['Key'])[0];keyNod1=list(cho['Key'])[-1]
		keylis=",".join(cho['Key'])
		index1=cho.index[0]
		index2=cho.index[-1]
		now=now+1
		nowins=df2hap.loc[index1]
		if index1!=0:
			if index1-1!=0:
				reflef_mark=df2hap.loc[index1-1,'mark']
				reflef_seg=df2hap[df2hap['mark']==reflef_mark].copy()
				reflef_seg['T2Tpos'] = reflef_seg['Key'].apply(lambda x: node_to_start_dict.get(x))
				lefseq=list(reflef_seg['intefa'])[0]
				reflef_seg['end'] = reflef_seg['T2Tpos'] + reflef_seg['fa'].str.len()
				if reflef_mark==0:
					lefseq="N"
				# else:
					# if any(cenlines <= int(point) <= cenlinee for point in reflef_seg['T2Tpos']): ##之后用户输入这个可以加进去 
					# 	lefseq="N"
			else:
				lefseq="N"
		else:
			lefseq="N"
		if index2!=(len(df2hap)-1):
			refrig_mark=df2hap.loc[index2+1,'mark']
			refrig_seg=df2hap[df2hap['mark']==refrig_mark].copy()
			refrig_seg['T2Tpos']=refrig_seg['Key'].apply(lambda x: node_to_start_dict.get(x))
			rightseq=list(refrig_seg['intefa'])[0]
			refrig_seg['end'] = refrig_seg['T2Tpos'] + refrig_seg['fa'].str.len()
			if refrig_mark==0:
				rightseq="N"
			# else:
				# if any(cenlines <= int(point) <= cenlinee for point in refrig_seg['T2Tpos']): ##之后用户输入这个可以加进去 
				# 	rightseq="N"
		else:
			rightseq="N"
		insfa=nowins['newfa_concat']
		remask="NO"
		print("lef")
		if len(insfa)<5000:
			print("OK")
			if len(lefseq)>=len(insfa):
				lefseqnew=lefseq[-(len(insfa)+20):] ##取最后一段的序列
			if len(rightseq)>=len(insfa) : ## 这里我希望写一种情况，就是当插入序列很长的时候，我希望真正可以找到对应的区域
				rightseqnew=rightseq[:(len(insfa)+20)]  ##取最后一段的序列
			if lefseq=="N":
				lefseqnew=lefseq
			if rightseq=="N":
				rightseqnew=rightseq
			if len(lefseq)<len(insfa) and lefseq!="N":
				cumulative_sum=0 ###********************************************************************
				filtered_reflef= df2hap[df2hap.index < reflef_seg.index[0]]
				if len(filtered_reflef)!=0:
					stillnd=len(insfa)-len(lefseq)
					result_rows=[]
					for idx, row in filtered_reflef[::-1].iterrows():  # 倒序遍历
						cumulative_sum += row['base_count']
						result_rows.append(row)
						if cumulative_sum >= stillnd:
							break
					result_df = pd.DataFrame(result_rows)
					result_df = result_df.sort_index()
					result_df['T2Tpos'] = result_df['Key'].apply(lambda x: node_to_start_dict.get(x))
					result_df['end'] = result_df['T2Tpos'] + result_df['fa'].str.len()
					reflef_segnew = pd.concat([result_df, reflef_seg], axis=0)
				else:
					reflef_segnew=reflef_seg
				lefseqnew=''.join(reflef_segnew['newfa']) ###********************************************************************
				reflef_seg=reflef_segnew
			if len(rightseq)<len(insfa) and rightseq!="N":
				cumulative_sum=0 ###********************************************************************
				filtered_refrig= df2hap[df2hap.index >refrig_seg.index[-1]]
				if len(filtered_refrig)!=0:
					stillnd=len(insfa)-len(rightseq)
					result_rows=[]
					for idx, row in filtered_refrig[::1].iterrows():  # 倒序遍历
						cumulative_sum += row['base_count']
						result_rows.append(row)
						if cumulative_sum >= stillnd:
							break
					result_df = pd.DataFrame(result_rows)
					result_df = result_df.sort_index()
					result_df['T2Tpos'] = result_df['Key'].apply(lambda x: node_to_start_dict.get(x))
					result_df['end'] = result_df['T2Tpos'] + result_df['fa'].str.len()
					refrig_segnew = pd.concat([refrig_seg,result_df], axis=0)
				else:
					refrig_segnew=refrig_seg
				rightseqnew=''.join(refrig_segnew['newfa']) ###********************************************************************
				refrig_seg=refrig_segnew
			input_file=CHRn+target_value+".fa"
			allseq=lefseqnew+insfa+rightseqnew
			splice=[len(lefseqnew),len(lefseqnew+insfa)]
			refseq=lefseqnew+rightseqnew
			if lefseq!="N":
				result_rows=[]
				cumulative_sum=0 
				for idx, row in reflef_seg[::-1].iterrows():  # 倒序遍历
					cumulative_sum += row['base_count']
					result_rows.append(row)
					if cumulative_sum >= len(lefseqnew):
						break
				result_dfnew1 = pd.DataFrame(result_rows)
				reflef_segn=result_dfnew1
			if rightseq!="N":
				result_rows=[]
				cumulative_sum=0 
				for idx, row in refrig_seg[::1].iterrows():  # 倒序遍历
					cumulative_sum += row['base_count']
					result_rows.append(row)
					if cumulative_sum >= len(rightseqnew):
						break
				result_dfnew2 = pd.DataFrame(result_rows)
				refrig_segn=result_dfnew2
			result=makefasta(input_file,allseq)
			os.remove(input_file)  # 删除文件
			print("OK1")
			for ran in range(0,len(result)):
				trflis=result[ran].split(' ')
				if not (int(trflis[1]) <= splice[0] or int(trflis[0]) >= splice[1] or (int(trflis[0]) >= splice[0] and int(trflis[1]) <= splice[1])):
					motif=trflis[13]
					cn=trflis[3]	
					if lefseq!="N" and rightseq!="N":
						a=(refrig_segn['T2Tpos'].dropna().tolist() + reflef_segn['T2Tpos'].dropna().tolist())
						print(a)
						b=(refrig_segn['end'].dropna().tolist() + reflef_segn['end'].dropna().tolist())
						print(b)
						t2tpos=[min(a),max(b)]
					if lefseq=="N" and rightseq!="N":
						t2tpos=[min(refrig_segn['T2Tpos'].dropna().tolist()),max(refrig_segn['end'].dropna().tolist())]
					if lefseq!="N" and rightseq=="N":
						t2tpos=[min(reflef_seg['T2Tpos'].dropna().tolist()),max(reflef_seg['end'].dropna().tolist())]
					if lefseq=="N" and rightseq=="N":
						break
					#t2tpos=[min(min(refrig_segn['T2Tpos']),min(reflef_segn['T2Tpos'])),max(max(refrig_segn['end']),max(reflef_segn['end']))]
					name =CHRn + ":" + str(int(t2tpos[0])) + "-" + str(int(t2tpos[1]))
					Seq1lef = (fa.fetch(region=name)).upper()
					Seq1lefrev=str(Seq(Seq1lef).reverse_complement())
					res1 = edlib.align(motif, refseq, mode='HW', task="path")
					if res1['editDistance'] <= 0.3*len(motif) and (sorted(res1['locations'])[-1][1]-sorted(res1['locations'])[0][0])/len(motif)>=0.8 and (sorted(res1['locations'])[-1][1]-sorted(res1['locations'])[0][0])/len(motif)<2:
						findseq=refseq[sorted(res1['locations'])[0][0]:sorted(res1['locations'])[-1][1]+1]
						xxx=findt2tloc(findseq,Seq1lef,Seq1lefrev,t2tpos[0])
						if xxx==None:
							continue
						else:
							addlist=[CHRn,xxx[0],xxx[1],trflis[2],motif,cn,target_value,lefseqnew+"*******"+insfa+"*******"+rightseqnew,insfa,nanclus,splice,keyNod,keylis]
							othervntr.append(addlist)
				if (int(trflis[0]) >= splice[0] and int(trflis[1]) <= splice[1]):
					motif=trflis[13]
					cn=trflis[3]
					if lefseq!="N":
						addlist=[CHRn,list(reflef_seg['end'])[-1],"INS",trflis[2],motif,cn,target_value,lefseqnew+"*******"+insfa+"*******"+rightseqnew,insfa,nanclus,splice,keyNod,keylis]
						othervntr.append(addlist)
					if rightseq!="N" and lefseq=="N":
						addlist=[CHRn,list(refrig_seg['T2Tpos'])[0],"INS",trflis[2],motif,cn,target_value,lefseqnew+"*******"+insfa+"*******"+rightseqnew,insfa,nanclus,splice,keyNod,keylis]
						othervntr.append(addlist)
		else:
			input_file=CHRn+target_value+".fa"
			insfanew1=insfa[:5000]
			lefseqnew=lefseq[-5000:]
			insfanew2=insfa[-5000:]
			rightseqnew=rightseq[:5000]
			allseq1=lefseqnew+insfanew1
			splice1=[len(lefseqnew)]
			allseq2=insfanew2+rightseqnew
			# refseq=lefseqnew+rightseqnew
			splice2=[len(insfanew2)]
			if lefseq!="N":
				result_rows=[]
				cumulative_sum=0 
				for idx, row in reflef_seg[::-1].iterrows():  # 倒序遍历
					cumulative_sum += row['base_count']
					result_rows.append(row)
					if cumulative_sum >= len(lefseqnew):
						break
				result_dfnew1 = pd.DataFrame(result_rows)
				reflef_segn=result_dfnew1
			if rightseq!="N":
				result_rows=[]
				cumulative_sum=0 
				for idx, row in refrig_seg[::1].iterrows():  # 倒序遍历
					cumulative_sum += row['base_count']
					result_rows.append(row)
					if cumulative_sum >= len(rightseqnew):
						break
				result_dfnew2 = pd.DataFrame(result_rows)
				refrig_segn=result_dfnew2
			result1=makefasta(input_file,allseq1)
			os.remove(input_file)  # 删除文件
			result2=makefasta(input_file,allseq2)
			os.remove(input_file)  # 删除文件
			print("OK1")
			for ran in range(0,len(result1)):
				trflis=result1[ran].split(' ')
				if ((int(trflis[0]) <= splice1[0] and int(trflis[1]) >= splice1[0])):
					motif=trflis[13]
					cn=trflis[3]
					if lefseq!="N":
						a=reflef_segn['T2Tpos'].dropna().tolist()
						b=reflef_segn['end'].dropna().tolist()
						t2tpos=[min(a),max(b)]
					if lefseq=="N":
						break
					fa = pysam.FastaFile(fafile)
					name =CHRn + ":" + str(int(t2tpos[0])) + "-" + str(int(t2tpos[1]))
					addlist=[CHRn,t2tpos[0],t2tpos[1],trflis[2],motif,cn,target_value,lefseqnew+"*******"+insfanew1+"*******",insfanew1,nanclus,splice1,keyNod,keylis]
					othervntr.append(addlist)
				if ((int(trflis[0]) >= splice1[0])):
					motif=trflis[13]
					cn=trflis[3]
					if lefseq!="N":
						a=reflef_segn['T2Tpos'].dropna().tolist()
						b=reflef_segn['end'].dropna().tolist()
						t2tpos=[min(a),max(b)]
					if lefseq=="N":
						break
					fa = pysam.FastaFile(fafile)
					name =CHRn + ":" + str(int(t2tpos[0])) + "-" + str(int(t2tpos[1]))
					addlist=[CHRn,t2tpos[0],"INS",trflis[2],motif,cn,target_value,lefseqnew+"*******"+insfanew1+"*******",insfanew1,nanclus,splice1,keyNod,keylis]
					othervntr.append(addlist)
			for ran in range(0,len(result2)):
				trflis=result2[ran].split(' ')
				if ((int(trflis[0]) <= splice2[0] and int(trflis[1]) >= splice2[0])):
					motif=trflis[13]
					cn=trflis[3]
					if rightseq!="N":
						a=(refrig_segn['T2Tpos'].dropna().tolist())
						b=(refrig_segn['end'].dropna().tolist())
						t2tpos=[min(a),max(b)]
					if rightseq=="N":
						break
					fa = pysam.FastaFile(fafile)
					name =CHRn + ":" + str(t2tpos[0]) + "-" + str(t2tpos[1])
					addlist=[CHRn,t2tpos[0],t2tpos[1],trflis[2],motif,cn,target_value,"*******"+insfanew2+"*******"+rightseqnew,insfanew2,nanclus,splice2,keyNod1,keylis]
					othervntr.append(addlist)
				if ((int(trflis[1]) <= splice2[0])):
					motif=trflis[13]
					cn=trflis[3]
					if rightseq!="N":
						a=(refrig_segn['T2Tpos'].dropna().tolist())
						b=(refrig_segn['end'].dropna().tolist())
						t2tpos=[min(a),max(b)]
					if rightseq=="N":
						break
					fa = pysam.FastaFile(fafile)
					name =CHRn + ":" + str(t2tpos[0]) + "-" + str(t2tpos[1])
					addlist=[CHRn,t2tpos[0],"INS",trflis[2],motif,cn,target_value,"*******"+insfanew2+"*******"+rightseqnew,insfanew2,nanclus,splice2,keyNod1,keylis]
					othervntr.append(addlist)
	return othervntr

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
	with mp.Pool(processes=30) as pool:
		results = pool.map(process_other, chrpathsam)
	for co in results:
		for res in co:
			graphvntr.append(res)
	# del aldicmulti
	# del aldicmultipath
	# del interdict
	graphvntrdf=(pd.DataFrame(graphvntr))
	allPAdf=(pd.DataFrame(allPA))
	graphvntrdf.to_csv("reg_"+str(samlisid)+'other.copynumber.csv', index=True,sep="\t")  # 不保存索引列
	# del allPAdf
	# del graphvntrdf
	# del allPA
	# del graphvntr
	
