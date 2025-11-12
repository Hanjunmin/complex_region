import numpy as np
import re
import math
import Levenshtein
import pandas as pd
import pysam
import multiprocessing as mp
import edlib
from collections import defaultdict
import argparse
import os
import subprocess

def expand_motif(result_dict,motif, ref, interval,infector,alllenlist):
	start,end = interval
	seq = (ref[start:(end + 1)]).upper()
	res1 = edlib.align(motif, seq, mode='HW', task="path")  # align with motif 
	if res1['editDistance'] <= infector*len(motif):	
		res1['locations']=test_intersect(res1['locations'])
		ninterval=[(start + x[0], start + x[1]) for x in res1['locations']]
		alllenlist.append(ninterval)
		_ = [result_dict[ref[x[0]:x[1]+1]].append(x) for x in ninterval]
		intervallist=find_non_continuous_intervals(start,end,ninterval,len(motif))
		for interval in intervallist:
			expand_motif(result_dict,motif, ref, interval,infector,alllenlist)

def test_intersect(lis):
	new_lis=list(range(lis[0][0],lis[0][1]+1))
	new_tup=[(lis[0])]
	for li in range(1,len(lis)):
		ne=list(range(lis[li][0],lis[li][1]+1))
		if len(set(ne)&set(new_lis))/len(set(ne))>=0.7:
			continue
		else:
			new_lis.extend(ne)
			new_tup.append((lis[li]))
	return new_tup


def tuple_process(sorted_tuples,motif_len):
	result_dict = defaultdict(list)
	id=0
	result_dict[id].append(sorted_tuples[0])
	for tupid in range(1,len(sorted_tuples)):
		if abs(sorted_tuples[tupid][0]-sorted_tuples[tupid-1][1])<=math.ceil(motif_len * 0.20):
			result_dict[id].append(sorted_tuples[tupid])
		else:
			id=id+1
			result_dict[id].append(sorted_tuples[tupid])
	return result_dict

def find_non_continuous_intervals(refstart,refend,points,len_motif):
    intervals = []
    if points[0][0]!=refstart:
    	intervals.append((refstart,points[0][0]-1))
    for i in range(1, len(points)):
        if abs(points[i][0] - points[i - 1][1])>len_motif*0.5: ### 这里的10后续需要根据VNTR（STR）的长度调整
            intervals.append((points[i - 1][1]+1,points[i][0]-1)) ### 这里的结果都是左闭右闭的，在构建序列的时候需要右边加1 
    if points[-1][1]!=refend:
    	intervals.append((points[-1][1]+1,refend))
    return intervals



def cal_cn1(alllenlist,motif_len,tolerance):
	all_tuples = [tup for sublist in alllenlist for tup in sublist]
	sorted_tuples = sorted(all_tuples)
	endtup=tuple_process(sorted_tuples,motif_len)
	endtupall=sorted(integrate_fragmentedmotif(endtup,motif_len,tolerance))
	cn=len(endtupall)
	news=endtupall[0][0]
	newe=endtupall[-1][1]
	return cn,news,newe,endtupall


def process_vamos(store):
    storelistsort=sorted(store,reverse=True)
    m=storelistsort[0]
    havestore=list(range(m[1],m[2]))
    save=[]
    save.append(m)
    for m in storelistsort[1:]:
        now=set(list(range(m[1],m[2]+1)))
        if len(set(havestore) & now)/len(now)>=0.9:
            continue
        else:
            havestore.extend(now)
            save.append(m)
    return save

def integrate_fragmentedmotif(tuplist,motiflen,tolerance):
	x=tuplist
	longest_key = max(x, key=lambda k: len(x[k]))
	initlongest_key=longest_key
	nowlen=len(x[longest_key])
	nowlist=[]
	nowlist.extend(x[longest_key])
	nowkey=initlongest_key
	while longest_key + 1 in x:
	    next_start = x[longest_key + 1][0][0]
	    last_end = x[nowkey][-1][1]
	    if (next_start - last_end)/motiflen < min(3,math.ceil(tolerance * nowlen)) and (next_start - last_end)>=0:
	        nowlen=len(nowlist)
	        nowlist.extend(x[longest_key + 1])
	        nowlist=sorted(nowlist)
	        longest_key += 1
	        nowkey=longest_key
	    elif (next_start - last_end)/motiflen >-0.2 and (next_start - last_end)<=0:
	        nowlen=len(nowlist)
	        nowlist.extend(x[longest_key + 1])
	        nowlist=sorted(nowlist)
	        longest_key += 1
	        nowkey=longest_key
	    else:
	        longest_key += 1
	longest_key=initlongest_key
	nowkey=initlongest_key
	while longest_key - 1 in x:
	    before_start = x[longest_key - 1][-1][1]
	    last_end = x[nowkey][0][0]
	    if (last_end  - before_start)/motiflen>-0.2 and (last_end  - before_start)/motiflen<=0:
	        nowlen=len(nowlist)
	        nowlen=nowlist[-1][1]-nowlist[0][0]
	        nowlist.extend(x[longest_key - 1])
	        nowlist=sorted(nowlist)
	        longest_key -= 1
	        nowkey=longest_key
	    elif (last_end  - before_start)/motiflen<min(3,math.ceil(tolerance * nowlen)) and (last_end  - before_start)/motiflen>=0:
	        nowlen=len(nowlist)
	        nowlen=nowlist[-1][1]-nowlist[0][0]
	        nowlist.extend(x[longest_key - 1])
	        nowlist=sorted(nowlist)
	        longest_key -= 1
	        nowkey=longest_key
	    else:
	        longest_key -= 1
	return nowlist






def tuple_process(sorted_tuples,motif_len):
	result_dict = defaultdict(list)
	id=0
	result_dict[id].append(sorted_tuples[0])
	for tupid in range(1,len(sorted_tuples)):
		if abs(sorted_tuples[tupid][0]-sorted_tuples[tupid-1][1])<=math.ceil(motif_len * 0.20):
			result_dict[id].append(sorted_tuples[tupid])
		else:
			id=id+1
			result_dict[id].append(sorted_tuples[tupid])
	return result_dict



def makefasta(input_file,allseq):
	with open(input_file, 'w') as f:
		f.write('>sequence\n')  # 描述行
		f.write(allseq+ '\n')
	trf_command = f'trf {input_file}  2 7 7 80 10 50 500 -h -l 20 -ngs'
	result = subprocess.run(['bash', '-c', trf_command], check=True, text=True, capture_output=True)
	result=result.stdout.split('\n')
	result=[item for item in result[1:] if item]
	return result


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
# 		print(sam)
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
# 				expand_motif(resdict,motif, aa, interval,0.2)
# 				if alllenlist!=[]:
# 					motiflist=[]
# 					enall=cal_cn1(alllenlist,len(motif),0.3)
# 					for mo in enall[3]:
# 					    motiflist.append(ref[mo[0]:(mo[1]+1)])
# 					cn=enall[0]
# 					if cn>cnmax:
# 						print(motif)
# 						print("a")
# 						cnmax=cn
# 						vntridcopy_number[sam]=cn
# 						enstatis=enall
# 						enmotiflist=motiflist
# 			t2tmotiflist.append([item,t2ts+enall[1]-1,t2ts+enall[2],cnmax,",".join(set(enmotiflist))])
# 			vntridcopy_number[sam]=cnmax
# 			if cnmax==0:
# 				vntridcopy_number[sam]="NONONO"
# 		else:
# 			print("aaaa")
# 			if len(ref)>100000:
# 				continue
# 			cnmax=0
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
# 						vntridcopy_number[sam]=cn1
# 				if alllenlist1==[] and alllenlist2!=[]:
# 					if cn2>cnmax:
# 						cnmax=cn2
# 						vntridcopy_number[sam]=cn2
# 			if cnmax==0:
# 				vntridcopy_number[sam]="NONONO"
# 			copy_number[vntrid]=vntridcopy_number
# 	t2tmotiflist_unique = [list(item) for item in set(tuple(sublist) for sublist in t2tmotiflist)]
# 	delnode=set(delnode)
# 	return copy_number,storep,delnode,t2tmotiflist_unique



