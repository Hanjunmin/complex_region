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
from functools import reduce

###################################



def inte_rowsname(end,samndic): ###这个函数是把column进行合并
	merged_columns = []
	temp_start, temp_end = None, None
	temp_base = None
	allcom=end.columns[:-7]
	data = []
	for col in allcom:
		base, rows = col.split('_rows_')
		start, end_row = map(int, rows.split('_to_'))
		data.append((base, start, end_row,col))
	dfdata = pd.DataFrame(data, columns=['base', 'start', 'end','colname']).sort_values(['base', 'start'])
	merged_columns = []
	for base, group in dfdata.groupby('base'):
		group = group.sort_values('start').reset_index(drop=True)
		group['clus']=0
		clus=0
		for idx, row in group.iloc[1:].iterrows():
			if row['start'] == group.iloc[idx - 1]['end'] + 1:
				group.loc[idx, 'clus'] = group.loc[idx - 1, 'clus']
			else:
				clus = clus + 1
				group.loc[idx, 'clus'] = clus
		for clus_id, subgroup in group.groupby('clus'):
			base_val = base
			min_start = subgroup['start'].min()
			max_end = subgroup['end'].max()
			colname_list = subgroup['colname'].tolist()
			merged_columns.append((base_val, min_start, max_end, colname_list))
	for mer in merged_columns:
		if len(mer[3])!=1:
			nname=mer[0]+"_rows_"+str(mer[1])+"_to_"+str(mer[2])
			end[nname] = end[mer[3]].apply(
			lambda row: [x for x in row if pd.notna(x)], axis=1
			)
			for idx, row in enumerate(end.iterrows(), start=0): ###对合并后的列进行处理
				if row[1][nname]!=[]:
					if len(row[1][nname])==1:
						end[nname].values[idx] = row[1][nname][0]
					else:
						print("test")
						end[nname].values[idx]=reduce(lambda x, y: x + y, row[1][nname])
			end[nname] = end[nname].apply(lambda x: None if x == [] else x)	
			if nname in end.columns:  # 检查列是否存在
				column = end.pop(nname)  # 将 'nname' 列移除
				# end.drop(nname, axis=1, inplace=True)  # 删除已存在的列
				end.insert(0, nname, column)  # 将列插入到最前面
				end =end.drop(mer[3], axis=1)  # 删除名为 'col1' 和 'col2' 的两列
				df_cleaned =end[end[nname].notna()]
				samndic[nname]=pd.Series(df_cleaned[nname], index=df_cleaned.index).to_dict()
	return end,samndic

    



def process_line2(line,uniqtup,samdic):
	line_list=line.strip().split('\t',)
	df = pd.DataFrame(columns=['SDregion', 'ps', 'pe',"pc","sim"])
	df2 = pd.DataFrame(columns=['SDregion', 'Sample', 'Sim'])
	my_list=[]
	line_list=line.strip().split('\t',)
	letters_to_find =[part for part in re.split(r'[<>]', line_list[3]) if part]
	rilis=re.findall(r'[<>]', line_list[3])
	samples_with_intersection = list(uniqtup)
	positions = {}
	
	with envpos.begin(write=False) as txn:
		for sample in samples_with_intersection:
			positions[sample] = list(map(samdic[sample].get, letters_to_find))
	
	end=pd.DataFrame(positions)
	end.index=letters_to_find
	with env.begin(write=False) as txn:
		cursor = txn.cursor()
		values = cursor.getmulti([key.encode('utf-8') for key in letters_to_find])
		end['fa'] = [msgpack.unpackb(value[1]) for value in values if value[1] is not None]
	
	end['fa_length'] = end['fa'].apply(len)
	end['start_pos'] = end['fa_length'].cumsum() - end['fa_length'] + 1
	end['end_pos'] = end['fa_length'].cumsum()
	end['ori']=rilis
	start_range=int(line_list[1])
	end_range=int(line_list[2])
	end['marked'] = np.where(
	(end['start_pos'] < start_range) & (end['end_pos'] > end_range), 'Yes3',  # 同时满足两个条件
	np.where(end['start_pos'] < start_range, 'Yes1',  # 仅小于 start_range
	np.where(end['end_pos'] > end_range, 'Yes2', 'No'))  # 仅大于 end_range
	)
	end['newfa'] = end.apply(
		lambda row: row['fa'][(start_range):min(end_range,row['fa_length'])]
		if row['start_pos'] < start_range and row['end_pos'] <= end_range else row['fa'] ,
		axis=1
	)
	end['newfa'] = end.apply(
		lambda row: row['fa'][:(end_range - row['start_pos'] + 1)]
		if row['end_pos'] > end_range and row['start_pos'] >= start_range else row['newfa'] ,
		axis=1
	)
	end['newfa'] = end.apply(
		lambda row: row['fa'][(start_range):min(end_range,row['fa_length'])]
		if row['start_pos'] < start_range and row['end_pos'] > end_range else row['newfa'] ,
		axis=1
	)
	end,samndic=inte_rowsname(end,{})
	merged_dict = samdic | samndic
	for col in end.columns[:-7]:
		first_two_columnsdf = pd.DataFrame(columns=['Column1', 'Column2'])
		rows_without_na = end[(~end[col].isnull())]
		rows_with_na = end[(end[col].isnull())]
		length_sum = rows_with_na['newfa'].apply(len).sum()
		init_len=int(line_list[2])-int(line_list[1])+1
		new_row_df2=[line_list[0],col,1-(length_sum/init_len)]
		new_row_df2 = pd.DataFrame([new_row_df2], columns=df2.columns)
		df2 = pd.concat([df2,new_row_df2], ignore_index=True)
		if rows_without_na.shape[0]!=0:
			rows_without_na_no = end[(~end[col].isnull()) & (end['marked'] == 'No')]
			if rows_without_na_no.shape[0]!=0:
				indices = rows_without_na_no.index.tolist()
				mapped_values =pd.Series(indices).apply(lambda x: merged_dict[col].get(x))
				multiple_lists_count = (mapped_values.apply(len) > 1).sum()  ##这些大于1的说明有多个node，就需要重新去算
				if multiple_lists_count>1:
					print("no!!!!!")
					print(col)
					my_list.extend([line_list[0],0,col])
					odf=multi(line_list[0],end,col,start_range,end_range,sam,merged_dict)
					print(odf)
					df = pd.concat([df, odf], ignore_index=True)
					continue
				all_lists = [item for sublist in mapped_values for item in sublist]
				all_lists_array = np.array(all_lists)
				first_two_columns = all_lists_array[:, :2]
				first_two_columnsdf=pd.DataFrame(first_two_columns)
			rows_without_naother1=end[(~end[col].isnull()) & (end['marked'] == 'Yes1')]
			rows_without_naother2=end[(~end[col].isnull()) & (end['marked'] == 'Yes2')]
			rows_without_naother3=end[(~end[col].isnull()) & (end['marked'] == 'Yes3')]
			xx1=rows_without_naother1.loc[:,[col,"fa","newfa","fa_length","start_pos","end_pos","ori","marked"]]
			xx2=rows_without_naother2.loc[:,[col,"fa","newfa","fa_length","start_pos","end_pos","ori","marked"]]
			xx3=rows_without_naother3.loc[:,[col,"fa","newfa","fa_length","start_pos","end_pos","ori","marked"]]
			for index, row in xx1.iterrows():
				pose=(min(end_range,row['fa_length']))
				poss=(start_range)
				befores=merged_dict[col][index][0][0]
				beforee=merged_dict[col][index][0][1]
				if merged_dict[col][index][0][4]==row['ori']:
					afts=befores+poss
					afte=befores+pose
				else:
					aftsrev=row['fa_length']-pose
					afterev=row['fa_length']-poss
					afts=befores+aftsrev
					afte=befores+afterev-1
				new_row = pd.DataFrame([[afts, afte]], columns=first_two_columnsdf.columns)
				first_two_columnsdf = pd.concat([first_two_columnsdf, new_row], ignore_index=True)
			for index, row in xx2.iterrows():
				poss=0
				pose=end_range - row['start_pos'] + 1
				befores=merged_dict[col][index][0][0]
				beforee=merged_dict[col][index][0][1]
				if merged_dict[col][index][0][4]==row['ori']:
					afts=befores+poss
					afte=befores+pose
				else:
					aftsrev=row['fa_length']-pose
					afterev=row['fa_length']-poss
					afts=befores+aftsrev
					afte=befores+afterev-1
				new_row = pd.DataFrame([[afts, afte]], columns=first_two_columnsdf.columns)
				first_two_columnsdf = pd.concat([first_two_columnsdf, new_row], ignore_index=True)
			for index, row in xx3.iterrows():
				if df.shape[0]!=0:
					poss=start_range
					pose=end_range
					befores=merged_dict[col][index][0][0]
					beforee=merged_dict[col][index][0][1]
					if merged_dict[col][index][0][4]==row['ori']:
						afts=befores+poss
						afte=befores+pose
					else:
						aftsrev=row['fa_length']-pose
						afterev=row['fa_length']-poss
						afts=befores+aftsrev
						afte=befores+afterev-1
					new_row = pd.DataFrame([[afts, afte]], columns=first_two_columnsdf.columns)
					first_two_columnsdf = pd.concat([first_two_columnsdf, new_row], ignore_index=True)
				if df.shape[0]==0 and end.shape[0]==1:
					poss=start_range
					pose=end_range
					befores=merged_dict[col][index][0][0]
					beforee=merged_dict[col][index][0][1]
					if merged_dict[col][index][0][4]==row['ori']:
						afts=befores+poss
						afte=befores+pose
					else:
						aftsrev=row['fa_length']-pose
						afterev=row['fa_length']-poss
						afts=befores+aftsrev
						afte=befores+afterev-1
					new_row = pd.DataFrame([[afts, afte]], columns=first_two_columnsdf.columns)
					first_two_columnsdf = pd.concat([first_two_columnsdf, new_row], ignore_index=True)
			first_two_columnsdf.iloc[:, 0] = pd.to_numeric(first_two_columnsdf.iloc[:, 0], errors='coerce')
			first_two_columnsdf.iloc[:, 1] = pd.to_numeric(first_two_columnsdf.iloc[:, 1], errors='coerce')
			ps=first_two_columnsdf.iloc[:,0].min()
			pe=first_two_columnsdf.iloc[:,1].max()
			pc=col
			ratio=(min(pe-ps,end_range-start_range))/(max(pe-ps,end_range-start_range))
			filt=True
			if ratio<=0.5:
				# dbscan = DBSCAN(eps=(end_range-start_range), min_samples=1)
				if first_two_columnsdf.shape[0]==1:
					filt=True
				else:
					merout=bedmerge(first_two_columnsdf.iloc[:,0:2],line_list[0],sam)
		 			#label=dbscan.fit_predict(first_two_columnsdf.iloc[:,0:2])
					if len(merout)!=1:
						clusterer = hdbscan.HDBSCAN(min_cluster_size=2, min_samples=1, cluster_selection_epsilon=(end_range-start_range)/2,core_dist_n_jobs=1)
						first_two_columnsdf=bedmerge2(merout,clusterer,first_two_columnsdf) #label=clusterer.fit_predict(first_two_columnsdf.iloc[:,0:2])
						label=first_two_columnsdf['cluster']
						for lab in set(label):
							x=first_two_columnsdf[first_two_columnsdf['cluster']==lab]
							ps1=x.iloc[:,0].min()
							pe1=x.iloc[:,1].max()
							ratio=(min(pe-ps,end_range-start_range))/(max(pe-ps,end_range-start_range))
							if ratio>0.5:
								filt=False
								new_row_df=[line_list[0],ps1,pe1,pc,ratio]
								new_row_df = pd.DataFrame([new_row_df], columns=df.columns)
								df = pd.concat([df, new_row_df], ignore_index=True)
					else:
						filt=True
			if filt:
				new_row_df=[line_list[0],ps,pe,pc,ratio]
				new_row_df = pd.DataFrame([new_row_df], columns=df.columns)
				df = pd.concat([df, new_row_df], ignore_index=True)
	return df, my_list,df2





def get_nested_value_from_lmdb(db_path,key):
    with db_path.begin() as txn:
        value = msgpack.unpackb(txn.get(key.encode('utf-8')),raw=False, strict_map_key=False,use_list=False)  # 外层键
        return value
    

	

def bedmerge(mergebef,possd,sam):
	mergebef = mergebef.copy()
	mergebef.columns=["start","end"]
	mergebef['chr']="chrN"
	new_df = mergebef[['chr', 'start', 'end']]
	input_file=possd+sam+".txt"
	new_df.to_csv(input_file, sep='\t', index=False,header=False)
	result = subprocess.run(
		['bash', '-c', f'/home/Public/software/bedtools-2.30.0/bedtools sort -i {input_file} | /home/Public/software/bedtools-2.30.0/bedtools merge'],  # 使用管道
		stdout=subprocess.PIPE, 
		stderr=subprocess.PIPE,
		text=True
	)
	os.remove(input_file)
	data_io = StringIO(result.stdout)
	merout = pd.read_csv(data_io, sep='\t', header=None, names=['chr', 'start', 'end'])
	return merout

def bedmerge2(merout,clusterer,mergebef):
	label=clusterer.fit_predict(merout.iloc[:,1:3])
	merout['label']=label
	for lab in set(label):
		cho=merout[merout['label']==lab]
		rana=list(cho['start'])[0]
		ranb=list(cho['end'])[-1]
		mergebef.loc[(mergebef.iloc[:, 0] >= rana) & (mergebef.iloc[:, 1] <= ranb), 'cluster'] = lab
	return mergebef


def multi(possd,end,col,start_range,end_range,sam,merdic):
	print(possd)
	print(end)
	mul = pd.DataFrame(columns=['SDregion', 'ps', 'pe',"pc","sim"])
	mulall=pd.DataFrame(columns=['SDregion', 'ps', 'pe',"pc","sim"])
	rows_without_na = end[(~end[col].isnull())]
	if rows_without_na.shape[0]!=0:
		indices = rows_without_na.index.tolist()
		mapped_values =pd.Series(indices).apply(lambda x: merdic[col].get(x))
		all_lists = [list(item) + [index] for index, sublist in zip(indices, mapped_values) for item in sublist]
		all_lists_array = np.array(all_lists)
		first_two_columns = all_lists_array[:, :6]
		first_two_columnsdf=pd.DataFrame(first_two_columns)
		first_two_columnsdf.iloc[:, 0] = pd.to_numeric(first_two_columnsdf.iloc[:, 0], errors='coerce')
		first_two_columnsdf.iloc[:, 1] = pd.to_numeric(first_two_columnsdf.iloc[:, 1], errors='coerce')
		merout=bedmerge(first_two_columnsdf.iloc[:,0:2],possd,sam)
		if len(merout)!=1:
			dbscan = DBSCAN(eps=(end_range-start_range), min_samples=1)
			first_two_columnsdf=bedmerge2(merout,dbscan,first_two_columnsdf) #label=dbscan.fit_predict(first_two_columnsdf.iloc[:,0:2])
		else:
			first_two_columnsdf['cluster']=0
		label=first_two_columnsdf['cluster']
	filtall=True
	for clus in set(label):
		end1id=first_two_columnsdf[first_two_columnsdf['cluster']==clus]
		end1=end.loc[end1id.iloc[:,5]]
		filtered_end1 = end1[end1['marked'] == "No"]
		valid_indices = filtered_end1.index.values
		filtered_df = end1id[end1id.iloc[:, 5].isin(valid_indices)]
		result2 = filtered_df.iloc[:, :2]
		newdic= defaultdict(list)
		end1id.columns=['start', 'end', 'chr_pos', 'chr_val','chr_ori','node','cluster']
		for start, endno, chr_pos, chr_val,chr_ori,node,cluster  in zip(end1id['start'],end1id['end'],end1id['chr_pos'],end1id['chr_val'],end1id['chr_ori'],end1id['node'],end1id['cluster']):
			newdic[node].append([start, endno, chr_pos, chr_val,chr_ori])	
		rows_without_naother1=end1[(~end1[col].isnull()) & (end1['marked'] == 'Yes1')]
		rows_without_naother2=end1[(~end1[col].isnull()) & (end1['marked'] == 'Yes2')]
		rows_without_naother3=end1[(~end1[col].isnull()) & (end1['marked'] == 'Yes3')]
		xx1=rows_without_naother1.loc[:,[col,"fa","newfa","fa_length","start_pos","end_pos","ori","marked"]]
		xx2=rows_without_naother2.loc[:,[col,"fa","newfa","fa_length","start_pos","end_pos","ori","marked"]]
		xx3=rows_without_naother3.loc[:,[col,"fa","newfa","fa_length","start_pos","end_pos","ori","marked"]]
		for index, row in xx1.iterrows():
			pose=(min(end_range,row['fa_length']))
			poss=(start_range)
			befores=int(newdic[index][0][0])
			beforee=int(newdic[index][0][1])
			if newdic[index][0][4]==row['ori']:
				afts=befores+poss
				afte=befores+pose
			else:
				aftsrev=row['fa_length']-pose
				afterev=row['fa_length']-poss
				afts=befores+aftsrev
				afte=befores+afterev-1
			new_row = pd.DataFrame([[afts, afte]], columns=result2.columns)
			result2 = pd.concat([result2, new_row], ignore_index=True)
		for index, row in xx2.iterrows():
			poss=0
			pose=end_range - row['start_pos'] + 1
			befores=int(newdic[index][0][0])
			beforee=int(newdic[index][0][1])
			if newdic[index][0][4]==row['ori']:
				afts=befores+poss
				afte=befores+pose
			else:
				aftsrev=row['fa_length']-pose
				afterev=row['fa_length']-poss
				afts=befores+aftsrev
				afte=befores+afterev-1
			new_row = pd.DataFrame([[afts, afte]], columns=result2.columns)
			result2 = pd.concat([result2, new_row], ignore_index=True)
		for index, row in xx3.iterrows():
			if xx3.shape[0]!=0:
				poss=start_range
				pose=end_range
				befores=int(newdic[index][0][0])
				beforee=int(newdic[index][0][1])
				if newdic[index][0][4]==row['ori']:
					afts=befores+poss
					afte=befores+pose
				else:
					aftsrev=row['fa_length']-pose
					afterev=row['fa_length']-poss
					afts=befores+aftsrev
					afte=befores+afterev-1
				new_row = pd.DataFrame([[afts, afte]], columns=result2.columns)
				result2 = pd.concat([result2, new_row], ignore_index=True)
		result2.iloc[:, 0] = pd.to_numeric(result2.iloc[:, 0], errors='coerce')
		result2.iloc[:, 1] = pd.to_numeric(result2.iloc[:, 1], errors='coerce')
		ps=result2.iloc[:,0].min()
		pe=result2.iloc[:,1].max()
		pc=col
		ratio=(min(pe-ps,end_range-start_range))/(max(pe-ps,end_range-start_range))
		filt=True
		if ratio<=0.5:
			if first_two_columnsdf.shape[0]==1:
				filt=True
			else:
				merout=bedmerge(result2.iloc[:,0:2],possd,sam)
				if len(merout)!=1:
					dbscan = DBSCAN(eps=(end_range-start_range), min_samples=1)
					result2=bedmerge2(merout,dbscan,result2) #label=clusterer.fit_predict(first_two_columnsdf.iloc[:,0:2])
					label1=result2['cluster']
				else:
					result2['cluster']=0
					label1=[0]
				for lab in set(label1):
					x=result2[result2['cluster']==lab]
					ps1=x.iloc[:,0].min()
					pe1=x.iloc[:,1].max()
					ratio=(min(pe1-ps1,end_range-start_range))/(max(pe1-ps1,end_range-start_range))
					if ratio>0.5:
						filt=False
						filtall=False
						new_row_df=[possd,ps1,pe1,pc,ratio]
						new_row_df = pd.DataFrame([new_row_df], columns=mul.columns)
						mul = pd.concat([mul, new_row_df], ignore_index=True)
		else:
			filt=False
			filtall=False
			new_row_df=[possd,ps,pe,pc,ratio]
			new_row_df = pd.DataFrame([new_row_df], columns=mul.columns)
			mul = pd.concat([mul, new_row_df], ignore_index=True)
		if filt:
			new_row_df=[possd+"_complex",ps,pe,pc,ratio]
			new_row_df = pd.DataFrame([new_row_df], columns=mul.columns)
			mulall = pd.concat([mulall, new_row_df], ignore_index=True)
	if filtall:
		mul= pd.concat([mul, mulall], ignore_index=True)
	return(mul)


parser = argparse.ArgumentParser(description='Pan-SD-idencordinate')
parser.add_argument('--W', help='graph path')
parser.add_argument('--S', help='corresponding sample name')
parser.add_argument('--i', help='transform region')
parser.add_argument('--o', help='output file')
parser.add_argument('--p', help='work folder absolute path')
args = parser.parse_args()


sam=args.S
win=args.W
inputpos=args.i
outputfile=args.o
outputfile1=outputfile+".sim"
workdir=args.p




# sam="C002-CHA-E02-Mat"
# workdir="/home/jmhan/iterall"
envfinaldic = lmdb.open(workdir+'/db/'+sam+'/'+'finaldb', create=False,readonly=True, lock=False)  # 请根据需求修改路径和大小
envpos = lmdb.open(workdir+'/db/'+sam+'/'+'segpos', create=False,readonly=True, lock=False)  # 10 GB
envori = lmdb.open(workdir+'/db/'+sam+'/'+'segori', create=False,readonly=True, lock=False)  # 10 GB
env = lmdb.open(workdir+'/seg_db', create=False,readonly=True, lock=False)  # 10 GB
    
##分个组
with open(win, 'r') as f:
	lines = f.readlines()


with envpos.begin() as txn:
    cursor = txn.cursor()  # 获取数据库的游标
    nodedic = {}  # 用来存储所有解包后的字典数据
    for key, value in cursor:
        unpacked_data = msgpack.unpackb(value)
        nodedic[key.decode('utf-8')] = unpacked_data.keys()  # 将解包后的字典存储在 all_data 中，使用键作为字典的键


path_sets = {sample: set(path) for sample, path in nodedic.items()}


#bedtools  subtract -b /home/jmhan/pan-SD/acrocentric.region  -a <(less -S alldfe.txt |tr ":" "\t" |tr "-" "\t" |less -S) -A |less -S |awk '{print $1":"$2"-"$3"\t"$4"\t"$5"\t"$6}' |less -S >tessd.bed
#bedtools  subtract -b acro.region  -a /dssg/home/acct-clsmyf/clsmyf-user1/pangenome/APG_seg/input.rev -A  >input.rev.delacro

#"/home/jmhan/pangenome/WGSfea/allW/reproce.input"
# inputpos="/home/jmhan/pangenome/WGSfea/allW/APG3allW/x"
with open(inputpos, 'r') as f:
	lines = f.readlines()


indexed_lines = [f"{i}\t{line}" for i, line in enumerate(lines)]


# for line in indexed_lines:
def pro_intersect(chun):
	reschun=[]
	for line in chun:
		line_list=line.strip().split('\t',)
		my_list=[]
		line_list=line.strip().split('\t',)
		print(line_list[1]) #line=lines[38762]
		letters_to_find =[part for part in re.split(r'[<>]', line_list[4]) if part]
		rilis=re.findall(r'[<>]', line_list[4])
		setletters_to_find=set(letters_to_find)
		samples_with_intersection = [sample for sample, path in path_sets.items() if setletters_to_find & path]
		reschun.append([int(line_list[0]),samples_with_intersection])
	# empty_df = pd.DataFrame()	
	# empty_df.to_csv("/home/jmhan/iterall/test/"+str(line_list[0]), sep='\t', index=False)
	return reschun


import itertools

def chunked_data(iterable, chunk_size):
    it = iter(iterable)
    for chunk in iter(lambda: list(itertools.islice(it, chunk_size)), []):
        yield chunk

chunk_size=10000
chunksline = list(chunked_data(indexed_lines, chunk_size))  # 创建数据块
with Pool(processes=15) as pool:  # 设置并行核心数
	results = pool.map(pro_intersect, chunksline)


test=[]
for res in results:
	for ree in res:
		test.append(ree)



# with open("length_dict.json", "r") as f:
#     fastalen_dict = next(ijson.items(f, '', use_float=True))  # 从根节点读取整个 JSON 对象


testdf=pd.DataFrame(test)
tuple_series = (testdf)[1].apply(tuple)
unique_items_count = tuple_series.value_counts()
import pickle
with open(sam+'unique_items_count.pkl', 'wb') as f:
    pickle.dump(unique_items_count, f)

unique_items = set(tuple_series)
testdf['tuple'] = tuple_series
del test
del tuple_series


processid = 1
total = len(unique_items)
all_results = []  # List to collect results from all uniqtups
items_greater_than_100 = set(unique_items_count[unique_items_count > 300].index)
items_smaller_than_100 = set(unique_items_count[unique_items_count <= 50].index)
filtered_set = set(unique_items_count[(unique_items_count > 50) & (unique_items_count <= 300)].index)

print("HELLO!!!!")
print(len(items_greater_than_100))
print("HELLO!!!!")
print(len(items_smaller_than_100))
def process_chunk(uniqtup):
    print(uniqtup)
    global processid
    print(f"Processing##########################################: {processid}/{total}")  # 打印当前进度
    processid = processid + 1  # 更新进度
    if uniqtup == ():
        return []  # 如果当前元组为空，返回空列表
    uniqtupdf = testdf[testdf['tuple'] == uniqtup]
    sel = list(uniqtupdf[0])  # 获取需要处理的索引
    samdic = {}
    for tupfac in uniqtup:
        samdic[tupfac] = get_nested_value_from_lmdb(envfinaldic, tupfac)
    results=[]
    for line in [lines[i] for i in sel]:
        resnow=process_line2(line,uniqtup,samdic)
        results.append(resnow)
    # empty_df = pd.DataFrame()
    # random_code = ''.join(random.choices(string.ascii_letters + string.digits, k=20))
    # empty_df.to_csv("/home/jmhan/iterall/test/"+str(random_code), sep='\t', index=False)
    return results



with Pool(processes=20) as pool:  # 主进程中的并行
    all_results = pool.map(process_chunk, list(items_smaller_than_100))


# my_list=[]
df = pd.DataFrame(columns=['SDregion', 'ps', 'pe',"pc"])
dff =  pd.DataFrame(columns=['SDregion', 'Sample', 'Sim'])
for res in all_results:
	for df1,my_list1,df2 in res:
		# my_list.append(my_list1)
		df = pd.concat([df, df1], ignore_index=True)
		dff = pd.concat([dff, df2], ignore_index=True)


with Pool(processes=20) as pool:  # 主进程中的并行
    all_results = pool.map(process_chunk, list(filtered_set))

for res in all_results:
	for df1,my_list1,df2 in res:
		# my_list.append(my_list1)
		df = pd.concat([df, df1], ignore_index=True)
		dff = pd.concat([dff, df2], ignore_index=True)

import functools
all_results = []  # List to collect results from all uniqtups
for uniqtup in items_greater_than_100:
    print(f"Processing##########################################: xx/{total}")  # 打印当前进度
    if uniqtup == ():
        continue
    uniqtupdf = testdf[testdf['tuple'] == uniqtup]
    sel=list(uniqtupdf[0])
    samdic={}
    for tupfac in uniqtup:
        samdic[tupfac]=get_nested_value_from_lmdb(envfinaldic,tupfac)
    partial_process_line2 = functools.partial(process_line2, uniqtup=uniqtup, samdic=samdic)
    with Pool(processes=15) as pool:  # 设置并行核心数
        results = pool.map(partial_process_line2, [lines[i] for i in sel])
    all_results.extend(results)
        # processid=processid+1
        # empty_df = pd.DataFrame()
        # random_code = ''.join(random.choices(string.ascii_letters + string.digits, k=20))
        # empty_df.to_csv("/home/jmhan/iterall/test/"+str(random_code), sep='\t', index=False)

# time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())



###注意一下这里的ps 和pe指的是原先的path（最初始的W那个path）对应下的坐标 所以只看最开始的就好啦
#### eg:merged_dict[X_chr1_157181305_224264158_rows_3580432_to_5352659']['4190468']
#### ((3580432, 3580432, 'C002-CHA-E02-Mat_0_C002-CHA-E02_Mat_chr1_157181305_224264158_3580432', 'C002-CHA-E02-Mat_0_C002-CHA-E02_Mat_chr1_157181305_224264158', '>'),)
#### merged_dict['C002-CHA-E02-Mat_0_C002-CHA-E02_Mat_chr1_157181305_224264158_rows_3580432_to_5352659']['4267147']
#### ((5352590, 5352659, 'C002-CHA-E02-Mat_0_C002-CHA-E02_Mat_chr1_157181305_224264158_5352590', 'C002-CHA-E02-Mat_0_C002-CHA-E02_Mat_chr1_157181305_224264158', '>'),)

for df1,my_list1,df2 in all_results:
    df = pd.concat([df, df1], ignore_index=True)
    dff = pd.concat([dff, df2], ignore_index=True)

#"can_sd.corregion.bed"
df.to_csv(outputfile, sep='\t', index=False)
dff.to_csv(outputfile1, sep='\t', index=False)
# my_list = [item for item in my_list if item != []]




