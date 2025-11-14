import argparse
import csv
import re
import numpy as np
import pandas as pd
import json
import os
from concurrent.futures import ProcessPoolExecutor
import msgpack
import lmdb

parser = argparse.ArgumentParser(description='Pan-SD')
parser.add_argument('--W', help='sample W')
parser.add_argument('--d', help='The path of database')

# parser.add_argument('--mercan', help='ther merged result of sdpair')

# parser.add_argument('--sdpath', help='sd path')
args = parser.parse_args()


# with open('/home/jmhan/pan-SD/CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
#     fasta_dict = json.load(gz_file)

# env = lmdb.open('/share/home/zhanglab/user/maoyafei/project/pansd/db/CHM13-APGp1-HPRCp1-HGSVCp3_MC_db_msgpack', create=False,readonly=True, lock=False) 
env = lmdb.open(args.d, create=False,readonly=True, lock=False) 


Wpath =args.W
# chooseW =args.mercan
# Wpath ="/share/home/zhanglab/user/maoyafei/decompose/splinew/C001-CHA-E01-Mat.W"
chooseW ="process.W"
mercan="../uniqpair.sort"



# Wpath ="/home/jmhan/pangenome/WGSfea/allW/APG3allWnew/C001-CHA-E01-Mat.W"
# chooseW ="../process.W"
dfch = pd.read_csv(chooseW,sep="\t",header=None)
dfch2 = pd.read_csv(mercan,sep="\t",header=None)
with open(Wpath, 'r') as f:
	lines = f.readlines()


# for line in lines:
# def process_line(line):
for line in lines:
	line_list = line.strip().split('\t')
	#print(line_list[0])
	intersdpair=((dfch2[0] == line_list[3]) & (dfch2[1] <= int(line_list[5])) & (dfch2[2] >= int(line_list[4])))
	matching_rows = dfch2[intersdpair]
	if line_list[3] in dfch[0].values and int(line_list[4]) in dfch[1].values and int(line_list[5]) in dfch[2].values:
		print(line_list[3])
		letters_to_find =[part for part in re.split(r'[<>]', line_list[6]) if part]
		orilis=re.findall(r'[<>]', line_list[6])
		end = pd.DataFrame()
		end['node']=letters_to_find 
		end['ori']=orilis
		with env.begin(write=False) as txn:
			cursor = txn.cursor()
			values = cursor.getmulti([key.encode('utf-8') for key in letters_to_find])
			end['fa'] = [msgpack.unpackb(value[1]) for value in values if value[1] is not None]
		# end['fa']=list(map(fasta_dict.get, letters_to_find))  #end.isnull().values.any()
		end['start_position'] = end['fa'].apply(len).cumsum() - end['fa'].apply(len) + 1
		end['end_position'] = end['start_position'] + end['fa'].str.len()-1
		end['chr_position'] = line_list[1]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5] +"@"+ end['start_position'].astype(str)
		end['chr']=line_list[3]
		end['pathstart']=int(line_list[4])
		end['genome_start']=int(line_list[4])+end['start_position']
		end['genome_end']=int(line_list[4])+end['end_position'] #### chr genomestart genome_end
		genome_start_array = end['genome_start'].values
		genome_end_array = end['genome_end'].values
		new_df = pd.DataFrame(columns=end.columns)  # 初始化新数据框，列名与 end 数据框相同
		for index, row in matching_rows.iterrows():
			sdpairscan=row[1]
			sdpairecan=row[2]
			if len(np.where(genome_start_array <= sdpairscan)[0])==0:
				continue
			if len(np.where(genome_end_array >= sdpairecan)[0])==0:
				continue
			index_start = np.where(genome_start_array <= sdpairscan)[0][-1]
			index_end = np.where(genome_end_array >= sdpairecan)[0][0]
			splicedf=end.loc[index_start:index_end]
			new_df = pd.concat([new_df, splicedf], ignore_index=True)
		print(new_df)
		new_df.to_csv(line_list[1]+line_list[3]+"_"+line_list[4]+"_"+line_list[5]+'chooWpos.tsv', sep='\t', index=False)
# 		return new_df
# 	else:
# 		return pd.DataFrame()

# with ProcessPoolExecutor(max_workers=30) as executor:
#     results = executor.map(process_line, lines)




# ###下面为所有样本一起跑的并行：

# with open("/home/jmhan/pangenome/WGSfea/allW/APG3allWnew/allsamtest", 'r') as file:
# 	for sam in file:
# 		sam=sam.strip() # 使用 strip() 去除行末的换行符
# 		Wpath="/home/jmhan/pangenome/WGSfea/allW/APG3allWnew/"+sam+".W"
# 		chooseW="/home/jmhan/iterall/uniqpair/"+sam+"/process.W"
# 		mercan="/home/jmhan/iterall/uniqpair/"+sam+"/uniqpair.mer"
# 		if os.path.exists(chooseW):
# 			dfch = pd.read_csv(chooseW,sep="\t",header=None)
# 			dfch2 = pd.read_csv(mercan,sep="\t",header=None)
# 			if os.path.exists(chooseW):
# 				with open(Wpath, 'r') as f:
# 					lines = f.readlines()
# 				os.chdir("/home/jmhan/iterall/uniqpair/"+sam+"/samPathpos/")
# 				print("当前工作目录:", os.getcwd())
# 				with ProcessPoolExecutor(max_workers=30) as executor:
# 					results = executor.map(process_line, lines)






# Wpath ="/home/jmhan/pangenome/WGSfea/allW/APG3allWnew/C001-CHA-E01-Mat.W"
# chooseW ="process.W"




# with ProcessPoolExecutor(max_workers=30) as executor:
#     results = executor.map(process_line, lines)