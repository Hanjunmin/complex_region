import pandas as pd
import argparse
import subprocess

parser = argparse.ArgumentParser(description='INTEGRATE THE TEINS')
parser.add_argument('--chr', help='chromosome')
args = parser.parse_args()


CHR=args.chr
sdpa=pd.read_csv("ref."+CHR+".anchoral.rm",sep="\t",header=None)
grouped = sdpa.groupby([10, 11, 12])
REF_SDR={}  ###生成SDR的dict
for name, group in grouped:
	#REF_SDR[name]=("SDR_"+group[6]).value_counts().to_dict()
	REF_SDR[name]=group.iloc[:,[0,1,2,3,4,5,6,7,8,9]]


sdpa=pd.read_csv("win.al.rm",sep="\t",header=None)
RMset=set(list(sdpa[9]))

grouped = sdpa.groupby([0, 1, 2])
REF_WIN={}  ###生成REF的dict
for name, group in grouped:
	REF_WIN[name]=("REF_"+group[9]).value_counts().to_dict()



import pandas as pd
from collections import Counter
samlis=pd.read_csv("sam.list",sep="\t",header=None)



# for sam in list(samlis[0]):
import multiprocessing as mp
def process_item(sam):
	print(sam)
	sdpa=pd.read_csv(sam+".win.outn",sep="\t",header=None)
	grouped = sdpa.groupby([0, 1, 2])
	alldict={}
	# x="a"
	for name, group in grouped:
		alldict[name]=("QUE_"+group[12]).value_counts().to_dict()

	###
	sdpa=pd.read_csv(sam+".win.outposn",sep="\t",header=None)
	grouped = sdpa.groupby([0, 1, 2])
	SDRalldict={}
	for name, group in grouped:
		# if name==("chrY",23983500,23997500):
		# 	break
		uniqlis=group.iloc[:, [3, 4, 5]].drop_duplicates()
		subdiclis=[]
		for  index,row in uniqlis.iterrows():
			row_tuple = (row[3], row[4], row[5])
			if row_tuple in REF_SDR.keys():
				subdiclis.extend(REF_SDR[row_tuple].values.tolist())
		subdiclisdf=pd.DataFrame(subdiclis).drop_duplicates()
		if len(subdiclisdf)!=0:
			SDRalldict[name]=("SDR_"+subdiclisdf[6]).value_counts().to_dict()
		# SDRdict= Counter()
		# # for d in subdiclis:
		# # 	SDRdict.update(d)
		# if SDRdict!=Counter():
		# 	SDRalldict[name]=SDRdict


	flattened_data = []
	win=pd.read_csv(sam+".win.txt",sep="\t",header=None)
	for _, (chr, start, end) in win[[0,1,2]].iterrows():
		id=(chr, start, end)
		row = {'chr': chr, 'start': start, 'end': end}
		if (chr, start, end) in REF_WIN.keys():  ####INS的dict
			row.update(REF_WIN[id])
		if (chr, start, end) in alldict.keys():  ####INS的dict
			row.update(alldict[(chr, start, end)])
		if (chr, start, end) in SDRalldict.keys():  ###SDR(REF)的dict
			row.update(SDRalldict[(chr, start, end)])
		flattened_data.append(row)
	input_bed=sam+".win.nointer.outpos"
	input_bed2=sam+".win.txt"
	output_file=sam+".temp"
	cmd = f"bedtools coverage -a {input_bed} -b {input_bed2} | awk '!($7 == 1 && $4==1)' |cut -f 1-3 > {output_file}"
	subprocess.run(cmd, shell=True, check=True)
	win2=pd.read_csv(sam+".temp",sep="\t",header=None)
	for _, (chr, start, end) in win2[[0,1,2]].iterrows():
		id=(chr, start, end)
		row = {'chr': chr, 'start': start, 'end': end}
		if (chr, start, end) in REF_WIN.keys():  ####INS的dict
			row.update(REF_WIN[id])
		if (chr, start, end) in alldict.keys():  ####INS的dict
			row.update(alldict[(chr, start, end)])
		if (chr, start, end) in SDRalldict.keys():  ###SDR(REF)的dict
			row.update(SDRalldict[(chr, start, end)])
		flattened_data.append(row)
	df = pd.DataFrame(flattened_data)
	df = df.fillna(0)
	# allcolumn=[]
	# for RMtype in RMset:
	# 	allcolumn.extend(["REF_"+RMtype,"SDR_"+RMtype,"QUE_"+RMtype])
	RMset={'Retroposon',  'LTR', 'DNA', 'SINE', 'RC',  'LINE'}
	allcolumn=[]
	for RMtype in RMset:
		allcolumn.extend(["REF_"+RMtype,"SDR_"+RMtype,"QUE_"+RMtype])
	for col in allcolumn:
		if col not in df.columns:
			df[col] = 0
	dfn=df.copy()
	for RMtype in RMset:
		dfn[RMtype]=dfn["REF_"+RMtype]-dfn["SDR_"+RMtype]+dfn["QUE_"+RMtype]
	dfn[sam]=dfn['Retroposon']+dfn['LTR']+dfn['DNA']+dfn['SINE']+dfn['RC']+dfn['LINE']
	import pybedtools
	dfn.loc[:,['chr','start','end',sam]].to_csv(sam+"temp.txt", sep="\t", index=False, header=False)
	b = pybedtools.BedTool(sam+"temp.txt")
	a = pybedtools.BedTool(CHR+".window")
	result = a.intersect(b, wa=True, wb=True,f=1)
	resultdf = result.to_dataframe()
	resultdf.columns=[0,1,2,3,4,5,"TE"]
	grouped =resultdf.groupby([3,4,5])
	# for name, group in grouped:
	#     if len(group)>1:
	#     	break
	group_stats = resultdf.groupby([3, 4, 5])["TE"].agg(["count", "first"]).reset_index()
	group_stats[sam]=group_stats['first']/group_stats['count']
	df = resultdf.merge(
		group_stats[[3, 4, 5, sam]],
		on=[3, 4, 5],
		how="left"
	)
	df=df.loc[:,[0,1,2,sam]]
	df.columns=['chr','start','end',sam]
	df.to_csv(sam+"TE.txt", sep="\t", index=False, header=True)
	return "OK"

with mp.Pool(processes=50) as pool:
    results = pool.map(process_item, list(samlis[0]))

# ###CHM13
# flattened_data = []
# win=pd.read_csv(CHR+".window",sep="\t",header=None)
# for _, (chr, start, end) in win[[0,1,2]].iterrows():
# 	id=(chr, start, end)
# 	row = {'chr': chr, 'start': start, 'end': end}
# 	if (chr, start, end) in REF_WIN.keys():  ####INS的dict
# 		row.update(REF_WIN[id])
# 	flattened_data.append(row)

# df = pd.DataFrame(flattened_data)
# df = df.fillna(0)
# # allcolumn=[]
# # for RMtype in RMset:
# # 	allcolumn.extend(["REF_"+RMtype,"SDR_"+RMtype,"QUE_"+RMtype])
# RMset={'Retroposon',  'LTR', 'DNA', 'SINE', 'RC',  'LINE'}
# allcolumn=[]

# for RMtype in RMset:
# 	allcolumn.extend(["REF_"+RMtype])

# for col in allcolumn:
# 	if col not in df.columns:
# 		df[col] = 0
# dfn=df.copy()
# for RMtype in RMset:
# 	dfn[RMtype]=dfn["REF_"+RMtype]

# dfn["REF"]=dfn['Retroposon']+dfn['LTR']+dfn['DNA']+dfn['SINE']+dfn['RC']+dfn['LINE']
	

# alwin=pd.read_csv(CHR+".window",sep="\t",header=None)
# samlis=pd.read_csv("sam.list",sep="\t",header=None)
# alwin.columns=['chr','start','end']

# for sam in list(samlis[0]):
# 	print(sam)
# 	# if sam=="C002-CHA-E02-Mat":
# 	# 	break
# 	samdf=pd.read_csv(sam+"TE.txt",sep="\t",header=0)
# 	alwin=alwin.merge(samdf.drop_duplicates(),on=["chr", "start","end"],how="left")
# 	print(len(alwin))

# alwin=alwin.merge(dfn[['chr','start','end','REF']].drop_duplicates(),on=["chr", "start","end"],how="left")
	
# alwin.to_csv("TEal.txt", sep="\t", index=False, header=True)