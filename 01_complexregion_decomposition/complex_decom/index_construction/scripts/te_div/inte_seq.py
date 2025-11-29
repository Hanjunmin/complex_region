##python
import pandas as pd
import re
import json
import argparse
parser = argparse.ArgumentParser(description='INTEGRATE THE TEINS')
parser.add_argument('--chr', help='chromosome')
parser.add_argument('--json', help='json path')
args = parser.parse_args()


chr=args.chr


with open(args.json, 'rt', encoding='utf-8') as gz_file:
    fasta_dict = json.load(gz_file)


with open("seg.temid" , 'r') as f:
    lines = f.readlines()

aldic={}
for line in lines:
	lis=re.split(r'>|\s', line)
	aldic[lis[1]]=lis[2]


sdpa=pd.read_csv("output_n.fasta_rm.bed",sep="\t",header=None)
sdpa['name']=sdpa[0].apply(lambda x: aldic.get(x))
split_cols = sdpa.iloc[:, 10].str.split(r'[@,]', expand=True)
split_cols.columns = ["sam","x","conitig","contigs","contige","nodes","nodee","pos","pose"]
split_cols=pd.concat([sdpa, split_cols], axis=1)
split_cols['fas'] = list(map(fasta_dict.get, split_cols['nodes']))
split_cols['fae'] = list(map(fasta_dict.get, split_cols['nodee']))
split_cols['len1']=split_cols['fas'].str.len()
split_cols['len2']=split_cols['fae'].str.len()
split_cols['pos1s'] = split_cols['pos'].astype(float)
split_cols['pos1e'] = split_cols['pos'].astype(float)+split_cols['len1']
split_cols['pos2s'] = split_cols['pose'].astype(float)
split_cols['pos2e'] = split_cols['pose'].astype(float)+split_cols['len2']
# split_cols[split_cols['nodes']=="-1"]['pos1e']=split_cols[split_cols['nodes']=="-1"]['pos2s'] ##前面没有比对
# split_cols[split_cols['nodee']=="-2"]['pos1e']=split_cols[split_cols['nodes']=="-1"]['pos2s'] ##后面没有比对
###把上面两种过滤掉吧
split_cols=split_cols[split_cols['nodes']!="-1"]
split_cols=split_cols[split_cols['nodee']!="-2"]
x=split_cols[split_cols['pos1e']<=split_cols['pos2s']]
x['s']=x['pos1e']
x['e']=x['pos2s']
y=split_cols[split_cols['pos1e']>split_cols['pos2s']]
y['s']=y['pos2e']
y['e']=y['pos1s']
df = pd.concat([x, y], axis=0)
x=df[df['nodes']==df['nodee']]
y=df[df['nodes']!=df['nodee']]
x['s']=x['pos1s']
x['e']=x['pos1e']
df = pd.concat([x, y], axis=0)
df['CHR']=str(chr)
df["s"]=df["s"].astype("int")
df["e"]=df["e"].astype("int")
end=df.loc[:,["CHR","s","e","sam",1,2,3,4,5,6,7,"conitig","contigs","contige"]]
end.to_csv("reg.fa_rm.out", sep="\t", index=False, header=False)
extra=end[end['s']>end['e']]
x=pd.DataFrame()
if len(extra)!=0:
	x.to_csv("reg.err", sep="\t", index=False, header=False)




with open("seg.temid" , 'r') as f:
    lines = f.readlines()

alis=[]
for line in lines:
	lis=re.split(r'>|\s', line)
	sam=lis[1].split("NN")[0]
	douli=lis[2].split(",")
	regs=douli[3]
	rege=douli[4]
	nodes=douli[1]
	nodee=douli[2]
	atlis=douli[0].split("@")
	alis.append([sam,regs,rege,nodes,nodee,atlis[2],atlis[3],atlis[4]])


split_cols=pd.DataFrame(alis)
split_cols.columns=["sam","pos","pose","nodes","nodee","conitig","contigs","contige"]
split_cols['fas'] = list(map(fasta_dict.get, split_cols['nodes']))
split_cols['fae'] = list(map(fasta_dict.get, split_cols['nodee']))
split_cols['len1']=split_cols['fas'].str.len()
split_cols['len2']=split_cols['fae'].str.len()
split_cols['pos1s'] = split_cols['pos'].astype(float)
split_cols['pos1e'] = split_cols['pos'].astype(float)+split_cols['len1']
split_cols['pos2s'] = split_cols['pose'].astype(float)
split_cols['pos2e'] = split_cols['pose'].astype(float)+split_cols['len2']
split_cols=split_cols[split_cols['nodes']!="-1"]
split_cols=split_cols[split_cols['nodee']!="-2"]
x=split_cols[split_cols['pos1e']<=split_cols['pos2s']]
x['s']=x['pos1e']
x['e']=x['pos2s']
y=split_cols[split_cols['pos1e']>split_cols['pos2s']]
y['s']=y['pos2e']
y['e']=y['pos1s']
df = pd.concat([x, y], axis=0)
x=df[df['nodes']==df['nodee']]
y=df[df['nodes']!=df['nodee']]
x['s']=x['pos1s']
x['e']=x['pos1e']
df = pd.concat([x, y], axis=0)
df['CHR']=str(chr)
df["s"]=df["s"].astype("int")
df["e"]=df["e"].astype("int")
end=df.loc[:,["CHR","s","e","sam","conitig","contigs","contige"]]
end.to_csv("reg.fa_rm.alseq", sep="\t", index=False, header=False)
extra=end[end['s']>end['e']]
x=pd.DataFrame()
if len(extra)!=0:
	x.to_csv("reg.err1", sep="\t", index=False, header=False)