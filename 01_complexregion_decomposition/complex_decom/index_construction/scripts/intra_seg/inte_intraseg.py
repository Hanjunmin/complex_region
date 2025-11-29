from Bio.Seq import Seq
import edlib
import subprocess


def inter_del(df): ###delete the intersect lines(>0.9 similarity)
	iterdict={}
	df['len']=df['RP_samE']-df['RP_samS']
	df_sorted = df.sort_values(by='len', ascending=False)
	for index,row in df.iterrows():
		iterdict[index]=[row['SAM_CHR'],row['RP_samS'],row['RP_samE']]
	for index,row in df_sorted.iterrows():
		if index not in iterdict.keys():
			continue
		nindex=index
		nos=row['RP_samS']
		noe=row['RP_samE']
		for indexot,value in list(iterdict.items()):
			if indexot==index or value[0]==row['SAM_CHR']:
				continue
			else:
				ratio=(min(noe,value[2])-max(nos,value[1]))/(value[2]-value[1]) ###intersect/compare
				if (ratio>=0.95) & (value[0]==row['SAM_CHR']):
					iterdict.pop(indexot) 
	df=df.loc[iterdict.keys()]
	df=df.drop(['len'],axis=1)
	return df


def motif_sim(string_list):
    distances = {}
    n = len(string_list)
    for i in range(n):
        for j in range(i+1, n):
            distance1 = edlib.align(string_list[i], string_list[j],mode='global',task="path")['editDistance']
            distance2 = edlib.align(string_list[i], Seq(string_list[j]).reverse_complement(),mode='global',task="path")['editDistance']
            distance=min(distance1,distance2)
            distances[i, j] = distance/max(len(string_list[i]), len(string_list[j]))
    return distances

def net_compari(pairdict):
        G = nx.Graph()
        for row in pairdict:
                G.add_edge(row[0], row[1])
        connected_components = list(nx.connected_components(G))
        region_to_group = {}
        for i, component in enumerate(connected_components, start=1):
                for region in component:
                        region_to_group[region] = i
        return region_to_group



### py script
import pandas as pd
import networkx as nx
final_vntrdf = pd.DataFrame()
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import argparse
parser = argparse.ArgumentParser(description='Pan-VNTR-coordi')
parser.add_argument('--pos', help='dot analysis of position')
args = parser.parse_args()


df = pd.read_csv("END.out", sep="\t",header=0)
dfcolumns = df.columns.tolist() 
dfcolumns[2] ="CN"
dfcolumns[3] ="SAM"
dfcolumns[4] ="SAM_CHR"
dfcolumns[5] ="RP_samS"
dfcolumns[6] ="RP_samE"
dfcolumns[7] ="MOTIF"
dfcolumns[8] ="MOTIF_sim"
df.columns=dfcolumns
df["submotif"] = df["MOTIF"].str.split(",").str[0]
df['motif_len'] = df['submotif'].str.len()
alist=[]
alldict={}
sudict={}
sudictmotif={}
sudictcn={}
sudictsim={}
modic={}
for clus in set(df['name2']):
	clusdata=df[df['name2']==clus]
	vntrpro= defaultdict(list)
	for idx, row in clusdata.iterrows():
		vntrpro[row['SAM']].append(row['CN'])
	alldict[clus]=vntrpro
	sudict[clus]=max(clusdata['CN']*clusdata['motif_len'])
	sudictmotif[clus]=max(clusdata['CN'])
	sudictcn[clus]=max(clusdata['motif_len'])
	sudictsim[clus]=clusdata['MOTIF_sim'].mean()
	modic[clus]=",".join(set(",".join(clusdata['MOTIF']).split(",")))



for region, samples in alldict.items():
	for sample, value in samples.items():
		# if isinstance(value, list):  # 如果值是一个列表
		alldict[region][sample] = value[0]  # 求和 

othervntr=pd.DataFrame(alldict).T
othervntr['vntr']=othervntr.index 
othervntr['motiflen'] = othervntr['vntr'].map(sudictcn)
othervntr['filt'] = othervntr['vntr'].map(sudict)
othervntr['cn'] = othervntr['vntr'].map(sudictmotif)
othervntr['motifsim'] = othervntr['vntr'].map(sudictsim)
othervntr['motif'] = othervntr['vntr'].map(modic)


# final_vntrdf = pd.read_csv("FINAL_OTHERALL.csv", sep="\t")
ll=othervntr[(othervntr['filt']>=150) & (othervntr['motiflen']>6) & (othervntr['filt']<=10000)]
tes=ll[ll['cn']>2]



samcor = pd.read_csv(args.pos, sep="\t",header=None)  ##对每个样本进行纠正
# "/home/jmhan/project/APG/github_final/01_1p36/index_construction/dot_refine/allsampos"
for index,row in tes.iterrows():
	poslis=row['vntr'].split("|")
	if poslis[2]=="INS":
		poslis[2]=int(poslis[1])+1
	cho=samcor[(samcor[0]==poslis[0]) & (samcor[1]<=int(float(poslis[1])))  & (samcor[2]>=int(float(poslis[2])))]
	tes.loc[index, 'chrom'] = poslis[0]          # 染色体
	tes.loc[index, 'start_pos'] = int(float(poslis[1]))  # 起始位置
	tes.loc[index, 'end_pos'] = int(float(poslis[2]))    # 结束位置
	havelis=list(set(cho[3]))
	for col in havelis:
		if col not in tes.columns:
			tes[col] = None  # 或 tes[col] = 0
	tes.loc[index, havelis] = tes.loc[index, havelis].fillna(0)

tes =tes.fillna('-1') ###这里相当于是vcf中的dot


output_df = tes[['chrom','start_pos','end_pos']].copy()
output_df['start_pos'] = output_df['start_pos'].astype(int)
output_df['end_pos'] = output_df['end_pos'].astype(int)
output_df.to_csv('tes.tem', index=False, sep='\t',header=None)
cmd = "bedtools coverage -a tes.tem -b x | cut -f 1-3,7 > y"
subprocess.run(cmd, shell=True, check=True)
cov = pd.read_csv("y", sep="\t",header=None)
tes['cov'] = cov[3].values  # 使用.values获取值数组
tes=tes[tes['cov']<0.2]




tes.to_csv('FINAL_OTHERvntr.csv', index=False,sep="\t")

