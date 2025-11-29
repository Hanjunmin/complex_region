
import argparse
import pandas as pd
parser = argparse.ArgumentParser(description='merge-SDR')
parser.add_argument('--sam', help='sample name')
args = parser.parse_args()

sam=args.sam
import pybedtools
def iter(sdpa):
	grouped = sdpa.groupby([3, 4, 5])
	sdpacopy=sdpa.copy()
	for name, group in grouped:
		alcov=group[(group[4]>=group[1]) & (group[5]<=group[2])] ##SDR被完全覆盖
		if len(alcov)!=len(group):
		# break
			#subcov=group[~((group[4]>=group[1]) & (group[5]<=group[2]))] ##SDR没有被完全覆盖
			s=min(group[1])
			e=max(group[2])
			group[1]=s
			group[2]=e
			mask = (sdpacopy[1]>=s) & (sdpacopy[2]<=e)
			sdpacopy.loc[mask, 1] = s
			sdpacopy.loc[mask, 2] = e
	return sdpacopy.iloc[:,0:3].drop_duplicates()


sdpa=pd.read_csv(sam+".win.outpos",sep="\t",header=None)
end=iter(sdpa)
	
end.to_csv(sam+"path.txt", sep="\t", index=False, header=False)
lastlen=len(end)

x=0
while x==0:
	a = pybedtools.BedTool(sam+"path.txt")
	b = pybedtools.BedTool(sam+".outpos")
	result = a.intersect(b, wa=True, wb=True)
	df = result.to_dataframe()
	df.columns=[0, 1, 2, 3, 4, 5, 6]
	end=iter(df)
	id=len(end)
	if id==lastlen:
		# end.to_csv(sam+".win.txt", sep="\t", index=False, header=False)
		break
	else:
		lastlen=id


newen = end.sort_values(end.columns[1]).reset_index(drop=True)
newen['cluster']=0
clus=0
for i in range(len(end)):
	newen.iloc[i,3]=clus
	clus=clus+1
	if end.iloc[i,1]<=end.iloc[i-1,2] and end.iloc[i,1]>=end.iloc[i-1,1]:
		A=max(end.iloc[i,1],end.iloc[i-1,1])
		B=min(end.iloc[i,2],end.iloc[i-1,2])
		if B-A>500:
			newen.loc[i,"cluster"]=newen.loc[i-1,"cluster"]
	if end.iloc[i,1]>=end.iloc[i-1,1] and end.iloc[i,2]<=end.iloc[i-1,2]:
		newen.loc[i,"cluster"]=newen.loc[i-1,"cluster"]
		# print(newen)
		# print(B-A)
 
newen.columns=['chr','start','end','cluster']
merged_df = newen.groupby('cluster', as_index=False).agg({
	'chr':'first',
    'start': 'min',
    'end': 'max'
})


merged_df.iloc[:,[1,2,3]].to_csv(sam+".win.txt", sep="\t", index=False, header=False)
