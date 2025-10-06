from collections import defaultdict
import pandas as pd
from Bio.Seq import Seq
import re
import argparse
from datetime import datetime
import sys
import lmdb
import msgpack
import os
sys.path.append('./project/ENDVNTR')
from vntr_fun import *
parser = argparse.ArgumentParser(description='REFINE-LARGE')
parser.add_argument('--SAM', help='SAMPLE')
args = parser.parse_args()
sam=args.SAM

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
			if indexot==index:
				continue
			else:
				ratio=(min(noe,value[2])-max(nos,value[1]))/(value[2]-value[1]) ###intersect/compare
				if (ratio>=0.95) & (value[0]==row['SAM_CHR']):
					iterdict.pop(indexot) 
	df=df.loc[iterdict.keys()]
	df=df.drop(['len'],axis=1)
	return df



fold= os.getcwd()
smvntr_file = os.path.join(fold, sam+"smvntr.out")
larvntr_file = os.path.join(fold,sam+"lar.refine.out")

INScolumn=['REFCHR','REFPOS','TYPE','SAM_CHR','CN','SAM','MOTIF_LEN','RP_samS','RP_samE','MOTIF','motifsim']
if os.path.exists(smvntr_file):
	smVNTR = pd.read_csv(smvntr_file,sep="\t",header=0)
	INSsm=smVNTR[['REFCHR','REFPOS','TYPE','SAM_CHR','CN','SAM','MOTIF_LEN','RP_samS','RP_samE','MOTIF','motifsim']]
	INSsm=inter_del(INSsm)
else:
	INSsm=pd.DataFrame()


if os.path.exists(larvntr_file):
	larVNTR = pd.read_csv(larvntr_file,sep="\t",header=0)
	INSlar=larVNTR[larVNTR['END']==0]
	INSlar=INSlar[['REFCHR','REFPOS','TYPE','SAM_CHR','CN','SAM','MOTIF_LEN','RP_samS','RP_samE','MOTIF','motifsim']]
	INSlar.columns=INScolumn
	if len(INSlar)!=0:
		INSlar=inter_del(INSlar)
	else:
		INSlar=pd.DataFrame()
	OTHERlar=larVNTR[larVNTR['END']!=0]
	OTHERlar=OTHERlar[['REFCHR','STRAT','END','SAM_CHR','CN','SAM','MOTIF_LEN','RP_samS', 'RP_samE', 'MOTIF','motifsim']]
	OTHERlar=inter_del(OTHERlar)
else:
	OTHERlar=pd.DataFrame()
	INSlar=pd.DataFrame()


combined_df = pd.concat([INSlar, INSsm], ignore_index=True)
combined_df.to_csv(sam+".INS.END", index=False,sep="\t")  # Saves without the index column
OTHERlar.to_csv(sam+".OTHER.END", index=False,sep="\t")  # Saves without the index column

