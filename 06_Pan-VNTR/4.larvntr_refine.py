

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


def check_lef_rig(checkin):
	if checkin['INS_samS']==checkin['Seq_samS']:
		return "LEFT" ###ref在右边，ins左边不完整
	if checkin['INS_samE']==checkin['Seq_samE']:
		return "RIGHT" ###ref在左边，ins右边不完整


def leffun(index1,lef):
	cumulative_sum=0 
	filtered_reflef= df2hap[df2hap.index < index1]
	if len(filtered_reflef)!=0:
		stillnd=lef
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
	result_df_ref=result_df.dropna(subset=['Value'])
	newposs=min(result_df_ref['T2Tpos'])
	newpose=max(result_df_ref['end'])
	return newposs,newpose,sum(result_df_ref['base_count'])/sum(result_df['base_count'])


def rigfun(index2,rig):
	cumulative_sum=0
	filtered_refrig= df2hap[df2hap.index >index2]
	if len(filtered_refrig)!=0:
		stillnd=rig
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
	result_df_ref=result_df.dropna(subset=['Value'])
	newposs=min(result_df_ref['T2Tpos'])
	newpose=max(result_df_ref['end'])
	return newposs,newpose,sum(result_df_ref['base_count'])/sum(result_df['base_count'])

##################################################################################################################################################################
from collections import defaultdict
import pandas as pd
from Bio.Seq import Seq
import re
import argparse
from datetime import datetime
import sys
import lmdb
import msgpack
sys.path.append('./project/ENDVNTR')
from vntr_fun import *
parser = argparse.ArgumentParser(description='REFINE-LARGE')
parser.add_argument('--SAM', help='SAMPLE')
args = parser.parse_args()

allWfold="./decompose/splinew/"
posfold="./project/panVNTR/data/"
env = lmdb.open('./project/pansd/pandb/CHM13-APGp1-HPRCp1-HGSVCp3_MC_db_msgpack', create=False,readonly=True, lock=False)  # 10 GB

sam=args.SAM
SAMLARVNTR = pd.read_csv(sam+"larvntr.refine",sep="\t",header=0)
chrpathsam=set(list(SAMLARVNTR['sampath']))
Wfile=allWfold+sam+".W"
aldicmulti=process_file2(Wfile)
aldicmultipath=process_file3(Wfile)


SAMLARVNTR=SAMLARVNTR[SAMLARVNTR['RP_samE']>SAMLARVNTR['INS_samS']]
SAMLARVNTR['STRAT']=0
SAMLARVNTR['END']=0
SAMLARVNTR['SIM']=0
SAMLARVNTR['TYPE']="a"


for larindex,larline in SAMLARVNTR.iterrows():
	larline=SAMLARVNTR.loc[larindex]
	CHRn=larline['REFCHR']
	T2Tchmpos=posfold+CHRn+".posCHM.txt"
	T2Tpos = pd.read_csv(T2Tchmpos,sep="\t",header=None)
	T2Tpos.columns=["CHR","start_position","end_position","node","orient"]
	T2Tpos['node'] = T2Tpos['node'].astype(str)
	node_to_start_dict = T2Tpos.set_index('node')['start_position'].to_dict() ##存放了T2T的node对应的位置信息
	pos_to_start_dict = T2Tpos.set_index('start_position')['node'].to_dict() ##存放了T2T的node对应的位置信息
	allnode=T2Tpos['node']
	allnodeset=set(allnode)
	interdict={} ###存放每个样本与T2Tnode有交集的node以及对应的pos
	for target_value in set(chrpathsam):
		interset=allnodeset & set(aldicmulti[target_value]) ### target_value="HG01243_1_JAHEOY010000244.1_3031_869230"
		new_dict = {key: node_to_start_dict.get(key) for key in interset}
		interdict[target_value]=new_dict
	target_value=larline['sampath']
	othervntr=[]
	df = pd.DataFrame(list(interdict[target_value].items()), columns=['Key', 'Value'])
	insert_list=aldicmulti[target_value]
	pa=aldicmultipath[target_value]
	symbols = re.findall(r'[<>]',pa)
	insert_df = pd.DataFrame({'Key': insert_list})
	df_combined = pd.merge(insert_df, df, on='Key', how='left')
	with env.begin(write=False) as txn:
		cursor = txn.cursor()
		values = cursor.getmulti([key.encode('utf-8') for key in df_combined['Key']])
		df_combined['fa'] = [msgpack.unpackb(value[1]) for value in values if value[1] is not None]
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
	df_mark_zero = df_combined[df_combined['mark'] == 0].copy()
	df_mark_zero['cluster'] = (df_mark_zero.index.to_series().diff() != 1).cumsum()
	df2hap=df2hapdf(df_combined,df_mark_zero)
	nan_rows = df2hap[df2hap['Value'].isna()].copy()
	nan_rows.loc[:, 'fa_length'] = nan_rows['newfa_concat'].apply(lambda x: len(x) if pd.notna(x) else 0)
	non_nan_rows = df2hap[~df2hap['Value'].isna()].copy()###*******
	nanclus=larline['cluster'] 
	cho=nan_rows[nan_rows['cluster']==nanclus]
	keyNod=list(cho['Key'])[0];keyNod1=list(cho['Key'])[-1]
	index1=cho.index[0]
	index2=cho.index[-1]
	nowins=df2hap.loc[index1]
	insfa=nowins['newfa_concat']
	check_lef_rig(larline) ###为RIGHT 
	print(larline['RP_samS'],larline['RP_samE'],larline['INS_samS'],larline['INS_samE'])
	rms=larline['RP_samS']
	rme=larline['RP_samE']
	if check_lef_rig(larline)=="RIGHT": ##说明REF在右边，ins左边不完整
		inss=larline['INS_samS']
		inse=inss+len(insfa)
		lef=inss-rms
		rig=rme-inse
		if lef>0 and rig>0:
			newposs1,newpose1,sim1=leffun(index1,lef)
			newposs2,newpose2,sim2=rigfun(index2,rig)
			newposs=min(newposs1,newposs2)
			newpose=max(newpose1,newpose2)
			sim=min(sim1,sim2)
		if rig>0 and lef<=0:
			newposs,newpose,sim=rigfun(index2,rig)
		if rig<=0 and lef>0:
			newposs,newpose,sim=leffun(index1,lef)
		type=check_lef_rig(larline)
		# print(larindex,check_lef_rig(larline),lef,rig,newposs,newpose,sim)
	if check_lef_rig(larline)=="LEFT": ##说明REF在左边，ins右边不完整
		inse=larline['INS_samE']
		inss=inse-len(insfa)
		lef=inss-rms
		rig=rme-inse
		if lef>0 and rig>0:
			newposs1,newpose1,sim1=leffun(index1,lef)
			newposs2,newpose2,sim2=rigfun(index2,rig)
			newposs=min(newposs1,newposs2)
			newpose=max(newpose1,newpose2)
			sim=min(sim1,sim2)
		if rig>0 and lef<0:
			newposs,newpose,sim=rigfun(index2,rig)
		if rig<0 and lef>0:
			newposs,newpose,sim=leffun(index1,lef)
		# print(larindex,check_lef_rig(larline),lef,rig,newposs,newpose,sim)
		type=check_lef_rig(larline)
	if (lef <= 0) & (rig <= 0) & (check_lef_rig(larline) == "RIGHT"):
		newposs=df2hap.loc[index1-1]['Value']
		newpose=0
		sim=0
	if (lef<=0) & (rig<=0) & (check_lef_rig(larline)=="LEFT"):
		newposs=df2hap.loc[index2+1]['Value']
		newpose=0
		sim=0
	SAMLARVNTR.loc[larindex,'SIM']=sim
	SAMLARVNTR.loc[larindex,'STRAT']=newposs
	SAMLARVNTR.loc[larindex,'END']=newpose
	SAMLARVNTR.loc[larindex,'TYPE']=newpose


SAMLARVNTR.to_csv(sam+'lar.refine.out', sep="\t")  # 不保存索引列
