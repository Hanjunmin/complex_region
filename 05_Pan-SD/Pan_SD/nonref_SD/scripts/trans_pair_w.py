import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='trans_Pair')
parser.add_argument('--pos', help='position file')
parser.add_argument('--W', help='W bed')
parser.add_argument('--o', help='output file')
parser.add_argument('--trans', help='trans file')
args = parser.parse_args()
# sampos = pd.read_csv('C001-CHA-E01-Mat.pos',sep="\t",header=None)
# initW = pd.read_csv('C001-CHA-E01-Mat.init.W',sep="\t",header=None)


sampos = pd.read_csv(args.pos,sep="\t",header=None)
initW = pd.read_csv(args.W,sep="\t",header=None)

sampos.columns=["sd",'s','e','from','xx','sim']
initW.columns=["sam",'no','chr','ps','pe','from']
merged = pd.merge(sampos, initW, on='from')  # 内连接
merged['ns']=merged['ps']+merged['s']
merged['ne']=merged['ps']+merged['e']
mergedend=merged[['sd','sam','chr','ns','ne','sim']]

mergedendpair=mergedend.copy()
mergedendpair['pair']=mergedendpair['chr']+ ':'+mergedendpair['ns'].astype(str)+'-' +mergedendpair['ne'].astype(str)
transpair=mergedendpair[['sd','pair']]
#allpair=pd.read_csv('/share/home/zhanglab/user/maoyafei/project/pansd/iter/other/pairall/tranall.uniq',sep="\t",header=None)

allpair=pd.read_csv(args.trans,sep="\t",header=None)


pairs_dict = dict(zip(transpair['sd'], transpair['pair']))
endpairall=allpair.copy()
endpairall[3]=allpair[0].map(pairs_dict).fillna("NA")
endpairall[4]=allpair[1].map(pairs_dict).fillna("NA")
allpair[0] = allpair[0].map(pairs_dict).fillna("NA")
allpair[1] = allpair[1].map(pairs_dict).fillna("NA")
allpair = allpair[(allpair[0] != "NA") & (allpair[1] != "NA")]
allpair['sorted_pair'] = allpair.apply(lambda row: tuple(sorted([row[0], row[1]])), axis=1)
allpair_unique = allpair.drop_duplicates(subset='sorted_pair').drop(columns='sorted_pair')
# allpair.to_csv("C001-CHA-E01-Mat.uniqpair", index=False,sep="\t",header=None)
allpair.to_csv(args.o, index=False,sep="\t",header=None)
endpairall.to_csv(args.o+".all", index=False,sep="\t",header=None)
