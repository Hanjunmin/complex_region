import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='trans_Pair')
parser.add_argument('--pos', help='position file')
parser.add_argument('--W', help='W bed')
parser.add_argument('--o', help='output file')
parser.add_argument('--pair', help='PAIR PATH')
args = parser.parse_args()


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
transpair=mergedendpair[['sd','pair','sim']]

file_path = args.pair
allpair=pd.read_csv(file_path,sep="\t",header=None)
allpair.columns=['A', 'B']


transpair["pair_sim"] = transpair['pair'] + "@" + transpair['sim'].astype(str)
pairs_dict = dict(zip(transpair['sd'], transpair['pair_sim']))
endpairall=allpair.copy()
endpairall['A.sam']=""
endpairall['B.sam']=""
endpairall['A.sam']=allpair['A'].map(pairs_dict).fillna("NOSAM")
endpairall['B.sam']=allpair['B'].map(pairs_dict).fillna("NOSAM")
endpairall[["A.pos", "A.sim"]] = endpairall['A.sam'].str.split("@", expand=True)
endpairall[["B.pos", "B.sim"]] = endpairall['B.sam'].str.split("@", expand=True)


endpairall=endpairall[['A', 'B', 'A.pos', 'A.sim', 'B.pos', 'B.sim']]
endpairall.to_csv(args.o+".pos_sim", index=False,sep="\t")
