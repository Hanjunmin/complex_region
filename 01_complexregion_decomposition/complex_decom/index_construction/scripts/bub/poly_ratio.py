import subprocess
#from Bio import SeqIO
import re
from multiprocessing import Pool
import time
import argparse
# import numpy as np
from sklearn.cluster import DBSCAN

import pandas as pd 
import re
import json
from Bio.Seq import Seq
import lmdb
import numpy as np
import math
import numpy as np
from sklearn.cluster import DBSCAN
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='simnode inte file')
parser.add_argument('--reg', help='region')
args = parser.parse_args()

region=args.reg
def split_fasta(fasta_dict,splitbase):
    ####1.split fasta
    # splitbase=10 ###10bp as the split unit 
    fasta_dictsplit={}
    for index,value in fasta_dict.items():
        if len(value)>10:
            splilen=math.floor(len(value) / splitbase)
            fasta_dictsplit[index]=[str(index) for i in range(0, splilen)]
            # fasta_dictsplit[index]=[str(index) +"#"+ str(i) for i in range(0, splilen)]
        else:
            # fasta_dictsplit[index]=[str(index) +"#0"]
            fasta_dictsplit[index]=[str(index)]
    return fasta_dictsplit


def splitpath(in1):
    ###生成每条path对应的segement和对应的split后的多个顶点的位置 in1out in1df
    in1df=[]
    start=0
    in1out = []
    for item in in1:
        if item in keys_set:  # 快速查找
            in1out.extend(fasta_dictsplit[item])
            end=start+len(fasta_dictsplit[item])    
        else:
            in1out.append(item)
            end=start+1
        in1df.append([item,start,end])
        start=end
    in1df=pd.DataFrame(in1df)
    in1df['clus']=-1
    in1df[2]=in1df[2]-1
    return [in1df,in1out]


import sys
import os


file_path = "NODE.P."+region+"data"
if not os.path.exists(file_path):
    polyratio=1
    with open('polyratio'+region+'.txt', 'w', encoding='utf-8') as f:
        f.write(str(polyratio))
    sys.exit()




S_dict={}
dupelement=[]
with open(file_path, 'r') as file:
    for line in file:
        linelis=line.split("\t")
        S_dict[linelis[1]]=linelis[2].replace("+","").replace("-","").split(",")
        dupelement.extend(linelis[2].replace("+","").replace("-","").split(","))



dupelement=set(dupelement)
allpath = pd.read_csv("temp."+region+".allpath",sep="\t",header=None)

sdpa=pd.read_csv("temp."+region+".S.start",sep="\t",header=None)
sdpa[0]=sdpa[0].astype(str)
fasta_dict=dict(zip(sdpa.iloc[:, 0], sdpa.iloc[:, 1]))
splitbase=10
fasta_dictsplit=split_fasta(fasta_dict,splitbase)



# with open('/home/jmhan/pan-SD/CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
#     fasta_dict = json.load(gz_file)
keys_set = set(fasta_dictsplit.keys())



import numpy as np
dupelemental=[]
for index,row in allpath.iterrows():
    in1=str.split(row[2].replace("-",",-").replace("+",",+"),",")
    in1df,in1out=splitpath(in1)
    if len(dupelement)!=0:
        positions = np.where(np.isin(in1out, list(dupelement)))[0].tolist()
        print(len(positions))
        if len(positions)==0:
            continue
        arr_reshaped = np.array(positions).reshape(-1, 1)
        dbscan = DBSCAN(eps=100, min_samples=1)  ###如果位置差异大于20就分为另一类
        label=dbscan.fit_predict(arr_reshaped)
        for lab in set(label):
            res = np.array(positions)[label == lab]
            in1df.loc[(in1df[1]<=max(res)) & (in1df[2]>=min(res)),'clus']=lab
            xt=in1df.loc[(in1df[1]<=max(res)) & (in1df[2]>=min(res))]
            res=[xt.index[0],xt.index[-1]]
            pa="".join(in1[min(res):(max(res)+1)]).replace("+","+,").replace("-","-,")+"+"
            dupelemental.extend(in1[min(res):(max(res)+1)])
    dupelemental=list(set(dupelemental))


sdpa['clus']=-1
sdpa.loc[sdpa[0].isin(dupelemental), "clus"] = 1
sdpa['len']=sdpa[1].str.len()

polyratio=sum(sdpa.loc[sdpa['clus']==-1,"len"])/sum(sdpa["len"])
with open('polyratio'+region+'.txt', 'w', encoding='utf-8') as f:
    f.write(str(polyratio))