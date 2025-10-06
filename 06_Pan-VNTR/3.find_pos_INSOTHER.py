import pandas as pd
from collections import defaultdict
import lmdb
import re
import msgpack
import warnings# 忽略所有警告
import argparse
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description='Pan-VNTR-coordi')
parser.add_argument('--S', help='corresponding sample name')
parser.add_argument('--type', help='INS or other')
args = parser.parse_args()


sam=args.S
if args.type=="INS":
    nodecolumn="V12"
else:
    nodecolumn="V9"

Wfile="./decompose/splinew/"+sam+".W"
env = lmdb.open('./project/pansd/pandb/CHM13-APGp1-HPRCp1-HGSVCp3_MC_db_msgpack', create=False,readonly=True, lock=False)  # 10 GB


with open(Wfile, 'r') as f:
    lines = f.readlines()

local_dict = defaultdict(list)  # 使用局部字典减少锁竞争
for line in lines:
    line_list = line.strip().split('\t')
    samname=line_list[1]+"@"+line_list[2]+"@"+line_list[3]+"@"+line_list[4]+"@"+line_list[5]
    local_dict[samname] = line_list[6]


data=pd.read_csv(sam)
keys = set(data['V7'].astype(str))
local_dict = {k: local_dict[k] for k in keys if k in local_dict}



allist=[]
for key, group in data.groupby('V7'):
    print(key)
    value=local_dict[key]
    symbols = re.findall(r'[<>]',value)
    numbers = re.findall(r'\d+', value)
    end = pd.DataFrame()
    end['ori']=symbols
    end['node']=numbers
    with env.begin(write=False) as txn:
        cursor = txn.cursor()
        values = cursor.getmulti([key.encode('utf-8') for key in numbers])
        end['fa'] = [msgpack.unpackb(value[1]) for value in values if value[1] is not None]
    end['base_count'] = end['fa'].apply(lambda x: len(x))
    end['pos']=end['base_count'].cumsum()
    end['pos_s']=end['pos']-end['base_count']
    filtered_end = end[end['node'].isin(group[nodecolumn].astype(str).tolist())]
    for index,row in group.iterrows():
        x=filtered_end[filtered_end['node']==str(row[nodecolumn])]
        row['pos']="@".join(x['pos_s'].astype(str) + "_" + x['pos'].astype(str))
        # row['pos']=int(list(x['pos_s'])[0])
        # row['e']=int(list(x['pos'])[0]) ####这里只取了第一个位置，如果一个node对应多个位置要考虑的话之后要考虑，也就是说VNTR设置的那一个节点去确定样本位置并不可行  不唯一
        allist.append(list(row))


df=pd.DataFrame(allist)
df['subchr'] = df[6].str.split('@', expand=True)[2]
df['s'] = df[6].str.split('@', expand=True)[3]
# df['ns'] = pd.to_numeric(df['s'], errors='coerce') + df[12]
# df['ne'] = pd.to_numeric(df['s'], errors='coerce') + df[13]
# df.drop([12, 13,'s'], axis=1, inplace=True)

df.to_csv(sam+".out", sep='\t', index=False)
