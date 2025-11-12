####一个快速的脚本 
#### 图上某个样本和参考基因组相比的segment的差异


import pandas as pd
import json
import re
T2Tpos = pd.read_csv("/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.posCHM.txt",sep="\t",header=None)
T2Tpos.columns=["CHR","start_position","end_position","node","orient"]
T2Tpos['node'] = T2Tpos['node'].astype(str)
node_to_start_dict = T2Tpos.set_index('node')['start_position'].to_dict() ##存放了T2T的node对应的位置信息

file="/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/Wpath/C021-CHA-S01-Pat.W"
with open('/home/jmhan/project/APG/github_final/01_1p36/DATA/graph/reg.json', 'rt', encoding='utf-8') as gz_file:
    fasta_dict = json.load(gz_file)


with open(file, 'r') as f:
    lines = f.readlines()

dfal = pd.DataFrame()  # 初始化为空的DataFrame
for line in lines:
    line_list = line.strip().split('\t')
    symbols = re.findall(r'[<>]',line_list[6])
    numbers = re.findall(r'\d+', line_list[6])
    df = pd.DataFrame(numbers)
    df['fa'] = list(map(fasta_dict.get, df[0]))
    df['ori'] = symbols
    df['base_count'] = df['fa'].apply(lambda x: len(x))
    df['pos']=df['base_count'].cumsum()
    df['T2Tpos'] = df[0].apply(lambda x: node_to_start_dict.get(x))
    df['T2Tend'] = df['T2Tpos']+df['base_count']
    dfal = pd.concat([dfal, df], ignore_index=True)

dfal.to_csv('test.csv', index=False,sep="\t")