import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='trans_Pair')
parser.add_argument('--i', help='input file')
parser.add_argument('--o', help='output file')
parser.add_argument('--trans', help='trans file')
args = parser.parse_args()
# SDpair = pd.read_csv('C002-CHA-E02-Mat.SD.pair',sep="\t",header=None)
# transpair=pd.read_csv('/home/jmhan/iterall/pair/C002-CHA-E02-Mat.pair',sep="\t",header=None)
SDpair = pd.read_csv(args.i,sep="\t",header=None)
transpair=pd.read_csv(args.trans,sep="\t",header=None)
pairs_dict = dict(zip(transpair[1], transpair[0]))
SDpair[0] = SDpair[0].map(pairs_dict)
SDpair[1] = SDpair[1].map(pairs_dict)
df_cleaned = SDpair.dropna() ###为na的为在图上没有路径的区域
df_cleaned.to_csv(args.o, index=False,sep="\t",header=None)