import hashlib
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
from kmerhash import kmerhasher, hashes2seq
import pysam
import pickle
import pandas as pd
import hashlib
import time
from collections import OrderedDict
from collections import defaultdict
import argparse
from datasketch import MinHash
import scipy.sparse
import argparse
from datasketch import MinHash
from scipy.cluster.hierarchy import linkage, to_tree
from datasketch import MinHashLSH

parser = argparse.ArgumentParser("integrate pairs")
parser.add_argument("-fa", "--fa", type=str, required=True, help="the maskfa file")
parser.add_argument("-i", "--i", type=str, required=True, help="the process region file")
args = parser.parse_args()


def kmer_cal(seq,size):
    kmers = [str(seq[i:i+size]) for i in range(len(seq) - size + 1)]
    kmers = list(filter(lambda kmer: 'N' not in kmer, kmers))
    rkmers = [str(seq[i:i+size].reverse_complement()) for i in range(len(seq) - size + 1)]
    rkmers = list(filter(lambda kmer: 'N' not in kmer, rkmers))
    result = [min(kmer1, kmer2, key=lambda x: kmerhasher(x, size)) for kmer1, kmer2 in zip(kmers, rkmers)]
    result= {x.lower() for x in result}
    return(set(result))



faname=args.fa
fa = pysam.FastaFile(faname)
win=args.i
windf = pd.read_csv(win, delimiter='\t',header=None)
windf.columns = ['chr', 'start', 'end']
num_perm = 128



with open(win, 'r') as f:
    lines = f.readlines()


nowsamplemh = []
nowsamplemhname = []
i=0
for line in lines:
    print(i)
    i=i+1
    line_list=line.strip().split('\t',)
    chr1=line_list[0]
    start1=int(line_list[1])
    end1=int(line_list[2])
    name=chr1+":"+str(start1)+"-"+str(end1)
    Seq1 = fa.fetch(region=name)
    if set(Seq1) == {'N'}:
        continue
    kmer1=kmer_cal(Seq(Seq1),11)
    if len(kmer1)==0:
        continue
    m = MinHash(num_perm=num_perm)
    for kmer in kmer1:
        m.update(kmer.encode('utf8'))
    nowsamplemh.append(m)
    nowsamplemhname.append(name)


with open('records.pkl', 'wb') as file:
    pickle.dump({'var1':nowsamplemh, 'var2': nowsamplemhname}, file)


with open('records.pkl', 'rb') as file:
    loaded_vars = pickle.load(file)

nowsamplemh=loaded_vars['var1']
nowsamplemhname=loaded_vars['var2']
lsh = MinHashLSH(threshold=0.05, num_perm=128)  # 设置阈值和哈希函数数量
for i, m in enumerate(nowsamplemh):
    lsh.insert(str(i), m)




with open('samplemhnum', 'w') as file:
    file.write(str(len(nowsamplemh)))

snum=len(nowsamplemh)//10000*10000


otherid1=[]
pairset1={}
jacset1={}
for i, m in enumerate(nowsamplemh): 
    if i % 10000 == 0 and i!=0:
        with open(str(i)+'end.pkl', 'wb') as file:
            pickle.dump({'var1':otherid1, 'var2': pairset1,'var3': jacset1}, file)
        otherid1.clear()
        pairset1.clear()
        jacset1.clear()
    print(i)
    results = lsh.query(m)
    if len(results)>=100000:
        otherid1.append(nowsamplemhname[i])
        continue
    numbers = []
    jac=[]
    for j in results:
        if int(j)>int(i):
            numbers.append(int(j)) ## i line similary j
            jac.append(nowsamplemh[i].jaccard(nowsamplemh[int(j)])) # jac
    modified_list = [nowsamplemhname[num] for num in numbers]
    pairset1[nowsamplemhname[i]]=modified_list
    jacset1[nowsamplemhname[i]]=jac

with open(str(i)+'end.pkl', 'wb') as file:
    pickle.dump({'var1':otherid1, 'var2': pairset1,'var3': jacset1}, file)


