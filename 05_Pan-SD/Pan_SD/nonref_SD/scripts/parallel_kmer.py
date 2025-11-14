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
import h5py


from concurrent.futures import ProcessPoolExecutor, as_completed
from datasketch import MinHash

parser = argparse.ArgumentParser(description='Pan-SD-KMER')
parser.add_argument('--fa', help='pangenome fasta')
parser.add_argument('--i', help='input newwin')
parser.add_argument('--n', help='number of cores')
args = parser.parse_args()
fafile=args.fa
winsimin=args.i
nc=args.n


def process_line(line):
    fa = pysam.FastaFile(fafile)
    ### /home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/t2tsdminimap/HG00438_Pat.fa
    line_list = line.strip().split('\t')
    chr1 = line_list[0]
    start1 = int(line_list[1])
    end1 = int(line_list[2])
    name = chr1 + ":" + str(start1) + "-" + str(end1)
    Seq1 = fa.fetch(region=name)
    kmer1 = kmer_cal(Seq(Seq1), 11)
    m = MinHash(num_perm=num_perm)
    m.update_batch([kmer.encode('utf8') for kmer in kmer1])
    return m, name





def kmer_cal(seq,size):
    kmers = [str(seq[i:i+size]) for i in range(len(seq) - size + 1)]
    kmers = list(filter(lambda kmer: 'N' not in kmer, kmers))
    rkmers = [str(seq[i:i+size].reverse_complement()) for i in range(len(seq) - size + 1)]
    rkmers = list(filter(lambda kmer: 'N' not in kmer, rkmers))
    result = [min(kmer1, kmer2, key=lambda x: kmerhasher(x, size)) for kmer1, kmer2 in zip(kmers, rkmers)]
    result= {x.lower() for x in result}
    return(set(result))

from Bio.Seq import Seq



# faname=args.fa
# # fa = pysam.FastaFile("../t2tsdminimap/HG00438_Pat.fa")
# fa = pysam.FastaFile("/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/t2tsdminimap/HG00438_Pat.fa")

# win=args.i "/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/other/newwin"
# windf = pd.read_csv(winsimin, delimiter='\t',header=None)
# windf.columns = ['chr', 'start', 'end']
with open(winsimin, 'r') as f:
    lines = f.readlines()


num_perm = 128
nowsamplemh = []
nowsamplemhname = []
from concurrent.futures import ProcessPoolExecutor, as_completed
with ProcessPoolExecutor(max_workers=int(nc)) as executor: ##大概24万条
    futures = [executor.submit(process_line, line) for line in lines]

for i, future in enumerate(as_completed(futures)):
    result = future.result()
    if result:
        m, name = result
        nowsamplemh.append(m)
        nowsamplemhname.append(name)

print("All lines processed.")



with open('records.pkl', 'wb') as file:
    pickle.dump({'var1':nowsamplemh, 'var2': nowsamplemhname}, file)


# with open('/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/other/128filt0.3/records.pkl', 'rb') as file:
#     loaded_vars = pickle.load(file)

# nowsamplemh=loaded_vars['var1']
# nowsamplemhname=loaded_vars['var2']
lsh = MinHashLSH(threshold=0.1, num_perm=128)  # 设置阈值和哈希函数数量
for i, m in enumerate(nowsamplemh):
    lsh.insert(str(i), m)




# with open('samplemhnum', 'w') as file:
#     file.write(str(len(nowsamplemh)))

import pickle
from concurrent.futures import ProcessPoolExecutor

def process_chunk(chunk, start_idx):
    otherid1 = []
    pairset1 = {}
    jacset1 = {}
    for i, m in enumerate(chunk):
        # print(i)
        real_idx = start_idx + i  # 真实索引
        if real_idx % 10000 == 0 and real_idx != 0:
            with open(f"{real_idx}end.pkl", 'wb') as file:
                pickle.dump({'var1': otherid1, 'var2': pairset1, 'var3': jacset1}, file)
            otherid1.clear()
            pairset1.clear()
            jacset1.clear()
        results = lsh.query(m)
        if len(results) >= 100000:
            otherid1.append(nowsamplemhname[real_idx])
            continue
        numbers = []
        jac = []
        for j in results:
            if int(j) > real_idx:
                numbers.append(int(j))
                jac.append(m.jaccard(nowsamplemh[int(j)]))
        modified_list = [nowsamplemhname[num] for num in numbers]
        pairset1[nowsamplemhname[real_idx]] = modified_list
        jacset1[nowsamplemhname[real_idx]] = jac
    with open(f"{real_idx}end.pkl", 'wb') as file:
        pickle.dump({'var1': otherid1, 'var2': pairset1, 'var3': jacset1}, file)
    return otherid1, pairset1, jacset1

# 设置每个子列表的大小
chunk_size = len(nowsamplemh) // 50  # 假设使用10个进程
chunks = [nowsamplemh[i:i + chunk_size] for i in range(0, len(nowsamplemh), chunk_size)]

otherid1_final = []
pairset1_final = {}
jacset1_final = {}

# 使用多进程处理每个子列表
with ProcessPoolExecutor(max_workers=int(nc)) as executor:
    futures = [executor.submit(process_chunk, chunk, idx * chunk_size) for idx, chunk in enumerate(chunks)]
    for future in futures:
        otherid1_part, pairset1_part, jacset1_part = future.result()
        otherid1_final.extend(otherid1_part)
        pairset1_final.update(pairset1_part)
        jacset1_final.update(jacset1_part)

print("All chunks processed.")


# with open(str(i)+'end.pkl', 'wb') as file:
#     pickle.dump({'var1':otherid1, 'var2': pairset1,'var3': jacset1}, file)

