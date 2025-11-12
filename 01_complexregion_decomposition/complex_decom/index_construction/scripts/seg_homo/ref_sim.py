import pandas as pd
import re
from Bio.Seq import Seq
import pickle
from kmerhash import kmerhasher, hashes2seq
import json
import multiprocessing as mp
import pysam



import argparse
parser = argparse.ArgumentParser(description='SD-Sim')
parser.add_argument('--d', help='inter-seg kmer database')
parser.add_argument('--r', help='reference path')
parser.add_argument('--w', help='chromsome window file')
parser.add_argument('--o', help='output file')
parser.add_argument('--p', help='processes of running')
args = parser.parse_args()


def kmer_cal(seq,size):
	kmers = [str(seq[i:i+size]) for i in range(len(seq) - size + 1)]
	kmers = list(filter(lambda kmer: 'N' not in kmer, kmers))
	rkmers = [str(seq[i:i+size].reverse_complement()) for i in range(len(seq) - size + 1)]
	rkmers = list(filter(lambda kmer: 'N' not in kmer, rkmers))
	result = [min(kmer1, kmer2, key=lambda x: kmerhasher(x, size)) for kmer1, kmer2 in zip(kmers, rkmers)]
	result= {x.lower() for x in result}
	return(set(result))


### 1.segment-fasta file
# with open ('./CHM13-APGp1-HPRCp1-HGSVCp3_MC.S.json', 'rt', encoding='utf-8') as gz_file:
#     fasta_dict = json.load(gz_file)
import lmdb
import msgpack


# with open('./integrate/sedef/SD_lowerall829.pkl', 'rb') as f:
# 	loaded_set = pickle.load(f)


with open(args.d+"inter_seg_kmer.pkl", 'rb') as f:
	loaded_set = pickle.load(f)

print("OK")
#fasta_dict['00000']="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
### 2.file title
def process_line(line):
	print(line)
	allsim={}
	cor=[]
	line_list = line.strip().split('\t')
	CHR=line_list[0]
	wstart=int(line_list[1])
	wend=int(line_list[2])
	# FA="./breakpoints/chm13v2.0.fa"
	FA=args.r
	fa1 = pysam.FastaFile(FA)
	samplesim={}
	namewin=CHR + ":" + str(wstart) + "-" + str(wend)
	samseq=str(fa1.fetch(region=CHR + ":" + str(wstart) + "-" + str(wend)))
	samkmer=kmer_cal(Seq(samseq),31)
	sim = len(samkmer & loaded_set) / len(samkmer)
	allsim[namewin]=sim
	return allsim

#with open("chr1.window", 'r') as f:
with open(args.w, 'r') as f:
	lines = f.readlines()

with mp.Pool(processes=int(args.p)) as pool:
    res=pool.map(process_line, lines)



merged_dict = {}
for d in res:
    merged_dict.update(d)

simdf = pd.DataFrame({
    'win': merged_dict.keys(),
    'T2TCHM13init': merged_dict.values()
})

simdf.to_csv(args.o,sep='\t')
#simdf_transposed.to_csv("chr1.sim",sep='\t')



