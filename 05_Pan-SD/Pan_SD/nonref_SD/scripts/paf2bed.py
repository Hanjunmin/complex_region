## module load bedtools
# cat /home/jmhan/HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/asm_repeatmasker.out.bed | grep Satellite | cut -f 1-3 | bedtools sort -i - |  bedtools merge -i - | \
# bedtools coverage -header -a/home/jmhan/pangenome/SD/callSD/aligntest/end.txt -b - > {output.tmp}
# grep "^#" {input.sedef} > {output.bed}
# cat {output.tmp} >> {output.bed}
# sed -i '1{{s/$/\tcount_ovls\tsat_bases\ttotal_bases\tsat_coverage/}}' {output.bed}


# df = df[df.sat_coverage <= args.sat]
# args.sat=0.7

import subprocess
#from Bio import SeqIO
import re
from multiprocessing import Pool
import time
import argparse
parser = argparse.ArgumentParser(description='CIGAR-CALCULATE')
# parser.add_argument('--p', help='number of process')
# parser.add_argument('--r', help='Reference fasta')
# parser.add_argument('--q', help='query fasta')
parser.add_argument('--paf', help='alignment paf file')
parser.add_argument('--o', help='output file')
args = parser.parse_args()
# print(args.p)


# QUERY =args.q
# REF =args.r
# paf=args.paf
# #paf="/home/jmhan/CIGAR/T2T/cigartest.paf"
# subprocess.run(f"awk '{{print $1,$3,$4,\"\",\"\",$5}}' {paf}|tr ' ' '\t'  > query.bed", shell=True)
# subprocess.run(f"bedtools getfasta -fi {QUERY} -bed query.bed -s -fo que.fa", shell=True)

# subprocess.run(f"awk '{{print $6,$8,$9}}' {paf}|tr ' ' '\t'  > ref.bed", shell=True)  ##如果参考基因组对应的是正链不进行反转，如果是负链就进行反转
# subprocess.run(f"bedtools getfasta -fi {REF} -bed ref.bed -fo ref.fa", shell=True)



# refsequences = {record.id: record.seq for record in SeqIO.parse('ref.fa', 'fasta')}
# querysequences = {record.id: record.seq for record in SeqIO.parse('que.fa', 'fasta')}

paf=args.paf
# paf="all.paf"
with open(paf, 'r') as infile:
	all_lines = infile.readlines()


with open(args.o, 'w') as outfile:
	outfile.write(f"ref_chr\tref_start\tref_end\tque_chr\tque_start\tque_end\tmatch_base\tmismatch_base\torient\tidentity\talign_len\tmatchindel\trefsou\tquesou\n")
	for line in all_lines:
		columns = line.split('\t')
		parts = re.split(':', columns[0])
		refchr=parts[0]
		parts=re.split('-', parts[1])
		refsource=columns[0]
		refs=int(parts[0])+int(columns[2])-1
		refe=int(parts[0])+int(columns[3])-1
		reflenratio=(refe-refs)/(int(parts[1])-int(parts[0]))
		ori=columns[4]
		quesource=columns[5]
		parts = re.split(':', columns[5])
		quechr=parts[0]
		parts=re.split('-', parts[1])
		ques=int(parts[0])+int(columns[7])-1
		quee=int(parts[0])+int(columns[8])-1
		quelenratio=(quee-ques)/(int(parts[1])-int(parts[0]))
		m=int(columns[9])
		alignlen=int(columns[10])
		for xx in columns:
			if xx.startswith('cg'):
				CIGARstr = xx
				break
		index_of_first_digit = next((index for index, char in enumerate(CIGARstr) if char.isdigit()), None)
		CIGARstr =CIGARstr[index_of_first_digit:]
		matches = re.finditer(r'((?!1[DI])\d+[X])', CIGARstr)
		total_sum = 0
		for match in matches:
			matched_str = match.group(1)
			number_str = matched_str.replace("X", "")
			total_sum += int(number_str)
		mis=total_sum
		iden=m/(m+mis)
		matchindel=m/alignlen
		if refs==ques and refe==quee:
			continue
		outfile.write(f"{quechr}\t{ques}\t{quee}\t{refchr}\t{refs}\t{refe}\t{m}\t{mis}\t{ori}\t{iden}\t{alignlen}\t{matchindel}\t{refsource}\t{quesource}\t{reflenratio}\t{quelenratio}\n")
		#outfile.write(f"{refchr}\t{refs}\t{refe}\t{quechr}\t{ques}\t{quee}\t{m}\t{mis}\t{ori}\t{iden}\t{alignlen}\t{matchindel}\n")
	


# for match in matches:
# 	print(match)

# if(len>=10000):
# 	print(len)
