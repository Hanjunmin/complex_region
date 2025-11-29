import pandas as pd
import subprocess
import re as re
import subprocess
import random
import math
import argparse
parser = argparse.ArgumentParser(description='CIGAR-CALCULATE')
parser.add_argument('--i', help='input file')
parser.add_argument('--o', help='output file')
parser.add_argument('--fa', help='REF fasta')
args = parser.parse_args()
#FA="/dssg/home/acct-clsmyf/clsmyf-user1/data/chm13v2.0.fa"
FA=args.fa

def cal_iden(data):
	columns = data
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
	m=int(columns[9])
	alignlen=int(columns[10])
	mis=total_sum
	iden=m/(m+mis)
	matchindel=m/alignlen
	parts = re.split('[:-]', columns[0])
	refs=int(parts[1])+int(columns[2])-1
	refe=int(parts[1])+int(columns[3])-1
	ori=columns[4]
	parts = re.split('[:-]', columns[5])
	ques=int(parts[1])+int(columns[7])-1
	quee=int(parts[1])+int(columns[8])-1
	return iden,alignlen,matchindel,m,mis,refs,refe,ques,quee,ori


#df = pd.read_csv(,sep=",",header=None,keep_default_na=False,na_values=[])
# with open("end.txt", 'r') as f:
# 	lines = f.readlines()


### try一下模拟退火-------------------------------------------------try start
import subprocess
import random
import math



def minimap_align(rchr, rs, ree, qchr, qs, qe, next_situation):
	shifted_rs = int(rs) + next_situation[0]
	shifted_re = int(ree) + next_situation[1]
	shifted_qs = int(qs) + next_situation[2]
	shifted_qe = int(qe) + next_situation[3]
	command = f"minimap2 -c --eqx <(samtools faidx {FA} {rchr}:{shifted_rs}-{shifted_re}) <(samtools faidx {FA} {qchr}:{shifted_qs}-{shifted_qe}) -A 5 -B 4 -O 40 -E 1 -s 3000"
	result = subprocess.run(command, shell=True, executable='/bin/bash', capture_output=True, text=True)
	align = result.stdout
	if align == "":
		return None
	alignlist = align.strip().split('\t')
	e = cal_iden(alignlist)
	return e



with open(args.i, 'r') as f:
	lines = f.readlines()


kkk=0
output_file=[]
output_id=[]
with open(args.o, 'w') as f:
	for line in lines:
		kkk=kkk+1
		line = line.strip().split('\t')
		print(line)
		best_identity = float(line[9])
		found = False
		step_list=[5000,4000,3000,2000,1000,500,200,100,0]
		step_list=[5000,4000,3000,2000,1000,500,200,100,0]
		target_value=int(line[6])
		step_listnew=step_list[step_list.index(max([x for x in step_list if x <= target_value])):]
		step_listnew.reverse()
		for step in range(len(step_listnew)-1):
			exten=random.randint(step_listnew[step],step_listnew[step+1])
			allsituation=[[-exten,0,0,0],[-exten,+exten,0,0],[0,+exten,0,0],[+exten,0,0,0],[0,-exten,0,0],[+exten,-exten,0,0],
						  [0,0,-exten,0],[0,0,-exten,+exten],[0,0,+exten,0],[0,0,+exten,0],[0,0,0,-exten],[0,0,+exten,-exten],
			              [-exten,0,-exten,0],[-exten,+exten,-exten,+exten],[0,+exten,0,+exten],[+exten,0,+exten,0],[0,-exten,0,-exten],[+exten,-exten,+exten,-exten]]
			for next_situation in allsituation:
				rchr = line[0]
				rs = int(line[1])
				ree = int(line[2])
				qchr = line[3]
				qs = int(line[4])
				qe = int(line[5])
				e = minimap_align(rchr, rs, ree, qchr, qs, qe, next_situation)
				if e is None:
					continue
				delta_identity = e[0] - best_identity
				if delta_identity >= 0:
					best_identity = e[0]
					rs = e[5]
					ree = e[6]
					qs = e[7]
					qe = e[8]
					if e[0] >= 0.9 and e[1]>=1000:
						print(e)
						output_line = "\t".join(line + list(map(str, e))) + "\n"
						f.write(output_line)
						output_file.append(output_line)
						output_id.append(line[0])
						found = True
						break
				if found:
					break
			if found:
				break
