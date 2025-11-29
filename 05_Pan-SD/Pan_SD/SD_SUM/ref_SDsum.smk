
# script_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/Pan_SD/"
# base_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/run/"
# DATA_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/example/"



sam_list=DATA_dir+"samlis.txt"

import pandas as pd

import glob
import os
import pickle
samdf = pd.read_csv(sam_list,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
lines = samname
name=list(samdf[0])
ALLSAM=name

FA=DATA_dir+"/FA/"+refname+".fa"
repeatmasker="/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/asm_repeatmasker.out.bed" ## /home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/asm_repeatmasker.out.bed
split_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/minimapsplit")
minimap_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/minimapspnow")
allclus_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/endfilt")
fasta_file1=FA
fasta_file2=FA

# rule all:
# 	input:
# 		expand("{base_dir}02_refSD/sam_cor/2_ALIGN/{line}/{line}.ref.SD.al", base_dir=base_dir,line=lines),
# 		expand("{base_dir}03_nonref/02_ITER/2_SD_GENE/{sam}/ALIGN/{sam}.ITER.SD.al", base_dir=base_dir,sam=ALLSAM),
# 		f"{base_dir}ref.SD.al",
# 		f"{base_dir}iter.SD.al"


# rule gene_refSD_sta:
# 	input:
# 		edge="{base_dir}02_refSD/sam_cor/1_SD_CANPATH/sd.edge",
# 		# link="{base_dir}02_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/all.end.statistics.filt",  
# 		graphsim="{base_dir}02_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/T2Thave.bed",  
# 		# graphsim2="{base_dir}02_refSD/sam_cor/2_ALIGN/{line}/processnosam.txt"
# 	output:
# 		"{base_dir}02_refSD/sam_cor/2_ALIGN/{line}/{line}.ref.SD.al"
# 	shell:
# 		'''
# 		cd {base_dir}02_refSD/sam_cor/2_ALIGN/{wildcards.line}
# 		cat process.txt end.txt | grep -v "sam" >sim.txt
# 		cat {wildcards.line}.sd | tail -n +2 |awk 'OFS="\\t"{{print $2,$9":"$10+$3"-"$10+$4,$5}}' >pos.dic
# 		Rscript {script_dir}/SD_SUM/scripts/t2t_mer_al.r  {wildcards.line} {base_dir}/02_refSD/sam_cor/1_SD_CANPATH/sd.edge
# 		'''

rule merge_refSD:
	input:
		expand("{base_dir}02_refSD/sam_cor/2_ALIGN/{line}/{line}.ref.SD.al", base_dir=base_dir,line=lines)
	output:
		f"{base_dir}ref.SD.al"
	run:
		import pandas as pd
		import glob
		input_files = list(input)
		print(input_files)
		print(output)
		dfs = []
		for file in input_files:
			print(file)
			df = pd.read_csv(file, sep='\t', header=0)
			dfs.append(df)
		merged_df = pd.concat(dfs, axis=0)
		merged_df=merged_df[['chr1_init', 'chr2_init', 'chr1_source', 'chr2_source','chr2_init.sim','chr1_init.sim','chr1', 'start1', 'end1', 'chr2', 'start2', 'end2','match_base', 'mismatch_base', 'orient', 'identity', 'align_len','matchindel', 'SAM']]
		merged_df.to_csv("ref.SD.al", sep='\t', index=False)


# rule gene_nonrefSD_sta:
# 	input:
# 		# ALPAIR="{base_dir}03_nonref/02_ITER/2_SD_GENE/tranall.uniq",
# 		INIT= "{base_dir}/03_nonref/non_ref.OK",
# 		OUT="{base_dir}03_nonref/02_ITER/2_SD_GENE/{sam}/ALIGN/{sam}.eebed.rev.out",
# 		CAN="{base_dir}03_nonref/02_ITER/2_SD_GENE/{sam}/ALIGN/{sam}.all1.end.statistics",
# 	output:
# 		"{base_dir}03_nonref/02_ITER/2_SD_GENE/{sam}/ALIGN/{sam}.ITER.SD.al"
# 	shell:
# 		'''
# 		cd {base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.sam}/ALIGN/
# 		python {script_dir}/SD_SUM/scripts/trans_pair2.py --pos "{base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.sam}/ALIGN/{wildcards.sam}.pos" --W {wildcards.sam}.init.W --o {wildcards.sam} --pair "{base_dir}03_nonref/02_ITER/2_SD_GENE/tranall.uniq"
# 		python {script_dir}/SD_SUM/scripts/construct_dict.py --pair "{base_dir}03_nonref/02_ITER/2_SD_GENE/tranall.uniq" --sam {wildcards.sam} --scp {script_dir}/ref_SD/sam_sd/scripts/filter_sam_SD.r
# 		'''



rule merge_nonrefSD:
	input:
		expand("{base_dir}03_nonref/02_ITER/2_SD_GENE/{sam}/ALIGN/{sam}.ITER.SD.al", base_dir=base_dir,sam=ALLSAM)
	output:
		f"{base_dir}iter.SD.al"
	run:
		import pandas as pd
		import glob
		input_files = list(input)
		print(input_files)
		print(output)
		dfs = []
		for file in input_files:
			print(file)
			df = pd.read_csv(file, sep='\t', header=0)
			dfs.append(df)
		merged_df = pd.concat(dfs, axis=0)
		merged_df=merged_df[['chr1_init', 'chr2_init', 'chr1_source', 'chr2_source','chr2_init.sim','chr1_init.sim','chr1', 'start1', 'end1', 'chr2', 'start2', 'end2','match_base', 'mismatch_base', 'orient', 'identity', 'align_len','matchindel', 'SAM']]
		merged_df = merged_df.fillna('.')
		merged_df.to_csv("iter.SD.al", sep='\t', index=False)


