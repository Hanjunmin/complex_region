
# script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"
# base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"
# DATA_dir="/home/jmhan/project/APG/github_final/PanSD_end/example/"


# sam_list=DATA_dir+"samlis.txt"

# import pandas as pd

# import glob
# import os
# import pickle
# samdf = pd.read_csv(sam_list,sep="\t",header=None)
# refname=list(samdf[samdf[1]=="reference"][0])[0]
# samname=list(samdf[samdf[1]!="reference"][0])
# lines = samname
# name=list(samdf[0])

# rule all:
# 	input:
# 		f"{base_dir}03_refSD/ref_cor/para_gene.OK"
### generate windows
rule gene_windows:
	input:
		f"{base_dir}prepro.OK"
	output:
		f"{base_dir}03_refSD/ref_cor/para_gene.OK"
	shell:
		"""
		mkdir -p "{base_dir}/03_refSD/ref_cor/"
		cd "{base_dir}/03_refSD/ref_cor/"
		python {script_dir}/ref_SD/SD_cordi/scripts/simcountnew.py -fa "{base_dir}01_intersim/refNmask.fa"  -i {base_dir}01_intersim/input.rev ## 1 cores
		python {script_dir}/ref_SD/SD_cordi/scripts/writxt.py ## 30 cores
		touch para_gene.OK
		"""


