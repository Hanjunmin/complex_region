#cd /home/jmhan/project/APG/github_final/PanSD
### snakemake -s PanSD_anno.smk --keep-going --keep-incomplete -j 1
### software
# GFABASE="/home/jmhan/pangenome/software/gfabase"
# RM="/home/Public/software/RepeatMasker-4.1.4/"
# TRF="/home/Public/software/trf-4.09.1/trf"

# script_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/Pan_SD/"
# base_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/run/"
# DATA_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/example/"


VCF=DATA_dir+"example_mc.vcf"
GFA=DATA_dir+"example_mc.gfa"
sam_list=DATA_dir+"samlis.txt"

import pandas as pd
samdf = pd.read_csv(sam_list,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
lines = samname
name=list(samdf[0])


### generate windows
rule make_windows:
	input:
		{sam_list}
	output:
		f"{base_dir}00_preproc/WGS.window.bed"
	shell:
		"""
		cd "{base_dir}/00_preproc/"
		ref=$(cat {sam_list} |grep "reference" |cut -f1)
		samtools faidx {DATA_dir}/FA/${{ref}}.fa
		bedtools makewindows -g <(less -S {DATA_dir}/FA/${{ref}}.fa.fai  |awk 'OFS="\\t"{{print $1,$2}}') -w 1000 -s 500 |awk 'OFS="\\t"{{print $1,$2+1,$3}}'  >WGS.window.bed
		echo "#####################SUCCESS_00-GENERATE THE WINDOWS#####################"
		cd {script_dir}/inter_sim/DB/
		"""

### generate ref_segment_pos
rule gene_prepare:
	input:
		{sam_list}
	output:
		f"{base_dir}00_preproc/reg.posCHM.csv"
	shell:
		'''
		cd "{base_dir}/00_preproc/"
		gfabase load  -o example_mc.gfab {GFA}
		cat  {GFA} |awk '$1=="S"' >reg.S
		cat  {GFA} |awk '$1=="W"' |grep {refname} >ref.W
		{{   echo "{{";  awk '{{if (NR > 1) print prev ","; prev = "\\"" $2 "\\"" ":" "\\"" $3 "\\""}} END {{print prev}}' reg.S;    echo "}}"; }} >reg.json
		python {script_dir}/inter_sim/scripts/gene_pos.py	
		cat > temp_script.R << 'EOF'
		library(data.table)
		library(rjson)
		data<-fread("reg.posCHM.txt")
		node<-fromJSON(file = "reg.json")
		data[, V4 := as.character(V4)]
		data[, fa := node[V4]]
		colnames(data)<-c("chr","start_position","end_position","node","ori","fa")
		fwrite(data,"reg.posCHM.csv",sep="\\t")
EOF
		Rscript temp_script.R
		rm temp_script.R
		'''


rule gene_interseg_homo_var:
		input:
			rules.make_windows.output,
			rules.gene_prepare.output
		output:
			f"{base_dir}01_intersim/new.win"
		shell:
			'''
		    cd "{base_dir}/01_intersim/"
		    grep -v '^#' {VCF} >exam_clear.vcf
		    cat {VCF} |grep "#CHROM" >header
		    ref=$(cat {sam_list} |grep "reference" |cut -f1)
		    python {script_dir}/inter_sim/scripts/vcf_cor.py  --vcf exam_clear.vcf --h header --o reg.vcfcor  --seg {base_dir}/00_preproc/reg.posCHM.csv
		    python {script_dir}/inter_sim/scripts/chm13_sim.py   --w "{base_dir}/00_preproc/WGS.window.bed" --o chm13.sim --r {DATA_dir}/FA/${{ref}}.fa --p {script_dir} ##后面这些脚本内部的相关数据都要处理一下
			bedtools intersect -a "{base_dir}/00_preproc/WGS.window.bed" -b reg.vcfcor  -wa -wb >x
			Rscript -e '
			library(data.table)
			data<-fread("x")
			library(dplyr)
			result <- data %>%
			  group_by(V1, V2, V3) %>%
			  summarise(
			    min_col5 = min(V5),
			    max_col6 = max(V6)
			  )
			result <- result %>%
			  ungroup() %>%  
			  mutate(
			    group_id = group_indices(., min_col5, max_col6)  # 按 min_col5 和 max_col6 重新分组编号
			  )

			new_df <- result %>%
			group_by(group_id) %>%
			slice(1) %>%  
			ungroup()
			fwrite(new_df,"new.win",sep="\\t",col.names=FALSE)
			fwrite(result,"new.win.al",sep="\\t",col.names=FALSE)
			'
			'''


rule gene_interseg_calsim:
		input:
			rules.gene_interseg_homo_var.output
		output:
			# "{base_dir}interseg_homo_sim/reg.sdsim.2k", ##>2k
			# "{base_dir}interseg_homo_sim/reg.out"       ##<2k
			f"{base_dir}01_intersim/sim.out" 
		shell:
			'''
		    cd "{base_dir}/01_intersim/"
			##<2k
			ref=$(cat {sam_list} |grep "reference" |cut -f1)
			cat new.win.al |awk '$5-$4<=2000' >new.win.short
			bedtools coverage -a <(cat {base_dir}/00_preproc/reg.posCHM.csv |awk 'OFS="\\t"{{print $1,$2,$3}}' |tail -n +2) -b <(cat new.win.short |cut -f 1,4,5 |awk 'OFS="\\t"{{print $1,$2-1000,$3+1000}}' |awk 'OFS="\\t"{{$2=($2<0)?1:$2}} 1') |nl -  |awk '$8>0'  |cut -f1  >reg.posid
			awk 'NR==FNR {{lines[$1]; next}} FNR==1 || FNR in lines'  <(cat reg.posid |awk '{{print $1+1}}' )  {base_dir}/00_preproc/reg.posCHM.csv >reg.pos
			python {script_dir}/inter_sim/scripts/seg1_correc_sim_addbord.py --vcf exam_clear.vcf --r {DATA_dir}/FA/${{ref}}.fa --h header --w new.win.short  --o "reg.out" --seg reg.pos --p {script_dir}  ##好像有点慢

			##>2k
			less -S new.win.al  |awk '$5-$4>2000' |cut -f 1,4,5 |awk 'OFS="\\t"{{print $1,$2-1000,$3+1000}}' |awk 'OFS="\\t"{{$2=($2<0)?1:$2}} 1'>x
			less -S new.win.al  |awk '$5-$4>2000' >new.win.long
			cat reg.vcfcor|awk 'OFS="\\t"{{print $1,$2,$3,NR}}' >y
			bedtools coverage -a y -b x >z
			less -S z |awk '$8>0' |cut -f4 >tem
			awk 'NR==FNR {{lines[$1]; next}} FNR in lines'   tem exam_clear.vcf >vcf.extr

			less -S {base_dir}/00_preproc/reg.posCHM.csv |awk 'OFS="\\t"{{print $1,$2,$3,NR}}' |tail -n +2 >y
			bedtools coverage -a y -b x >z
			less -S z |awk '$8>0' |cut -f4 >tem1
			awk 'NR==FNR {{lines[$1]; next}} FNR==1 || FNR in lines'  tem1 {base_dir}/00_preproc/reg.posCHM.csv >node.extr
			python {script_dir}/inter_sim/scripts/sim_correct_2k.py  --o reg.sdsim.2k --h header --r {DATA_dir}/FA/${{ref}}.fa --p {script_dir}
			python {script_dir}/inter_sim/scripts/inte_allsim.py
			'''


rule gene_can_interseghomo:
		input:
			rules.gene_interseg_calsim.output
		output:
			f"{base_dir}01_intersim/input.rev" ,
			f"{base_dir}01_intersim/candidate.merge" 
		shell:
			'''
			cd "{base_dir}/01_intersim/"
			ref=$(cat {sam_list} |grep "reference" |cut -f1)
			FA={DATA_dir}/FA/${{ref}}.fa
			python {script_dir}/inter_sim/scripts/win_simout.py --p  {base_dir} --r ${{FA}} --d {script_dir}/inter_sim/DB/inter_seg_kmer.pkl
			cat simmax.al.add |tail -n +2 |awk '$2>=0.05' |less -S  |tr ":" "\\t" |tr "-" "\\t"  |bedtools sort | bedtools merge -d 35000 >candidate.merge
			bedtools intersect -a  "{base_dir}00_preproc/WGS.window.bed" -b candidate.merge  -wa  >input.rev

			'''


rule GENE_TRF_RM:
		input:
			{VCF}
		output:
			f"{base_dir}01_intersim/rm_trf.bed",
			f"{base_dir}01_intersim/rm_trf_cut1M.bed" ,
			f"{base_dir}01_intersim/refNmask.fa" ,
			f"{base_dir}01_intersim/CHM13v2.fa_rm.bed"
		shell:
			'''
			cd "{base_dir}/01_intersim/"
			ref=$(cat {sam_list} |grep "reference" |cut -f1)
			FA={DATA_dir}/FA/${{ref}}.fa
			RepeatMasker -x -parallel 50 -species human ${{FA}} -e ncbi -dir ref_rm
			cp ./ref_rm/CHM13v2.fa.out ./
			if [ -s CHM13v2.fa.out ]; then
			    RM2Bed.py CHM13v2.fa.out  ##如果没有out文件说明没有TE在这里
			else
			    touch CHM13v2.fa_rm.bed
			fi
			
			trf ${{FA}}  2 7 7 80 10 50 500 -h -l 20 -ngs >trf.out
			python {script_dir}/inter_sim/scripts/trf2bed.py --i trf.out --o trf.bed
			cat CHM13v2.fa_rm.bed <(cat trf.bed 2>/dev/null | tail -n +2) 2>/dev/null | cut -f 1-3 |bedtools sort| bedtools merge > rm_trf.bed  ####?????????
			bedtools maskfasta -fi ${{FA}}  -bed rm_trf.bed -fo refNmask.fa -mc
			cat rm_trf.bed |awk '($3-$2)<1000000{{print $0}}' > rm_trf_cut1M.bed

			'''

rule clear_inter_sim:
	input:
		rules.gene_can_interseghomo.output,
		rules.GENE_TRF_RM.output
	output:
		f"{base_dir}prepro.OK" 
	shell:
		'''
		cd "{base_dir}/01_intersim/"
		find . -maxdepth 1 -type f ! -name 'rm_trf_cut1M.bed' ! -name 'trf.bed'  ! -name 'rm_trf.bed'  ! -name 'CHM13v2.fa_rm.bed' ! -name 'refNmask.fa' ! -name 'input.rev' ! -name 'candidate.merge'  ! -name 'sim.out' -delete
		touch {output}

		'''







# python /home/jmhan/project/APG/github_final/PanSD/script/reference_SD/simcountnew.py -fa refNmask.fa -i input.rev ## 1 cores
# python /home/jmhan/project/APG/github_final/PanSD/script/reference_SD/writxt.py ## 30 cores


# snakemake -s /home/jmhan/project/APG/github_final/PanSD/script/reference_SD/SD.smk  -j 30 ## 2min
# fold="/home/jmhan/project/APG/github_final/PanSD/output/reference/"
# cp ./filt/beforegroup.txt ./
# ##是不是要加一下split文件夹新建
# Rscript /home/jmhan/project/APG/github_final/PanSD/script/reference_SD/1.group.r $fold"/filt/split/" $fold  ##1min
# cd $fold"/filt/split/"
# snakemake -s /home/jmhan/project/APG/github_final/PanSD/script/reference_SD/SD1.smk  -j 30 ###注意一下至少需要300G的内存

# snakemake -s /home/jmhan/project/APG/github_final/PanSD/script/reference_SD/endprocess.smk -j 1


