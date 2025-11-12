### snakemake -s index_construct_anno.smk --keep-going --keep-incomplete -j 1 -R
### software
# REPEATMASKER=config["software"]["RM"]
#REPEATMASKER="/home/Public/software/RepeatMasker"
# singularity_image=config["singularity"]["image"]
# singularity_args=config["singularity"]["args"]

# DOCKER=config["env"]
# TRF=config["software"]["TRF"]
# GFABASE=config["software"]["gfabase"]


base_dir=config['base_dir']
script_dir=config['script_dir']
data_dir=config['data_dir']
VCF=data_dir+"1p36.vcf"
GFA=data_dir+"1p36.gfa"
BED=data_dir+"1p36.bed" 
sam_list=data_dir+"samlis.txt"
REF=config['reference']

import pandas as pd
samdf = pd.read_csv(sam_list,sep="\t",header=None)
beddf = pd.read_csv(BED,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
name=list(samdf[0])
allchrom=list(beddf[0])
lines = name



# rule all:
# 	input:
# 		f"{base_dir}02_index_construction/seq_div/entropy_final.txt",
# 		f"{base_dir}02_index_construction/seq_div/simal",
# 		f"{base_dir}02_index_construction/seq_div/seq_div_gene.OK",
# 		f"{base_dir}02_index_construction/inter_seg/interseg_gene.OK",
# 		# f"{base_dir}02_index_construction/bub/bub_gene.OK",  ###基于bubble的最终结果
# 		# expand("{base_dir}02_index_construction/TE_div/{chrom}/ref.anchor.samal",base_dir=base_dir,chrom=allchrom),
# 		f"{base_dir}02_index_construction/TE_div/TE_gene.OK",
# 		f"{base_dir}02_index_construction/intra_seg/intra_gene.OK"


rule gene_seq_div:
	input:
		vcf={VCF},
		window=f"{base_dir}01_preproc/cor_window.bed",
		entroyin=f"{base_dir}01_preproc/conti_path/reg.vcf_corr_dot.en",
		pre=f"{base_dir}01_preproc/gene_prep.OK"
	output:
		f"{base_dir}02_index_construction/seq_div/entropy_final.txt",
		f"{base_dir}02_index_construction/seq_div/simal"
	# singularity_args: {singularity_args}
	shell:
		"""
		mkdir -p "{base_dir}02_index_construction/seq_div"
		cd "{base_dir}02_index_construction/seq_div"
		cat {input.vcf} |grep -v '^#' | awk -F '\\t' '{{OFS='\\t';maxlen=0;split($5,spres,",");same_vals = 1
		for(i in spres){{
		len=length(spres[i])
		if(len>maxlen) maxlen=len
		}}
		print $1"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"length($4)"\\t"maxlen
		}}' >spe.vcf.tem  
		cat {input.entroyin} |tail -n +2 >py_entrin.txt
		[ -f py_entrin.out ] && rm py_entrin.out
		python {script_dir}/index_construction/scripts/seq_div/entro_cal.py
		# [ -f py_entrin.txt ] && rm py_entrin.txt
		less -S spe.vcf.tem  |awk 'OFS="\\t"{{print $1,$2,$2+$6,$7}}' >lef.bed.tem
		paste <(less -S lef.bed.tem) <(cat py_entrin.out |awk '{{print $NF}}') >entropy_bed.txt
		bedtools intersect -a {input.window} -b entropy_bed.txt -wa -wb >x.tem
		Rscript -e '
		library(data.table)
		data<-fread("x.tem")
		library(dplyr)
		library(data.table)
		setDT(data)
		result <- data[, Sum_group := sum(V7), by = .(V1, V2, V3)]
		result$weight=result$V7/result$Sum_group
		result$weightentro=result$weight*result$V8
		end <- result[, endentro := sum(weightentro), by = .(V1, V2, V3)]
		fwrite(end,"testAL.temp",sep="\\t")
		fwrite(end[,c("V1","V2","V3","endentro")],"test.tem",sep="\\t",col.names = FALSE) '
		cat test.tem |uniq >entropy_final.txt


		less -S spe.vcf.tem  |awk 'OFS="\\t"{{print $1,$2,$2+$6,$6,$7}}' >1
		bedtools intersect -a {input.window}  -b 1  -wa -wb >2
		cat > temp_script.R << 'EOF'
		library(data.table)
		data<-fread("2")
		data$newwindos=data$V2
		data$newwindoe=data$V3
		data[data$V5<data$V2,]$newwindos=data[data$V5<data$V2,]$V5
		data[data$V6>data$V3,]$newwindoe=data[data$V6>data$V3,]$V6
		data$maxal=pmax(data$V7,data$V8)
		setDT(data)

		result <- data[, .(
		  alignlen   = sum(maxal),
		  diff       = sum(V7),
		  ns         = min(newwindos),
		  ne         = max(newwindoe),
		  range_flat = list(unlist(mapply(`:`, V5, V6))),
		  range_len  = length(unique(unlist(mapply(`:`, V5, V6)))) 
		), by = .(V1, V2, V3)]


		result[, nlen := length(unique(unlist(range_flat))), by = .(V1, V2, V3)]

		result$match=(result$ne-result$ns+1)-result$nlen
		result$sim=result$match/(result$match+result$alignlen)
		fwrite(result[,c("V1","V2","V3","sim")],"simal",sep="\\t")
EOF
		
		Rscript temp_script.R

		echo "#####################SUCCESS_01-GENERATE THE SEQ_DIV INDEX#####################"
		"""


rule clear_seqdiv_fold: 
		input:
			rules.gene_seq_div.output
		output:
			f"{base_dir}02_index_construction/seq_div/seq_div_gene.OK"
		shell:
			"""
			cd "{base_dir}02_index_construction/seq_div/"
			find . -maxdepth 1 -type f ! -name 'simal' ! -name 'entropy_final.txt' -delete
			touch "seq_div_gene.OK"
			"""


rule gene_interseg_homo_var:
		input:
			pre=f"{base_dir}01_preproc/gene_prep.OK"
		output:
			f"{base_dir}02_index_construction/inter_seg/new.win"
		shell:
			"""
		    cd "{base_dir}02_index_construction/inter_seg/"
		    less -S {VCF}  |grep -v "^#" >reg.vcf
			cat {VCF} | grep "#CHROM" >headern
		    python {script_dir}/preproc/scripts/vcf_cor.py --vcf reg.vcf   --o reg.vcfcor  --seg "{base_dir}01_preproc/reg.posCHM.txt"
		    python {script_dir}/index_construction/scripts/seg_homo/ref_sim.py --w "{base_dir}01_preproc/cor_window.bed"    --d {script_dir}/index_construction/scripts/seg_homo/DB/ --o chm13.sim --p 20 --r {REF} ##后面这些脚本内部的相关数据都要处理一下
			bedtools intersect -a "{base_dir}01_preproc/cor_window.bed" -b reg.vcfcor  -wa -wb >x
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
			"""


rule gene_interseg_calsim:
		input:
			rules.gene_interseg_homo_var.output
		output:
			f"{base_dir}02_index_construction/inter_seg/reg.sdsim.2k", ##>2k
			f"{base_dir}02_index_construction/inter_seg/reg.out"       ##<2k
		shell:
			"""
		    cd "{base_dir}02_index_construction/inter_seg/"
			##<2k
			cat new.win.al |awk '$5-$4<=2000' >new.win.short
			bedtools coverage -a <(cat "{base_dir}01_preproc/reg.posCHM.csv" |awk 'OFS="\\t"{{print $1,$2,$3}}' |tail -n +2) -b <(cat new.win.short |cut -f 1,4,5 |awk 'OFS="\\t"{{print $1,$2-1000,$3+1000}}' |awk 'OFS="\\t"{{$2=($2<0)?1:$2}} 1') |nl -  |awk '$8>0'  |cut -f1  >reg.posid
			awk 'NR==FNR {{lines[$1]; next}} FNR==1 || FNR in lines'  <(cat reg.posid |awk '{{print $1+1}}' )  "{base_dir}01_preproc/reg.posCHM.csv" >reg.pos
			python {script_dir}/index_construction/scripts/seg_homo/sim_correct_s2k.py   --vcf reg.vcf  --w new.win.short  --o "reg.out" --seg reg.pos --d1 {script_dir}/index_construction/scripts/seg_homo/DB/ --d2 {base_dir}01_preproc/seg_db --p 20 --r {REF}

			##>2k
			less -S new.win.al  |awk '$5-$4>2000' |cut -f 1,4,5 |awk 'OFS="\\t"{{print $1,$2-1000,$3+1000}}' |awk 'OFS="\\t"{{$2=($2<0)?1:$2}} 1'>x
			less -S new.win.al  |awk '$5-$4>2000' >new.win.long
			cat reg.vcfcor|awk 'OFS="\\t"{{print $1,$2,$3,NR}}' >y
			bedtools coverage -a y -b x >z
			less -S z |awk '$8>0' |cut -f4 >tem
			awk 'NR==FNR {{lines[$1]; next}} FNR in lines'   tem reg.vcf >vcf.extr

			less -S "{base_dir}01_preproc/reg.posCHM.csv" |awk 'OFS="\\t"{{print $1,$2,$3,NR}}' |tail -n +2 >y
			bedtools coverage -a y -b x >z
			less -S z |awk '$8>0' |cut -f4 >tem1
			awk 'NR==FNR {{lines[$1]; next}} FNR==1 || FNR in lines'  tem1 "{base_dir}01_preproc/reg.posCHM.csv" >node.extr
			python {script_dir}/index_construction/scripts/seg_homo/sim_correct_l2k.py  --o reg.sdsim.2k --d1 {script_dir}/index_construction/scripts/seg_homo/DB/ --d2 {base_dir}01_preproc/seg_db --p 20 --r {REF}

			"""

rule gene_interseg_procesdot:
		input:
			rules.gene_interseg_calsim.output
		output:
		    f"{base_dir}02_index_construction/inter_seg/sim.end_process.dot" 
		shell:
			"""
			cd "{base_dir}02_index_construction/inter_seg/"
			python {script_dir}/index_construction/scripts/seg_homo/inte_sim.py
			#### process the dot of file
			mkdir -p "{base_dir}02_index_construction/inter_seg/refin_dot"
			find "{base_dir}02_index_construction/inter_seg/refin_dot" -mindepth 1 -delete
			cd "{base_dir}02_index_construction/inter_seg/refin_dot"
			cat > temp_script.R << 'EOF'
			library(data.table)
			options(scipen=999)
			options(digits = 5)
			data<-fread(paste("../sim.out",sep=""))
			colname=colnames(data)
			colname=colname[c(-1,-543)] ##modify
			for(sam in colname){{
			        print(sam)
			        samlis=c("win",sam)
			        samdf=data[,..samlis]
			        samdf[, original_order := .I]
			        a <- samdf[samdf[[sam]] != "." & samdf[[sam]] != "", ]
			        a[[sam]] <- as.numeric(a[[sam]])
			        b <- samdf[!(samdf[[sam]] != "." & samdf[[sam]] != ""), ]
			        samdf <- rbind(a, b)[order(original_order), ]
			        samdf[, original_order := NULL]
			        fwrite(samdf,paste(sam,".bed",sep=""),sep="\\t")
			        scipen = 999
			}}

			column_names <- data.frame(colname)
			fwrite(column_names, file = "allsamname", col.names = FALSE)
EOF
		Rscript temp_script.R
		while IFS= read -r sam; do
		    samfile=$sam
		    echo $samfile
		    if [ "$sam" = "T2TCHM13" ]; then
		    samfile=={refname}
		    cat "$sam".bed |cut -f2 >$samfile".out"
		  else
		    cat {base_dir}01_preproc/conti_path/$samfile"_REFcor.pos" |cut -f2,9,10 |tail -n +2 >$samfile".temp"
		    bedtools coverage -a <(cat "$sam".bed |tr ":" "\\t" |tr "-" "\\t" |tail -n +2 ) -b $samfile".temp"  |awk '{{OFS="\\t"}}$4=="."{{$4="0"}}1' |awk '{{OFS="\t"}}$8==0{{$4="."}}1' |less -S  |awk 'OFS="\t"{{print $4}}' |awk -v sam=$sam 'BEGIN{{
		print sam }}1'  >$samfile".out"
		        fi
		done < allsamname
		paste <(cat "../sim.out" |cut -f 543) *.out  >sim.end_process.dot
		cd "{base_dir}02_index_construction/inter_seg/"
		cp ./refin_dot/sim.end_process.dot ./
		"""


rule gene_interseg_inte:
		input:
			rules.gene_interseg_procesdot.output
		output:
			f"{base_dir}02_index_construction/inter_seg/reg.simout" 
		shell:
			"""
			cd "{base_dir}02_index_construction/inter_seg/"
			cat > temp_script.R << 'EOF'
			library(dplyr)
			library(data.table)

			dfal <- data.frame()
			trimmed_var <- function(x, prop = 0.05) {{
			  x=x[!is.na(x)]
			  lower <- quantile(x, probs = prop, na.rm = TRUE)
			  upper <- quantile(x, probs = 1 - prop, na.rm = TRUE)
			  x_trimmed <- x[x >= lower & x <= upper]
			  if (length(x_trimmed) < 2) {{
			    NA
			  }} else {{
			    sd(x_trimmed)
			  }}
			}}
			library(data.table)
			library(dplyr)
			library(foreach)
			df <- fread("sim.end_process.dot")
			  win1 <- df$win
			  df <- df %>% mutate(across(everything(), as.numeric))
			  trimsd <- apply(df, 1, trimmed_var)
			  df <- df %>% mutate(sdrow_max = pmax(!!!select(., -c("win")), na.rm = TRUE))

			  df %>%
			    mutate(win = win1, sdtrimsd = trimsd) %>%
			    select(win, sdrow_max, sdtrimsd) %>%
			    fwrite(paste0("reg.simout"), sep = "\\t")

EOF
			Rscript temp_script.R
			"""


rule clear_interseg_fold: 
	input:
		rules.gene_interseg_inte.output
	output:
		f"{base_dir}02_index_construction/inter_seg/interseg_gene.OK"
	shell:
		"""
		cd "{base_dir}02_index_construction/inter_seg/"
		find . -maxdepth 1 -type f ! -name 'reg.simout' -delete
		[ -d "refin_dot" ] && rm -rf "refin_dot"
		touch "interseg_gene.OK"
		"""

rule gene_seq_bub:
	input:
		winin=f"{base_dir}01_preproc/cor_window.bed",
		pre=f"{base_dir}01_preproc/gene_prep.OK"
	output:
		f"{base_dir}02_index_construction/bub/WGS_bubble_APG.txt"
	shell:
		"""
		cd "{base_dir}02_index_construction/bub/"
		gfabase load -o my.gfab {GFA}
		bash {script_dir}/index_construction/scripts/bub/para.sh  "{base_dir}01_preproc/cor_window.bed" my.gfab gfabase  {refname}
		"""



rule gene_bub_path:
		input:
			rules.gene_seq_bub.output
		output:
			f"{base_dir}02_index_construction/bub/gene_path.OK"
		shell:
			"""
			cd "{base_dir}02_index_construction/bub/"
			cat > temp_script.R << 'EOF'
			library(data.table)
			data<-fread("WGS_bubble_APG.txt")
			library(dplyr)
			library(stringr)
			result <- data %>%
			  group_by(V4, V5, V1) %>%
			  summarise(
			    min_col1 = min(V2),
			    max_col3 = max(V3),
			    combined_V6 = paste(V6, collapse = ", "),  # 将 V6 列的值连接成一个字符串
			    .groups = 'drop'
			  )

			res <- result %>%
			  filter(str_detect(combined_V6, "extend"))
			res=res[,c(3,4,5)]
			fwrite(res,"bub.extend",sep="\\t")
EOF
			Rscript temp_script.R
			cat bub.extend |tail -n +2 |bedtools sort |bedtools merge >pathin
			[ -f Sregional.node ] && rm Sregional.node
			while read -r chr start end; do
			  gfabase sub my.gfab {refname}#0#$chr":"$start"-"$end --range --view --cutpoints 1 --guess-ranges -o $chr":"$start"-"$end.gfa
			  reg=$chr":"$start"-"$end
			  cat $chr":"$start"-"$end.gfa  |awk '$1=="S"' |less -S |cut -f 1-2 |awk  -v  reg=$reg 'OFS="\\t"{{print $1,$2,reg}}' |sort -k2,2n |uniq >>Sregional.node
			done < pathin
			python {script_dir}/index_construction/scripts/bub/subregion.py --pr  "{base_dir}01_preproc/Wpath/W.al"  ####--pr ./graph/2M/reg.path
			less -S pathin |awk '$3-$2<1000000' >al.subregional.1M
			[ -f CHM13.path ] && rm CHM13.path
			while read -r subchr substart subend; do
			  x=$subchr":"$substart"-"$subend
			  less -S ${{x}}.gfa |grep {refname} |grep "^S" |cut -f2,8 |sed "s/,//g" |sed "s/#/\\t/g" |cut -f 1,4 |tr ":" "\\t" |tr "-" "\\t" |sort -k3,3n |cut -f1|tr "\\n" "," |sed "s/,/+,/g" |awk -v reg=$x 'OFS="\\t"{{print "P",reg"@CHM13",$0,"*"}}' >>CHM13.path
			done < pathin
			python {script_dir}/index_construction/scripts/bub/high_freseg.py
			touch gene_path.OK
			"""

rule gene_bub_segnum:
	input:
		rules.gene_bub_path.output
	output:
		f"{base_dir}02_index_construction/bub/gene_segnum.OK"
	shell:
		"""
		cd "{base_dir}02_index_construction/bub/"
		while read -r chr start end; do
		  segm=$chr":"$start"-"$end
		  len=$(echo -e "$end-$start" | bc)
		  cat ${{segm}}.gfa |grep '^S' |cut -f 2-3 >"temp."${{segm}}.S.start
		  python {script_dir}/index_construction/scripts/bub/filter_path.py --i ${{segm}}.path.dou.filt ##30 cores based on similarity  ${{segm}}.path.dou.filt based on the complement of paths
		  cat ${{segm}}.path.dou.filt.filt  >"temp."${{segm}}.allpath
		  # python /home/jmhan/project/APG/github_final/01_1p36/index_construction/script/bub/POLY_complex_dupseg.py --sim 0.9 --mb 50 --dell 200  --region ${{segm}} ###普通模式
		  python {script_dir}/index_construction/scripts/bub/para.py --sim 0.9 --mb 50 --dell 200  --region ${{segm}} ###并行模式
		  cat ${{segm}}.gfa >end.gfa
		  cat ${{segm}}.gfa NODE.P.${{segm}}data >${{segm}}.test.gfa  || true
		  cat <(cat ${{segm}}.gfa |head -n 1) S.${{segm}}.data L.${{segm}}.data P.${{segm}}.data >${{segm}}.sim.gfa  || true
		  python {script_dir}/index_construction/scripts/bub/simnode_inte.py --reg ${{segm}}
		  python {script_dir}/index_construction/scripts/bub/poly_ratio.py --reg ${{segm}}
		  [ -f "temp."${{segm}}.S.start ] && rm "temp."${{segm}}.S.start
		  [ -f "temp."${{segm}}.allpath ] && rm "temp."${{segm}}.allpath

		done < pathin
		touch gene_segnum.OK


		"""


rule sum_bub_feature:
	input:
		rules.gene_bub_segnum.output
	output:
		f"{base_dir}02_index_construction/bub/bub.txt"
	shell:
		'''
		cd "{base_dir}02_index_construction/bub/"
		ls polycomplex* >reg.pm.bed
		[ -f reg.end.sta ] && rm reg.end.sta
		cat reg.pm.bed  |sed "s/polycomplex//g" |sed 's/.txt//g' >reg.pathregion.bed
		less -S reg.pathregion.bed |cut -f2 |tr ":" "\\t" |tr "-" "\\t" >reg.pathregion.se.bed
		paste reg.pm.bed reg.pathregion.bed reg.pathregion.se.bed>reg.pm_region.bed
		while read -r subfold subregion chr start end; do
		poly=$(cat $subfold)
		nodenum=$(cat  ${{subregion}}.sim.gfa|awk '$1=="S"'  |wc -l)
		intenode=$(cat simnode${{subregion}}.txt|cut -f1)
		polyratio=$(cat polyratio${{subregion}}.txt|cut -f1)
		polyratio=$(printf "%.20f" "$polyratio")
		subpath=$(cat ${{subregion}}.path.dou.filt|wc -l)
		allpath=$(cat ${{subregion}}.path.dou.bef|wc -l)
		pathdif=$(echo "scale=2; $subpath / $allpath" | bc)
		sumnode=$(cat ${{subregion}}.gfa  |grep '^S' |cut -f 2-3|  awk 'length($2) != 1'  |awk '{{sum += length($2)}} END {{print sum}}')
		truepos=$(echo "$end - $start" | bc)
		min_val=$(echo "if ($sumnode < $truepos) $sumnode else $truepos" | bc)
		max_val=$(echo "if ($sumnode >= $truepos) $sumnode else $truepos" | bc)
		refinelen=$(echo "scale=2; $min_val / $max_val" | bc)
		printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\n" "$chr" "$start" "$end" "$poly" "$polyratio" "$intenode" "$pathdif" "$refinelen" ${{subregion}}>> reg.end.sta
		done<reg.pm_region.bed
		bedtools intersect -b reg.end.sta -a "{base_dir}01_preproc/cor_window.bed"  -f 1 -wa  -wb >x
		cat > temp_script.R << 'EOF'
		library(data.table) ###原先的加权重
		data <- read.table("x", header = FALSE, sep = "\t")
		sigmoid <- function(x, k, c) {{
		  1 / (1 + exp(-k * (x - c)))
		}}
		# data$lenrefine=1-sigmoid(data$V6,300,0.85)
		data$polyrefine=sigmoid(data$V8,100,0.1) ###polyratio
		data$H=0
		data$H <- ifelse(data$V10 < 0.5, 1, 
		                ifelse(data$V10 > 0.5, 0, data$H)) ##path相似性
		data$V4new=2*sigmoid(data$polyrefine*data$V7,3,3) ###V7 poly
		data$V5new=sigmoid(data$V9,3,3) ####simplenode
		data$new=data$V4new+data$V5new
		data$new=data$new+1*data$H
		data <- na.omit(data)
		data$newnorm=(data$new-min(data$new))/(max(data$new)-min(data$new))
		fwrite(data,"bub.txt",sep="\\t",col.names=FALSE)
EOF
		Rscript temp_script.R
	'''



rule clear_bub_fold: 
	input:
		rules.sum_bub_feature.output
	output:
		f"{base_dir}02_index_construction/bub/bub_gene.OK"
	shell:
		"""
		cd "{base_dir}02_index_construction/bub/"
		find . -maxdepth 1 -type f ! -name 'bub.txt' -delete
		touch "bub_gene.OK"
		"""



rule gene_seq_TEpos:
	input:
		f"{base_dir}01_preproc/cor_window.bed",
		f"{base_dir}01_preproc/gene_prep.OK"
	output:
		"{base_dir}/02_index_construction/TE_div/{chrom}/ref.anchor.samal"
	shell:
		'''
		mkdir -p "{base_dir}02_index_construction/TE_div/{wildcards.chrom}"
		cd "{base_dir}02_index_construction/TE_div/{wildcards.chrom}"
		python {script_dir}/index_construction/scripts/te_div/graph_TE_addDEL.py --CHR {wildcards.chrom} --T2Tpos "{base_dir}01_preproc/reg.posCHM.txt" --s 0 --e 101 --r {REF} --w {base_dir}01_preproc/Wpath/ --json {base_dir}01_preproc/reg.json --p 20
		cat *.csv |grep "@" |grep -v "DEL" >seq.in
		cat *.csv |grep "@" |grep "DEL" >DEL.in
		python {script_dir}/index_construction/scripts/te_div/gene_seq.py
		cat output.fasta| sed "s/seq/\\n/g" >output_n.fasta
		samtools faidx output_n.fasta
		
		RepeatMasker -x -parallel 50 -species human output_n.fasta -e ncbi
		RM2Bed.py output_n.fasta.out

		bedtools getfasta -fi {REF} -bed {BED} -fo reg.fa
		
		RepeatMasker -x -parallel 50 -species human reg.fa -e ncbi
		RM2Bed.py reg.fa.out
		cat reg.fa_rm.bed | awk 'OFS="\\t"{{gsub(/[:,-]/,"\\t",$1); print}}' |awk 'OFS="\\t"{{print $1,$4+$2,$5+$2,$6,$7,$8,$9,$10,$11,$12}}' >reg.fa_rm.bed.END
		cat reg.fa_rm.bed.END | grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.rm
		cat output_n.fasta |grep ">" >seg.temid
		python {script_dir}/index_construction/scripts/te_div/inte_seq.py --chr {wildcards.chrom} --json {base_dir}01_preproc/reg.json
		cat reg.fa_rm.out |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >reg.end.out
		
		cat *fa_rm.alseq >reg.end.alSDR ###全部的SDR

		less -S reg.end.alSDR   |cut -f 1-3 |sort |uniq >ref.anchor
		less -S reg.end.alSDR   |cut -f 1-4 |sort |uniq >ref.anchor.sam
		cat DEL.in   |tr "@" "\\t" |awk -v chr={wildcards.chrom} 'OFS="\\t"{{print chr,$8,$9,$2}}'  |sed "s/\.0//g"  >DEL.in.sam
		cat DEL.in.sam |cut -f 1-3 >DEL.in
		bedtools intersect -b DEL.in.sam -a reg.fa_rm.bed.END -wa -wb  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.anchor.rm.del  || true
		bedtools intersect -b ref.anchor -a reg.fa_rm.bed.END -wa -wb -f 0.2  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.anchor.rm || true
		cat DEL.in ref.anchor |sort |uniq >ref.anchoral || true
		bedtools intersect -b ref.anchoral -a reg.fa_rm.bed.END -wa -wb -f 0.2  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.anchoral.rm || true
		cat DEL.in.sam ref.anchor.sam >ref.anchor.samal || true

		'''







rule gene_seq_dot_TE:
	input:
		rules.gene_seq_TEpos.output
	output:
		"{base_dir}/02_index_construction/TE_div/{chrom}/win.al.rm"
	shell:
		'''
		cd "{base_dir}02_index_construction/TE_div/{wildcards.chrom}"
		cat "{base_dir}01_preproc/cor_window.bed" |awk -v chr={wildcards.chrom} '$1==chr' >cor_window.bed
		cat > temp_script.R << 'EOF'
		library(data.table)
		library(parallel)
		data<-fread("reg.end.out")
		data1=fread("ref.anchor.samal")
		samlis=unique(data$V4)
		samlisdf=as.data.frame(samlis)
		fwrite(samlisdf,"sam.list",sep="\\t",col.names=FALSE)
		process_sample <- function(sam) {{
		data<-fread("reg.end.out")
		data1=fread("ref.anchor.samal")
		print(sam)
		samdf=data[data$V4==sam] ###SDR区域query的repeats组成
		samdf1=data1[data1$V4==sam] ###SDR区域坐标
		fwrite(samdf,paste(sam,".out",sep=""),sep="\\t",col.names=FALSE)
		fwrite(samdf1,paste(sam,".outpos",sep=""),sep="\\t",col.names=FALSE)
		common=paste("bedtools intersect -a cor_window.bed  -b ",sam,".out"," -wa -wb >",sam,".win.out",sep="")
		system(common)
		common=paste("bedtools intersect -a cor_window.bed  -b ",sam,".outpos"," -wa -wb >",sam,".win.outpos",sep="")
		system(common)
		common=paste("bedtools subtract -a cor_window.bed  -b ",sam,".outpos"," -A >",sam,".win.nointer.outpos",sep="")
		system(common)
		common=paste("python {script_dir}/index_construction/scripts/te_div/mergeSDR.py --sam ",sam,sep="")
		system(common)
		common=paste("bedtools intersect -a ",sam,".win.txt  -b ",sam,".out"," -wa -wb >",sam,".win.outn",sep="")
		system(common)
		common=paste("bedtools intersect -a ",sam,".win.txt  -b ",sam,".outpos"," -wa -wb >",sam,".win.outposn",sep="")
		system(common)
		return(paste("Completed:", sam))
		}}
		results <- mclapply(samlis, process_sample, mc.cores = 30)		
EOF
		Rscript temp_script.R
		bedtools intersect -a <(cat *win.txt cor_window.bed |sort |uniq) -b ref.rm -wa -wb |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >win.al.rm 

		'''




rule gene_seq_TE_final:
	input:
		rules.gene_seq_dot_TE.output
	output:
		"{base_dir}/02_index_construction/TE_div/{chrom}/TEEND"
	shell:
		'''
		cd "{base_dir}02_index_construction/TE_div/{wildcards.chrom}"
		python {script_dir}/index_construction/scripts/te_div/inte_allsam.py --chr {wildcards.chrom} --p 20
		python {script_dir}/index_construction/scripts/te_div/inte_sam2.py   --chr {wildcards.chrom}
		cat TEal.txt |tail -n +2 |cut -f 1-3 >TE.head
		Rscript {script_dir}/index_construction/scripts/te_div/filt_dot_cal.r	{base_dir}01_preproc/conti_path/allsampos	
		'''

rule inte_clear_TE:
	input:
		expand("{base_dir}/02_index_construction/TE_div/{chrom}/TEEND",base_dir=base_dir,chrom=allchrom)
	output:
		en="{base_dir}/02_index_construction/TE_div/TEEND",
		log="{base_dir}/02_index_construction/TE_div/TE_gene.OK",
	shell:
		'''
		cat {base_dir}/02_index_construction/TE_div/*/TEEND >{output.en}
		cd "{base_dir}/02_index_construction/TE_div/"
		find . -maxdepth 1 -type f ! -name 'TEEND' -delete
		find . -maxdepth 1 -type d ! -name '.' -exec rm -rf {{}} +
		touch {output.log}

		'''



rule gene_intraseg_initpos:
	input:
		f"{base_dir}01_preproc/gene_prep.OK",
		f"{base_dir}01_preproc/cor_window.bed"
		
	output:
		f"{base_dir}02_index_construction/intra_seg/intra_seg_extend.OK"
	shell:
		'''
		cd "{base_dir}02_index_construction/intra_seg/"
		bedtools getfasta -fi {REF} -bed {BED} -fo extracted_sequences.fa
		trf extracted_sequences.fa  2 7 7 80 10 50 500 -h -l 20 -ngs >trf.out
		python {script_dir}/index_construction/scripts/intra_seg/trf2bed.py --i trf.out --o trf.bed
		less -S trf.bed  |awk '$6>=7' |awk 'length($16)>100 && length($16)<10000 && $5>=3' >can.VNTR   ###candidate position
		less -S can.VNTR  |tr ":" "\\t" |tr "-" "\\t" |awk 'OFS="\\t"{{print $1,$2+$4,$2+$5,$7,$17}}' >can.VNTR.pos
		#### 感觉不是很好搞 先把对应的repeats在ref上的位置和motif搞出来
		python {script_dir}/index_construction/scripts/intra_seg/graph_VNTRend_onlypkl.py
		for i in {{0..99}};do
		echo $i
		python {script_dir}/index_construction/scripts/intra_seg/ref_extend.py --i ${{i}}pos.pkl --o ${{i}}pos.csv
		done
		touch intra_seg_extend.OK
		'''


rule gene_intraseg_ref_CN:
		input:
			rules.gene_intraseg_initpos.output
		output:
			f"{base_dir}02_index_construction/intra_seg/REF.CN"
		run:
			import glob   ###之后这里可以改成motif和CN具体的情况
			from collections import defaultdict
			import pandas as pd
			import os 
			import warnings
			warnings.filterwarnings('ignore')
			basedir=base_dir+"02_index_construction/intra_seg"
			print(input)
			print(basedir)
			allCON = defaultdict(dict)
			for file in glob.glob(f"{basedir}/*pos.csv"):
			    print(file)
			    try:
			        alldflisdf = pd.read_csv(file, sep="\t", header=0)
			        subdict=alldflisdf.groupby('0').apply(lambda x: dict(zip(x['2'], x['1']))).to_dict()
			        for key1, value1 in subdict.items():
			            allCON[key1].update(value1)
			    except pd.errors.EmptyDataError:
			        print("the file is empty")
			x = pd.DataFrame(allCON).T
			x=x[x[refname]>2] ##can be modify
			x.to_csv(basedir+"/REF.CN", index=True, sep="\t")


rule gene_intraseg_other_CN:
		input:
			rules.gene_intraseg_ref_CN.output
		output:
			f"{base_dir}02_index_construction/intra_seg/END.in"
		shell:
			'''
			cd "{base_dir}02_index_construction/intra_seg/"
			python {script_dir}/index_construction/scripts/intra_seg/graph_VNTR_other.py
			cat *.copynumber.csv |grep  "chr" |awk '$5!=1' >INSOTHER.al
			python {script_dir}/index_construction/scripts/intra_seg/process_otherpos1.py
			cat samlis |tail -n +2 | parallel -j 30 \
			  "python {script_dir}/index_construction/scripts/intra_seg/process_otherpos2.py --s {{}} --o {{}}.END"

			cat *END |grep "chr" |cut -f 5-7  >inter.in
			bedtools intersect -a inter.in  -b inter.in -wa -wb |awk '!($1==$4 && $2==$5 && $3==$6)' |less -S >x
			cat x |sort |uniq |awk '{{print $1":"$2"-"$3"\\t"$4":"$5"-"$6}}' >del_dup.in
			Rscript  {script_dir}/index_construction/scripts/intra_seg/del_dup.r
			cat *END |grep "chr"  >END.in

			'''


rule gene_intraseg_filter_otherCN:
		input:
			rules.gene_intraseg_other_CN.output
		output:
			f"{base_dir}02_index_construction/intra_seg/END.out"
		shell:
			'''
			cd "{base_dir}02_index_construction/intra_seg/"
			cat > temp_script.R << 'EOF'
			##R
			library(data.table)
			library(igraph)
			library(dplyr)
			data<-fread("END.in")
			delreg<-fread("del.reg",header=FALSE)
			data$name=paste(data$V5,":",data$V6,"-",data$V7,sep="")
			data$name2=paste(data$V10,"|",data$V11,"|",data$V12,"|",data$V13,sep="")
			delrepe=unique(data[data$name %in% delreg$V1,]$"name2") ##存在交集应该被删掉的repeats的情况 这样做不太对 直接通过repeats去删除会导致 比如 有的样本不应该被删掉 但是这里被删掉了 （逻辑应该改为去掉所有del.reg的样本的id后保留对应的repeats，这些repeats作为最终的所有的坐标）
			candata=data[!(data$name %in% delreg$V1),]
			saverepe=candata$name2
			data=data[(data$name2 %in% saverepe),]
			alis=c()
			for(group in unique(data$name)){{
			subdf=data[data$name==group,]
			can=subdf[1,]$name2
			for(x in subdf$name2){{alis<-rbind(alis,c(can,x))}}
			}}
			ou=as.data.frame(alis)
			colnames(ou) <- c("from", "to")
			g <- graph_from_data_frame(ou, directed = FALSE)
			comp <- components(g)
			subgraphs <- decompose.graph(g)
			delreg=c()
			for(i in 1:length(subgraphs)) {{
			  print(i)
			  subg <- subgraphs[[i]]
			  sub_df <- as.data.frame(vertex.attributes(subg))
			  delregap=unique(sub_df$name)[-1]
			  delreg<-c(delreg,delregap)
			  print(delreg)
			  }}
			unique_data <- data %>%
			  distinct(name, .keep_all = TRUE)
			fwrite(unique_data,"END.out",sep="\\t")  ###过滤了很多有重复的repeats等  还是有一些重复的 需要继续过滤
EOF
			Rscript temp_script.R
			'''

rule gene_intraseg_gini:
		input:
			rules.gene_intraseg_filter_otherCN.output
		output:
			f"{base_dir}02_index_construction/intra_seg/ENDVNTR.AL.cov"
		shell:
			'''
			cd "{base_dir}02_index_construction/intra_seg/"
			less -S trf.bed |cut -f 1-3,15 |less -S  |tr ":" "\\t" |tr "-" "\\t" |awk 'OFS="\\t"{{print $1,$4+$2,$5+$2,$6}}' >x
			python {script_dir}/index_construction/scripts/intra_seg/inte_intraseg.py
			Rscript {script_dir}/index_construction/scripts/intra_seg/gene_gini.r
			cp "{base_dir}01_preproc/cor_window.bed" ./
			bedtools intersect -b VNTRall.fill.VNTR.sd -a cor_window.bed -wa -wb >ENDVNTR.AL.inter
			bedtools coverage -b VNTRall.fill.VNTR.sd -a cor_window.bed >ENDVNTR.AL.cov
			'''

rule gene_intraseg_fea:
		input:
			rules.gene_intraseg_gini.output
		output:
			f"{base_dir}02_index_construction/intra_seg/intra_seg.fea"
		run:
			import pandas as pd
			sdpa=pd.read_csv(base_dir+"02_index_construction/intra_seg/ENDVNTR.AL.inter",sep="\t",header=None)
			end=sdpa.groupby([0, 1, 2])[6].agg(["max"]).reset_index()
			end1=sdpa.groupby([0, 1, 2])[7].agg(["max"]).reset_index()
			end=pd.concat([end,end1], axis=1)
			end=end.iloc[:,[0,1,2,3,7]]
			sdpa=pd.read_csv(base_dir+"02_index_construction/intra_seg/ENDVNTR.AL.cov",sep="\t",header=None)
			end1=sdpa.merge(end,on=[0,1,2],how="left")
			end1=end1.iloc[:,[0,1,2,6,7,8]]
			end1_filled = end1.fillna(0)
			end1_filled.columns=['chr','start','end','vntr_cov','vntr_sd_max','vntr_gini_max']
			end1_filled.to_csv(base_dir+"02_index_construction/intra_seg/intra_seg.fea", sep="\t", index=False, header=True)


rule clear_intra_fold: 
	input:
		rules.gene_intraseg_fea.output
	output:
		f"{base_dir}02_index_construction/intra_seg/intra_gene.OK"
	shell:
		"""
		cd "{base_dir}02_index_construction/intra_seg/"
		find . -maxdepth 1 -type f ! -name 'intra_seg.fea' -delete
		touch "intra_gene.OK"
		"""


rule find_al:
	input:
		rules.clear_intra_fold.output,
		rules.clear_bub_fold.output,
		rules.inte_clear_TE.output,
		rules.clear_interseg_fold.output,
		rules.clear_seqdiv_fold.output
	output:
		"{base_dir}/02_index_construction/gene_indices.OK"
	shell:
		'''
		touch {output}
		'''
##model

