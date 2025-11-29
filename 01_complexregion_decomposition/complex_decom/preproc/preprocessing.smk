### snakemake -s /home/jmhan/project/APG/github_final/1p_26_end/complex_decom/preproc/preprocessing.smk --keep-going --keep-incomplete -j 1 -R


# base_dir="/home/jmhan/project/APG/github_final/1p_26_end/run/"
# script_dir="/home/jmhan/project/APG/github_final/1p_26_end/complex_decom/"


# VCF="/home/jmhan/project/APG/github_final/1p_26_end/example/1p36.vcf"
# GFA="/home/jmhan/project/APG/github_final/1p_26_end/example/1p36.gfa"
# BED="/home/jmhan/project/APG/github_final/1p_26_end/example/1p36.bed" 
# sam_list="/home/jmhan/project/APG/github_final/1p_26_end/example/samlis.txt"

base_dir=config['base_dir']
script_dir=config['script_dir']
data_dir=config['data_dir']
VCF=data_dir+"1p36.vcf"
GFA=data_dir+"1p36.gfa"
BED=data_dir+"1p36.bed" 
sam_list=data_dir+"samlis.txt"


import pandas as pd
samdf = pd.read_csv(sam_list,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
name=list(samdf[0])
lines = name



# rule all:
# 	input:
# 		f"{base_dir}01_preproc/cor_window.bed",
# 		f"{base_dir}01_preproc/Wpath/gene_W.OK",
# 		f"{base_dir}01_preproc/reg.posCHM.txt",
# 		expand("{base_dir}01_preproc/conti_path/{li}_REFcor.pos",base_dir=base_dir,li=lines),
# 		f"{base_dir}01_preproc/conti_path/reg.vcf_corr_dot.en",
# 		f"{base_dir}01_preproc/conti_path/process_vcf.OK",
# 		f"{base_dir}01_preproc/gene_db.OK"


### generate windows
rule make_windows:
	input:
		{BED}
	output:
		f"{base_dir}01_preproc/cor_window.bed"
	shell:
		"""
		cd "{base_dir}"
		bedtools makewindows -b {BED} -w 1000 -s 500 |awk 'OFS="\\t"{{print $1,$2+1,$3}}'  >{output}
		echo "#####################SUCCESS_00-GENERATE THE WINDOWS#####################"
		"""

rule make_Wpath:
	input:
		{GFA}
	output:
		f"{base_dir}01_preproc/Wpath/gene_W.OK"
	shell:
		"""
		cd "{base_dir}01_preproc/Wpath/"
		cat {GFA} |awk '$1=="W"' >W.al
		awk '{{print > ($2 ".W")}}' W.al
		touch gene_W.OK
		echo "#####################SUCCESS_01-GENERATE THE PATHS######################"
		"""

rule make_json_pos:
	input:
		{GFA}
	output:
		f"{base_dir}01_preproc/reg.json",
		f"{base_dir}01_preproc/reg.posCHM.csv",
		f"{base_dir}01_preproc/reg.posCHM.txt"
	shell:
		'''
		cd "{base_dir}01_preproc/"
		REF=$(cat {sam_list} |grep "refe" |cut -f1)
		cat  {GFA}  |awk '$1=="S"' >reg.S
		cat  {GFA}  |awk '$1=="S"' |cut -f 2-3 >reg.S.in
		{{ echo "{{";  awk '{{if (NR > 1) print prev ","; prev = "\\"" $2 "\\"" ":" "\\"" $3 "\\""}} END {{print prev}}' reg.S;  echo "}}"; }} >reg.json
		less -S {GFA} |grep $REF |grep "^S" |cut -f2,8 |sed "s/,//g" |sed "s/#/\\t/g" |cut -f 1,4 |tr ":" "\\t" |tr "-" "\\t" |sort -k3,3n |less -S  >x
		less -S x | awk 'OFS="\\t"{{print $2,$3,$4,$1,">"}}' >reg.posCHM.txt
		rm x
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



### generate the continues region of each sample
rule gene_contireg: ##
	input:
		Wpath=rules.make_Wpath.output
	output:
		"{base_dir}01_preproc/conti_path/{li}_REFcor.pos"
	shell:
		"""
		cd "{base_dir}01_preproc/conti_path/"
		REF=$(cat {sam_list} |grep "refe" |cut -f1)
		python {script_dir}/preproc/scripts/align_W.py --sam {wildcards.li} --p "{base_dir}01_preproc/Wpath/" --r $REF
		python {script_dir}/preproc/scripts/align_W.pos.py --sam {wildcards.li} --p "{base_dir}01_preproc/Wpath/" --r $REF --reg "{base_dir}/01_preproc/reg.posCHM.txt"
		"""


	

### process the dot of VCF and the conflict regions
rule process_vcf: ##这里之后改为并行
	input:
		vcf={VCF},
		conseg=expand("{base_dir}01_preproc/conti_path/{li}_REFcor.pos",base_dir=base_dir, li=lines),
		pos=rules.make_json_pos.output
	output:
		f"{base_dir}01_preproc/conti_path/reg.vcf_corr_dot.en"
	shell:
		"""
		cd "{base_dir}01_preproc/conti_path/"
		cat *cor.pos  |awk 'OFS="\\t"{{print $2,$9,$10,$5}}' |grep -v "SAM" >allsampos

		cd {base_dir}01_preproc/conti_path/
		less -S {input.vcf}  |cut -f 1,2,10-1000 >tes
		less -S {input.vcf}  |grep -v "^#" >reg.vcf
		cat {input.vcf} | grep "#CHROM" >headern
		cat tes | grep "#CHROM" >header
		cat tes  |grep -v "^#" >test
		cat > temp_script.R << 'EOF'
		library(data.table)
		data<-fread("test",header=FALSE)
		colname<-fread("header")
		colnames(data)<-colnames(colname)
		colname=colnames(colname)[c(-1,-2)]
		for(sam in colname){{
		        print(sam)
		        samlis=c("#CHROM","POS",sam)
		        samdf=data[,..samlis]
		        fwrite(samdf,paste(sam,".samvcf",sep=""),sep="\\t")
		}}

		column_names <- data.frame(colname)
		fwrite(column_names, file = "allsamname", col.names = FALSE)
EOF

	Rscript temp_script.R
	##genevcfcore
	python {script_dir}/preproc/scripts/vcf_cor.py --vcf reg.vcf --o reg.vcfcor --seg "{base_dir}01_preproc/reg.posCHM.txt"
	while IFS= read -r sam; do
	    samfile=$sam
	    echo $samfile
	    cat {base_dir}01_preproc/conti_path/$samfile"_REFcor.pos" |cut -f2,9,10 |tail -n +2 >${{sam}}x.tem
	paste reg.vcfcor <(cat "$sam".samvcf |tail -n +2) | less -S | cut -f 1-3,6 >${{sam}}xx.tem
	bedtools coverage -a ${{sam}}xx.tem -b ${{sam}}x.tem >${{sam}}y.tem
	less -S ${{sam}}y.tem |awk '{{OFS="\\t"}}$8==1 && $4=="."{{$4=-1}}1'|cut -f 4| awk -v sam=$sam 'BEGIN{{print sam }}1'   > "$sam".samvcf.out
	done < allsamname
	paste *samvcf.out >reg.vcf_corr_dot

	less -S reg.vcf |nl -  |grep 'CONFLICT' |tr ";" "\\t" |less -S  |cut -f 1,13 >reg.conflict
	cat reg.conflict|tr "=" "\\n" |tr "," "\\n" |grep -v "CON" |sort |uniq >name
	while IFS= read -r sam; do
	cat reg.conflict |grep $sam  |cut -f1 >$sam.id
	awk 'NR==FNR {{lines[$1]; next}} FNR in lines {{$0="."}} 1' $sam.id  <(cat "$sam".samvcf.out |tail -n +2 ) | awk -v sam=$sam 'BEGIN{{print sam }}1'  >"$sam".samvcf.out.temp
	mv "$sam".samvcf.out.temp "$sam".samvcf.out
	done<name
	paste *samvcf.out >reg.vcf_corr_dot.en
	# rm -f !(reg.vcf_corr_dot.en|*REFcor.pos)
	"""

rule clear_fold: 
	input:
		rules.process_vcf.output
	output:
		f"{base_dir}01_preproc/conti_path/process_vcf.OK"
	shell:
		"""
		cd {base_dir}01_preproc/conti_path
		find . -maxdepth 1 -type f ! -name 'reg.vcf_corr_dot.en' ! -name '*REFcor.pos' ! -name 'allsampos' -delete
		touch process_vcf.OK
		"""


rule gene_seg_db:
	input:
		rules.make_Wpath.output
	output:
		f"{base_dir}01_preproc/gene_db.OK"
	run:
		import msgpack
		import lmdb
		import pandas as pd
		import shutil
		db_path = base_dir+'/01_preproc/seg_db'
		if os.path.exists(db_path):
		    shutil.rmtree(db_path)
		env = lmdb.open(db_path, map_size=100 * 1024 * 1024 * 1024)  # 请根据需求修改路径和大小
		nodedf = pd.read_csv(base_dir+"/01_preproc/reg.S.in",sep="\t",header=None)
		node_dict = dict(zip(nodedf[0], nodedf[1]))
		with env.begin(write=True) as txn:
		    for k, v in node_dict.items():
		        txn.put(str(k).encode('utf-8'), msgpack.packb(v))
		# with env.begin() as txn:
		#     [msgpack.unpackb(txn.get(key.encode('utf-8'))) for key in ['2','150','10']] #endfa['fa']=list(map(fasta_dict.get, letters_to_find))  #end.isnull().values.any()
		file_path = f"{base_dir}01_preproc/gene_db.OK"
		open(file_path, 'w').close()



rule prepr_al:
	input:
		rules.make_windows.output,
		rules.make_Wpath.output,
		rules.make_json_pos.output,
		rules.clear_fold.output,
		rules.gene_seg_db.output
	output:
		"{base_dir}/01_preproc/gene_prep.OK"
	shell:
		'''
		touch {output}
		'''