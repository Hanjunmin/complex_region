### ref:/home/jmhan/project/APG/github_final/PanSD/output/reference/filt/split/END/REF.end

# snakemake -s /home/jmhan/project/APG/github_final/PanSD/PanSD_anno1.smk -j 1

##构建数据库
import os
## cd /home/jmhan/iterall/
#snakemake -s /home/jmhan/iterall/interdb.smk --keep-going --keep-incomplete -j 10
import os
import glob
import msgpack
import lmdb
import pandas as pd
import shutil

script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"
base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"
DATA_dir="/home/jmhan/project/APG/github_final/PanSD_end/example/"
script=script_dir+"/ref_SD/sam_sd/scripts/",
T2Tsate="/home/jmhan/PANSDEND/DATA/CHM13.sat"
genomefa=DATA_dir+"FA/"
sam_list=DATA_dir+"samlis.txt"
rm="/home/jmhan/HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/asm_repeatmasker.out.bed"
RMsoft="/home/Public/software/RepeatMasker-4.1.4/"
import pandas as pd

import glob
import os
import pickle
samdf = pd.read_csv(sam_list,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
lines = samname
name=list(samdf[0])
FA=DATA_dir+"/FA/"+refname+".fa"

ALLSAM=name ###including reference

# 规则定义
# rule all:
#     input:
#     	f"{base_dir}03_refSD/sam_cor/Wpath/seg_db.OK",
#         expand("{base_dir}/03_refSD/sam_cor/Wpath/{sam}.W",base_dir=base_dir, sam=ALLSAM),
#         expand("{base_dir}/03_refSD/sam_cor/db/time/{sam}.txt",base_dir=base_dir, sam=ALLSAM),
#         f"{base_dir}03_refSD/sam_cor/1_SD_CANPATH/reproce.input",
#         f"{base_dir}03_refSD/sam_cor/1_SD_CANPATH/sd.edge",
#         expand("{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/time/{line}.txt",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/{line}.sd",base_dir=base_dir, line=lines),
#         f"{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/check_process_line.OK",
#         expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/other/{line}.other",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD/{line}.t2tsd",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/1.OK",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/2.OK",base_dir=base_dir, line=lines)



rule gene_Wpath:
    input:
        {sam_list}
    output:
        out2 = expand("{base_dir}/03_refSD/sam_cor/Wpath/{sam}.W", base_dir=base_dir, sam=ALLSAM),
        
    shell:
        '''
        mkdir -p {base_dir}/03_refSD/sam_cor/Wpath
        cd {base_dir}/03_refSD/sam_cor/Wpath
        cat {DATA_dir}/example_mc.gfa | awk '$1=="W" {{print > ($2 ".W")}}'
        cat {DATA_dir}/example_mc.gfa | awk '$1=="S"' | cut -f 2-3 > reg.S
        '''


rule gene_seg:
	input:
		rules.gene_Wpath.output
	output:
		f"{base_dir}03_refSD/sam_cor/Wpath/seg_db.OK"
	run:
		db_path = base_dir+'03_refSD/sam_cor/seg_db'
		if os.path.exists(db_path):
		    shutil.rmtree(db_path)
		env = lmdb.open(db_path, map_size=100 * 1024 * 1024 * 1024)  # 请根据需求修改路径和大小
		nodedf = pd.read_csv(base_dir+"03_refSD/sam_cor/Wpath/reg.S",sep="\t",header=None)
		node_dict = dict(zip(nodedf[0], nodedf[1]))
		with env.begin(write=True) as txn:
		    for k, v in node_dict.items():
		        txn.put(str(k).encode('utf-8'), msgpack.packb(v))
		# with env.begin() as txn:
		#     [msgpack.unpackb(txn.get(key.encode('utf-8'))) for key in ['2','150','10']] #endfa['fa']=list(map(fasta_dict.get, letters_to_find))  #end.isnull().values.any()
		file_path = f"{base_dir}03_refSD/sam_cor/Wpath/seg_db.OK"
		open(file_path, 'w').close()


rule gene_db:
    input:
        W="{base_dir}/03_refSD/sam_cor/Wpath/{sam}.W",
        S=rules.gene_seg.output,
    output:
        time_log="{base_dir}/03_refSD/sam_cor/db/time/{sam}.txt",  # 日志输出
        dirrec= directory("{base_dir}/03_refSD/sam_cor/db/{sam}")  # 数据库目录
    params:
        dbpath="{base_dir}/03_refSD/sam_cor/"  # 数据库路径
    shell:
        """
        mkdir -p {output.dirrec}  # 创建数据库目录
        mkdir -p {base_dir}/03_refSD/sam_cor/db/time
        /usr/bin/time -v python {script_dir}/ref_SD/sam_sd/scripts/inter_genedb.py --W {input.W} --S {wildcards.sam}  --p {params.dbpath}  2> {output.time_log}
        """


rule sta_gene_db:
    input:
        expand("{base_dir}/03_refSD/sam_cor/db/time/{sam}.txt",base_dir=base_dir, sam=ALLSAM),
    output:
        f"{base_dir}/03_refSD/sam_cor/db_check.OK"
    shell:
        '''
        touch {output}
        '''



rule gene_path:
    input:
        sd=f"{base_dir}03_refSD/ref_cor/filt/split/END/REF.end",
        can=f"{base_dir}01_intersim/candidate.merge"
    output:
        pa=f"{base_dir}03_refSD/sam_cor/1_SD_CANPATH/reproce.input",
        ed=f"{base_dir}03_refSD/sam_cor/1_SD_CANPATH/sd.edge"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/1_SD_CANPATH
		cd {base_dir}/03_refSD/sam_cor/1_SD_CANPATH
		cat {input.sd} |awk 'OFS="\\t"{{print $1":"$2"-"$3"\\t"$4":"$5"-"$6}}'  >sd.edge
        cat {input.sd} |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}' |less -S |sort |uniq |bedtools sort >endendcopy
        cat {input.can} |sort |uniq  |bedtools sort >endendcopy1
        # bedtools subtract -a endendcopy -b acro -A |grep -v "chrY" >endendcopy.filtacro
        # bedtools subtract -a endendcopy1 -b acro -A |grep -v "chrY" >endendcopy.filtacro1
        bedtools intersect -a endendcopy -b {base_dir}/00_preproc/reg.posCHM.txt -wa -wb|awk 'OFS="\\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' >sdinter.txt
        bedtools intersect -a endendcopy1 -b {base_dir}/00_preproc/reg.posCHM.txt -wa -wb|awk 'OFS="\\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' >sdinter1.txt
        Rscript {script_dir}/ref_SD/sam_sd/scripts/SDgenepath.r sdinter.txt alldf.txt
        Rscript {script_dir}/ref_SD/sam_sd/scripts/SDgenepath.r sdinter1.txt alldf1.txt
        cat alldf.txt |tail -n +2 |less -S  >alldfe.txt
        cat alldf1.txt |tail -n +2 |less -S  >alldfe1.txt
        cat alldfe1.txt alldfe.txt >{output.pa}
        """


rule process_line:
    input:
        W="{base_dir}/03_refSD/sam_cor/Wpath/{line}.W",
        inputallres=rules.gene_path.output.pa,
        infile=f"{base_dir}/03_refSD/sam_cor/db_check.OK"
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{line}can_sd.corregion.bed",  # 目标输出文件
        time_output="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/time/{line}.output.txt",  # 时间输出
        time_log="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/time/{line}.txt",  # 日志输出
    params:
        dbpath="{base_dir}/03_refSD/sam_cor/" 
    shell:
        """
        cd {base_dir}03_refSD/sam_cor/1_SD_CANPATH
        mkdir -p {base_dir}/03_refSD/sam_cor/1_SD_CANPATH/time
        /usr/bin/time -v python {script_dir}/ref_SD/sam_sd/scripts/iterpansd.py --W {input.W} --S {wildcards.line} --i {input.inputallres} --o {output.output_bed} --p {params.dbpath} > {output.time_output} 2> {output.time_log} ##这里面有个bedtools路径需要改
        """


rule check_process_line:
    input:
        expand("{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/time/{line}.txt",base_dir=base_dir, line=lines)
    output:
        f"{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/check_process_line.OK"
    shell:
        """
        touch {output}
        """


########NEXT STEP
rule init_sd_can:
    input:
        rules.process_line.output
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/{line}.sd"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        less -S {base_dir}/03_refSD/sam_cor/Wpath/{wildcards.line}.W | cut -f 2-6 | awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' > init.W
        cat "{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{wildcards.line}can_sd.corregion.bed" | sed \"s/_rows_/\\t/\" |tail -n +2 >{wildcards.line}.sd.pos
        Rscript {script}splidtcansd.r {base_dir}/03_refSD/sam_cor/1_SD_CANPATH/alldfe.txt {base_dir}/03_refSD/sam_cor/1_SD_CANPATH/alldfe1.txt {wildcards.line} {wildcards.line}.sd.pos init.W
        """


rule geneRMTRF:
    input:
        rules.init_sd_can.output,
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/RMTRFSAT.region"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        [ -d RMTRF ] && rm -r  RMTRF
        [ ! -d RMTRF ] && mkdir -p RMTRF
      cd RMTRF
      paste <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/{wildcards.line}.sd  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/{wildcards.line}.sd  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}') <(cat ../{wildcards.line}.sd|cut -f5)  |tail -n +2 >t2tsdrm.simn
      paste <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/{wildcards.line}.can  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/{wildcards.line}.can  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}') <(cat ../{wildcards.line}.can|cut -f5) |tail -n +2 >t2trm.simn
      cat t2tsdrm.simn  t2trm.simn |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7}}'>allrm
     awk '{{
    diff1 = $3 - $2; 
    diff2 = $6 - $5; 
    min = (diff1 < diff2) ? diff1 : diff2;  
    max = (diff1 > diff2) ? diff1 : diff2;  
    print $0"\\t"(max == 0 ? 0 : min / max);   
    }}' allrm >allrm.sim2

    less -S allrm.sim2 |awk '$7>=0.99 && $8>=0.99' >higncovn
    bedtools intersect -a {base_dir}/01_intersim/rm_trf.bed -b  <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -wa -wb >rm_inter_higncov
    bedtools subtract -a <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -b {base_dir}/01_intersim/rm_trf.bed -A  |cut -f 4-6|bedtools sort |bedtools merge >rm_inter_higncov.no
    Rscript /home/jmhan/PANSDEND/script/cal_RMTRF.r
    mv ALL.HIGHCOV.RMTRF ./RMTRF.END
    cat <(cat RMTRF.END |cut -f 1-3 |bedtools sort |bedtools merge) <(cat rm_inter_higncov.no) >RMTRF.region ###这个大样本区域
    cat RMTRF.END |cut -f 4-6 |bedtools sort |bedtools merge >RMTRF.rmtrf ###这个大样本区域对应的重复序列情况
    shopt -s nullglob
    bedtools intersect -a <(cat {T2Tsate}|cut -f 1-3) -b  <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -wa -wb >rm_inter_higncov
    bedtools subtract -a <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -b <(cat {T2Tsate}|cut -f 1-3) -A|cut -f 4-6 |bedtools sort |bedtools merge >rm_inter_higncov.no
    Rscript /home/jmhan/PANSDEND/script/cal_RMTRF.r
    mv ALL.HIGHCOV.RMTRF ./RMSAT.END
    cat <(cat RMSAT.END |cut -f 1-3 |bedtools sort |bedtools merge) <(cat rm_inter_higncov.no) >RMSAT.region ###这个大样本区域
    cat RMSAT.END |cut -f 4-6 |bedtools sort |bedtools merge >RMSAT.sat ###这个大样本区域对应的重复序列情况
    rm rm_inter_higncov.no
    rm allrm.sim2
    rm higncovn
    bedtools intersect -a RMSAT.region -b RMTRF.region >RMTRFSAT.region
    shopt -u nullglob
    """



rule process_line1_samSD:
    input:
        samnow=rules.process_line.output,
        trfrm="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/RMTRFSAT.region"
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/filin.bedend",
        updaten="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/rm_trf.s.bed",
        updaten1="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/T2Thave.bed",
        updaten2="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/filin.bedend"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        less -S {wildcards.line}.sd | tail -n +2 | awk 'OFS="\\t"{{print $7"#"$9,$10+$3,$10+$4}}' > allt2tsdcoorend.tsv
        less -S {wildcards.line}.sd |tail -n +2 |awk 'OFS="\\t"{{print $9,$10+$3,$10+$4}}' |less -S > allt2tsdcoorend_true.tsv
        less -S {wildcards.line}.can |tail -n +2 |awk 'OFS="\\t"{{print $7"#"$9,$10+$3,$10+$4}}' |less -S >candidateall.tsv
        less -S {wildcards.line}.can |tail -n +2 |awk 'OFS="\\t"{{print $9,$10+$3,$10+$4}}' |less -S >candidateall_true.tsv
        Rscript {script}T2Tmatch.r {base_dir}/03_refSD/sam_cor/1_SD_CANPATH/sd.edge {wildcards.line}.sd
        cat <(less -S process.txt |grep -v "nosam" |tail -n +2 |cut -f1,2 |tr "\\t" "\\n" |less -S |sort |uniq) <(cat end.txt |tail -n +2|cut -f1,2 |tr "\\t" "\\n") |less -S >reprocess.pair
        cat <(less -S process.txt|grep -v "nosam" |tail -n +2 |cut -f1,2) <(cat end.txt |tail -n +2 |cut -f1,2) |less -S >reprocess.pairall
        grep -Ff <(cut -f1 reprocess.pair) {wildcards.line}.sd > {wildcards.line}.multicoordi.inte.tsv
        less -S {wildcards.line}.multicoordi.inte.tsv |awk 'OFS="\\t"{{print $2,$7"@"$9"@"$10"@"$11"@"$3"@"$4}}' |less -S  >cor.id
        Rscript {script}intepair.r
        less -S pairall |cut -f 3,4 |less -S |tr "@" "\\t" |less -S |awk 'OFS="\\t"{{print $2,$3+$5,$3+$6,$8,$9+$11,$9+$12}}' |less -S |tail -n +2 >pairidall
        paste pairidall <(cat pairall | tail -n +2|cut -f 2) <(cat pairall|tail -n +2 |cut -f 1) >T2T2SAM
        less -S T2T2SAM  |awk 'OFS="\\t"{{print $1":"$2"-"$3,$7,$4":"$5"-"$6,$8}}' |awk 'OFS="\\t"{{print $1,$2"\\n"$3,$4}}' |less -S  |sort |uniq>T2T2SAM.pair
        echo "X"
        less -S process.txt |grep "nosam" >processnosam.txt ||true ##(这是第二种情况，是找不到对应的区域，相当于是变化极大)
        echo "X1"
        [ -d t2tsdminimap ] && rm -r  t2tsdminimap
        [ ! -d t2tsdminimap ] && mkdir -p t2tsdminimap
        cd t2tsdminimap
        fa=$(find {genomefa}{wildcards.line}*.fa)
        cat $fa >{wildcards.line}.fa
        samtools faidx {wildcards.line}.fa
        echo "A"

      [ -f all.paf ] && rm all.paf
      while IFS=$'\\t' read -r col1 col2 col3 col4 col5 col6
      do
      minimap2 -c --eqx <(samtools faidx  $fa "$col1":"$col2"-"$col3") <(samtools faidx  $fa "$col4":"$col5"-"$col6" ) -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >>all.paf
      done < ../pairidall
      python {script}paf2bed.py  --paf  all.paf  --o  all.end.statistics
      less -S all.end.statistics  |awk '$10>=0.9 && $11>=1000 && $12>=0.5' >all.end.statistics.filt
      less -S all.end.statistics.filt  |awk '$15>=0.8 && $16>=0.8' >T2Thave.bed
      less -S all.end.statistics.filt | awk '!($15>=0.8 && $16>=0.8)' > T2Tpolymor.bed 
      less -S T2Tpolymor.bed | awk '$15>=0.95 || $16>=0.95' >T2Tpolymor1.bed ##(√) ###这部分就不计算trf和repeatmasker了
      less -S T2Tpolymor.bed | awk '!($15>=0.95 || $16>=0.95)' > T2Tpolymor1other.bed
      cut -f 14 T2Tpolymor1other.bed >all_values.txt
      Rscript -e 'library(data.table)
        data<-fread("all_values.txt")
        data2<-fread("../T2T2SAM.pair",header=FALSE)
        colnames(data2)<-c("quesou","T2T")
        df <- merge(data, data2, by = "quesou", all.x = TRUE)
        fwrite(df,"other_T2T",sep="\\t")'
    echo "B"
    if [ -s "other_T2T" ] && [ $(wc -l < "other_T2T") -gt 1 ]; then
        less other_T2T |awk 'OFS="\\t"{{print $2,$1}}' | sort | uniq | tr ":" "\\t" | awk 'BEGIN {{OFS="\\t"}} {{gsub("-", "\\t", $2); gsub("-", "\\t", $4); print}}' | less -S |grep -v "T2T" |less -S >t2trm
        bedtools intersect -a t2trm -b <(cat {rm} |cut -f 1-3) -wa -wb >groupin
        bedtools subtract -a t2trm -b <(cat {rm} |cut -f 1-3) -A |cut -f 4-6>groupin.ana
        Rscript {script}expand_rmtrf.r "groupin" "{wildcards.line}.fa.fai" "rm_in" 2 ##这里的2表示左右扩的时候多扩两倍
        cat "rm_in" groupin.ana >rm_in.end
        bedtools coverage -a rm_in.end -b ../RMTRF/RMTRFSAT.region |awk '$7!=1' >rmin.in
        bash {script}run_rm_trf.sh rmin.in {wildcards.line} {RMsoft} {script} ##这步为跑trf和rm  
    else
        touch rm.out
        touch rm_trf.s.bed
        touch rmin.in
    fi
    echo "C"
       
      cat <(cat ../RMTRF/RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat rm.out |cut -f 1-3,5 )  >allsat
      cat rm_trf.s.bed ../RMTRF/RMTRF.rmtrf  >allrmtrf
      bash {script}transe.sh allrmtrf allrmtrf.trane
      bash {script}endfilt.sh allsat allrmtrf.trane T2Tpolymor1other.bed filin.bedend
      cd ..
      cd ..
      cp ./{wildcards.line}/t2tsdminimap/filin.bedend {output.output_bed}
       
       """


#####很重要的三个文件！！！RMTRFSAT.region RMSAT.sat RMTRF.rmtrf 
rule process_line2_samSD:
    input:
        rules.process_line1_samSD.output,
        samsim=rules.process_line.output
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/can/filin.bedend",
        updaten="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/can/rm_trf.s.bed",
        updaten1="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/can/groupin.ana",
        updaten2="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/can/rm_in"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        [ -d can ] && rm -r  can
        [ ! -d can ] && mkdir -p can
         
      cd can
      paste <(cat ../{wildcards.line}.sd  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.line}.sd  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}')  |less -S |tail -n +2 >t2tsdrm
      bedtools intersect -a t2tsdrm -b <(cat {base_dir}/01_intersim/rm_trf.bed |cut -f 1-3) -wa -wb >groupinsd
      echo "A"
      paste <(cat ../{wildcards.line}.can  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.line}.can  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}')  |less -S |tail -n +2 >t2trm
      bedtools intersect -a t2trm -b <(cat {base_dir}/01_intersim/rm_trf.bed |cut -f 1-3) -wa -wb >groupin
      bedtools subtract -a t2trm -b <(cat  {base_dir}/01_intersim/rm_trf.bed |cut -f 1-3 |bedtools sort |bedtools merge) -A  |cut -f 4-6>groupin.ana
      Rscript {script}expand_rmtrf.r "groupin" "../t2tsdminimap/{wildcards.line}.fa.fai" "rm_in" 1
      bash {script}transe.sh rm_in rm_in.transe
      mv rm_in.transe rm_in
      cat ../t2tsdminimap/T2Thave.bed   ../t2tsdminimap/T2Tpolymor1.bed ../t2tsdminimap/filin.bedend |cut -f 1-12 >t2tcoresbon.end ##*****************
      bedtools subtract -a <(cat rm_in groupin.ana) -b <(cat t2tcoresbon.end |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}') >candidat.sub
      bedtools subtract -a <(cat  candidat.sub |bedtools sort |bedtools merge) -b  ../allt2tsdcoorend_true.tsv >candidat.sub2
      [ -f allSDsample.fa ] && rm allSDsample.fa
      fa="../t2tsdminimap/{wildcards.line}.fa"
      echo "B"
      while IFS= read -r lin
      do
          chr=$(echo "$lin" | cut -f1)
          start=$(echo "$lin" | cut  -f2)
          end=$(echo "$lin" | cut  -f3)
          samtools faidx $fa "$chr:$start-$end" >> allSDsample.fa
      done < ../allt2tsdcoorend_true.tsv
      samtools faidx $fa $(cat candidat.sub2|awk '{{print $1":"$2"-"$3}}' |tr "\\n" "\\t" )  >candidat.sub2.fa
      minimap2 -c --eqx candidat.sub2.fa allSDsample.fa -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 > t2tsdsec2.paf
      python {script}paf2bed.py  --paf  t2tsdsec2.paf  --o  all.end.statistics
      less -S all.end.statistics  |awk '$10>=0.9 && $11>=1000 && $12>=0.5' >all.end.statistics.filt
    less -S all.end.statistics.filt |tail -n +2 |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}' >rminref
    bedtools coverage -a rminref -b ../RMTRF/RMTRFSAT.region |awk '$7!=1' >rmin.in
    cp ../t2tsdminimap/{wildcards.line}.fa ./
     echo "C"

    if [ -s "rmin.in" ]; then
        bash {script}run_rm_trf.sh rmin.in {wildcards.line} {RMsoft} {script} ##这步为跑trf和rm
    else
        touch rm.out
        touch rm_trf.s.bed
    fi

    echo "D"
    cat <(cat ../RMTRF/RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat rm.out |cut -f 1-3,5 )  >allsat
    cat rm_trf.s.bed ../RMTRF/RMTRF.rmtrf  >allrmtrf
    bash {script}transe.sh allrmtrf allrmtrf.trane
    bash {script}endfilt.sh allsat allrmtrf.trane all.end.statistics.filt filin.bedend

"""

rule merge_end_samSD:
    input:
        rules.process_line2_samSD.output,
    output:
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/SDend"
    threads: 1
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}
        mkdir -p endfilt && cd endfilt
        rm -f *
        if [ -s ../can/filin.bedend ]; then
            Rscript {script}ednfilt.r ../can/filin.bedend  {base_dir}03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/{wildcards.line}.fa {script}"paf2bed.py"
            ls |grep "csv"|less -S |tr "_" "\\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
            Rscript {script}run2filt.r
            touch nomatch.csv
            touch beforeminus0.9.csv
            cat *beforeminus0.9.csv *nomatch.csv *after.csv *nointe.csv  | tr "," "\\t" | tail -n  +2 |awk -F'\\t' '$10>=0.9 && $11>=1000 && $12>=0.5{{print $0}}'|less -S |grep -v "V1" |sort |uniq >ourend
            cd ..
            cp ./endfilt/ourend ./{output}



            mkdir -p endfilt1 && cd endfilt1
            rm -f *
            less -S  ourend |cut -f 1-12 >panSD_T2T.bed
            Rscript {script}endfilt2.r panSD_T2T.bed  {base_dir}03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/{wildcards.line}.fa {script}"paf2bed.py"
            ls |grep "csv"|less -S |tr "_" "\\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
            Rscript {script}run2filt.r
            touch nomatch.csv
            touch beforeminus0.9.csv
            cat *beforeminus0.9.csv *nomatch.csv *after.csv *nointe.csv  | tr "," "\\t" | tail -n  +2 |awk -F'\\t' '$10>=0.9 && $11>=1000 && $12>=0.5{{print $0}}'|less -S |grep -v "V1"|cut -f 1-12 |sort |uniq >ourend
            cd ..
            cp ./endfilt1/ourend ./{output}
        else
            cd ..
            touch SDend
        fi
        """

rule update_can_reg:
    input:
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/SDend",base_dir=base_dir,line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/can/rm_in",base_dir=base_dir,line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/can/groupin.ana",base_dir=base_dir,line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/T2Thave.bed",base_dir=base_dir,line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/t2tsdminimap/filin.bedend",base_dir=base_dir,line=lines)
    output:
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/other/{line}.other",
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD/{line}.t2tsd"
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/other
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD
        bedtools subtract -a <(cat {wildcards.line}/can/rm_in {wildcards.line}/can/groupin.ana) -b <(cat {wildcards.line}/t2tsdminimap/T2Thave.bed   {wildcards.line}/t2tsdminimap/T2Tpolymor1.bed {wildcards.line}/t2tsdminimap/filin.bedend {wildcards.line}/SDend |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}') |bedtools sort |bedtools merge >./other/{wildcards.line}".other"
        cat {wildcards.line}/t2tsdminimap/T2Thave.bed   {wildcards.line}/t2tsdminimap/T2Tpolymor1.bed {wildcards.line}/t2tsdminimap/filin.bedend {wildcards.line}/SDend >./T2TPANSD/{wildcards.line}".t2tsd"
        """





rule update_rmtrf_ref:
    input:
        rmnow=rules.process_line1_samSD.output
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/1.OK"
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/RMTRF
        cat <(cat {base_dir}03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/rmin.in |cut -f 1-3)  RMTRFSAT.region  >RMTRFSAT.region.temp
        cat {base_dir}03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/rm_trf.s.bed  RMTRF.rmtrf >RMTRF.rmtrf.temp
        cat <(cat RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/rm.out |cut -f 1-3,5 )  >RMSAT.sat.temp
        mv RMTRFSAT.region.temp  RMTRFSAT.region
        mv RMTRF.rmtrf.temp  RMTRF.rmtrf
        mv RMSAT.sat.temp  RMSAT.sat
        touch 1.OK
        """


rule update_rmtrf_can:
    input:
        rmnow=rules.process_line2_samSD.output
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/2.OK"
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/RMTRF
        cat <(cat {base_dir}03_refSD/sam_cor/2_ALIGN/{wildcards.line}/can/rmin.in |cut -f 1-3)  RMTRFSAT.region  >RMTRFSAT.region.temp
        cat {base_dir}03_refSD/sam_cor/2_ALIGN/{wildcards.line}/can/rm_trf.s.bed  RMTRF.rmtrf >RMTRF.rmtrf.temp
        cat <(cat RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/can/rm.out |cut -f 1-3,5 )  >RMSAT.sat.temp
        mv RMTRFSAT.region.temp  RMTRFSAT.region
        mv RMTRF.rmtrf.temp  RMTRF.rmtrf
        mv RMSAT.sat.temp  RMSAT.sat
        touch 2.OK
        """