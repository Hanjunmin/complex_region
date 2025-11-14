### ref:/home/jmhan/project/APG/github_final/PanSD/output/reference/filt/split/END/REF.end

# snakemake -s /home/jmhan/project/APG/github_final/PanSD/script/nonreference_SD/iter_END.smk -j 1

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
FA=DATA_dir+"/FA/"+refname+".fa"

ALLSAM=name ###including reference


# base_dir="/home/jmhan/project/APG/github_final/PanSD/output/"
# script_dir="/home/jmhan/project/APG/github_final/PanSD/script/"
SAMRMTRF=base_dir+"/03_refSD/sam_cor/2_ALIGN/" ###每个样本对应的RM和TRF的区域

# rule all:
#     input:
#         expand("{base_dir}03_nonref/02_ITER/2_SD_GENE/{line}/SPLIT.OK",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/expand.out",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapsplit/part_{i}", i=["{:02d}".format(i) for i in range(1, 51)],line=lines,base_dir=base_dir),
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/part_{i}.paf", i=["{:02d}".format(i) for i in range(1, 51)],line=lines,base_dir=base_dir),
#         f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/gene_paf.OK",
#         f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/gene_expand.OK",
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/all1.end.statistics",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/samuniq.pair",base_dir=base_dir, line=lines),
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/otherdfaddsam.txt",base_dir=base_dir, line=lines),
#         f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/df.pair",
#         f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/allpair.sort.END.e",
#         expand("{base_dir}03_nonref/02_ITER/3_PAIR_PATH/output/{sam}can_sd.corregion.bed",base_dir=base_dir, sam=ALLSAM),
#         expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/{line}.SD.trans",base_dir=base_dir, line=lines)


rule gene_samsim:
    input:
        f"{base_dir}01_intersim/sim.out"
    output:
        corespos=f"{base_dir}03_nonref/01_simpath/CHM13_APG_HPRC_HGSVC.pos_iter.txt",
        out=expand("{base_dir}03_nonref/01_simpath/{sam}.sim.sd.bed",base_dir=base_dir, sam=ALLSAM)
    shell:
        '''
        mkdir -p {base_dir}03_nonref/01_simpath/
        cd {base_dir}03_nonref/01_simpath/
        cat {input} |cut -f 2- >sim.temp
		IFS=$'\\t' read -r -a col_names < <(head -n 1 sim.temp)
        cols=${{#col_names[@]}}
        for ((j=1; j<=cols; j++)); do
            cut -f$j sim.temp > "${{col_names[j-1]}}.txt"
        done 
        for ((j=1; j<=cols; j++)); do
            paste win.txt "${{col_names[j-1]}}.txt" |tail -n +2 |tr ":" "\\t" |tr "-" "\\t"  |grep -v "e" |grep -v "win"  |awk '$4>=0.05'  |bedtools sort |bedtools merge > "${{col_names[j-1]}}.sim.sd.bed" 
        done 
        cat T2TCHM13.sim.sd.bed >CHM13v2.sim.sd.bed ###modify为对应的ref的样本的名字
        cat *sim.sd.bed  |bedtools sort |bedtools merge >cores.pos
        bedtools intersect -a cores.pos -b {base_dir}/00_preproc/reg.posCHM.txt -wa -wb |cut -f 4-8 >CHM13_APG_HPRC_HGSVC.pos_iter.txt
        '''

######ADD ACRO FILTER RULE

rule gene_sam_path:
    input:
        expand("{base_dir}03_nonref/01_simpath/{sam}.sim.sd.bed",base_dir=base_dir, sam=ALLSAM),
        rules.gene_samsim.output
    output:
        "{base_dir}03_nonref/01_simpath/{sam}alldf.e.txt"
    shell:
        '''
        cd {base_dir}03_nonref/01_simpath/
        bedtools intersect -a {wildcards.sam}.sim.sd.bed -b CHM13_APG_HPRC_HGSVC.pos_iter.txt -wa -wb | awk 'OFS="\\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' > {wildcards.sam}siminter.txt
        Rscript {script_dir}/nonref_SD/scripts/SDgenepath.r {wildcards.sam}siminter.txt {wildcards.sam}alldf.e.txt
        rm {wildcards.sam}siminter.txt
        '''

rule check_sam_path:
    input:
        expand("{base_dir}03_nonref/01_simpath/{sam}alldf.e.txt",base_dir=base_dir,sam=ALLSAM)
    output:
        f"{base_dir}03_nonref/01_simpath/simpath.OK"
    shell:
        '''
        touch {output}
        '''




rule gene_simapa:
    input:
        W="{base_dir}/03_refSD/sam_cor/Wpath/{line}.W",  # WCOM=rules.gene_Wpath.output,
        pain=rules.check_sam_path.output # inputallres=rules.gene_path.output.pa
    output:
        output_bed="{base_dir}03_nonref/02_ITER/1_SD_CANPATH/output/{line}can_sd.corregion.bed",  # 目标输出文件
        time_output="{base_dir}03_nonref/02_ITER/1_SD_CANPATH/time/{line}.output.txt",  # 时间输出
        time_log="{base_dir}03_nonref/02_ITER/1_SD_CANPATH/time/{line}.txt",  # 日志输出
    params:
        dbpath="{base_dir}/03_refSD/sam_cor/"  # 数据库路径
    shell:
        '''
        cd {base_dir}03_nonref/02_ITER/1_SD_CANPATH
        /usr/bin/time -v python {script_dir}/ref_SD/sam_sd/scripts/iterpansd.py  --W {input.W} --S {wildcards.line} --i "{base_dir}03_nonref/01_simpath/{wildcards.line}alldf.e.txt" --o {output.output_bed} --p {params.dbpath} > {output.time_output} 2> {output.time_log} ##这里面有个bedtools路径需要改 iterpansd.py 也可以改为两步骤 iterpansd1.py iterpansd2.py  
        '''

###################### SD calculation
rule process_line1:
    input:
        samwin="{base_dir}/03_refSD/sam_cor/Wpath/{line}.W",
        sammiin="{base_dir}/03_refSD/sam_cor/2_ALIGN/other/{line}.other",
        pa=rules.gene_simapa.output
    output:
        output_bed="{base_dir}03_nonref/02_ITER/2_SD_GENE/{line}/newin"
    shell:
        '''
        mkdir -p {base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.line}
        cd {base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.line}
        less -S {input.samwin} | cut -f 2-6 | awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' > init.W
        cat "{base_dir}03_nonref/02_ITER/1_SD_CANPATH/output/{wildcards.line}can_sd.corregion.bed" |sed \"s/_rows_/\\t/\" | awk '{{OFS="\\t"}}NR>1 {{print $1,$2,$3,$4,$6}} NR==1 {{print $0}}' >{wildcards.line}.sd.pos
        Rscript {script_dir}/nonref_SD/scripts/winsimprocess.r {wildcards.line}.sd.pos "{base_dir}03_nonref/02_ITER/1_SD_CANPATH/output/{wildcards.line}can_sd.corregion.bed.add"
        less -S "{base_dir}03_nonref/02_ITER/1_SD_CANPATH/output/{wildcards.line}can_sd.corregion.bed.add"  |awk 'OFS="\\t"{{print $6"#"$8,$9+$3,$9+$4}}' |less -S > wincoorend.tsv
        less -S "{base_dir}03_nonref/02_ITER/1_SD_CANPATH/output/{wildcards.line}can_sd.corregion.bed.add"  |awk 'OFS="\\t"{{print $8,$9+$3,$9+$4}}' |less -S > winend_true.tsv
        bedtools intersect -a {input.sammiin} -b winend_true.tsv |less -S |bedtools sort |bedtools merge >winbef.bed
        bedtools makewindows -b winbef.bed -w 1000 -s 500 > {output.output_bed}
        '''




rule prallel_kmercal:
    input:
        win_bed=rules.process_line1.output.output_bed
    output:
        output_bed="{base_dir}03_nonref/02_ITER/2_SD_GENE/{line}/filt/beforegroup.txt", #    threads: 7

    shell:
        """
        cd {base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.line}
        rm -r filt || true
        rm *.txt || true
        rm *.pkl || true
        samtools faidx {DATA_dir}/FA/{wildcards.line}.fa
        python  {script_dir}/nonref_SD/scripts/parallel_kmer.py --fa {DATA_dir}/FA/{wildcards.line}.fa --i newin --n 20
        python {script_dir}/ref_SD/SD_cordi/scripts/writxt.py
        rm *.pkl || true
        find . -maxdepth 1 -name "*.txt" -type f -size 0 -exec rm {{}} \;
        python {script_dir}/nonref_SD/scripts/para_kmer_wri.py
        """





rule cluster_expand:
    input:
        bef=rules.prallel_kmercal.output
    output:
        outbed="{base_dir}03_nonref/02_ITER/2_SD_GENE/{line}/SPLIT.OK"
    shell:
        """
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/filt
        cp {input.bef} ./beforegroup.init
        less -S beforegroup.txt  |cut -f 1-3 |less -s |uniq -c |less -S |awk '$1>100' >x
        less -S x |awk 'OFS="\\t"{{print $2,$3,$4}}' >rmin
        bedtools coverage -a rmin -b {SAMRMTRF}/{wildcards.line}/RMTRF/RMTRFSAT.region |awk '$7!=1' >rmin.in
        if [ -s "rmin.in" ]; then
        bash {script_dir}/ref_SD/sam_sd/scripts/run_rm_trf.sh rmin.in {DATA_dir}/FA/{wildcards.line}.fa  ##这步为跑trf和rm 只是那些低覆盖区域的
        else
        touch rm_trf.s.bed
        fi
        bedtools coverage -a rmin -b <(cat rm_trf.s.bed {SAMRMTRF}/{wildcards.line}/RMTRF/RMTRF.rmtrf)  |awk '$7>=0.9'|cut -f 1-3 >del.pair
        if [ -s del.pair ]; then
          awk 'NR==FNR {{a[$1,$2,$3]; next}} !($1,$2,$3) in a' del.pair beforegroup.txt >beforegroup.filt.txt
          cp beforegroup.filt.txt beforegroup.txt
        fi
        
        cp beforegroup.txt ../
        Rscript {script_dir}/nonref_SD/scripts/splitgroup.r {base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/filt/split {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/
        cd split
        python {script_dir}/nonref_SD/scripts/gene_chr_para.py -b {base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/filt/split/ -f {DATA_dir}/FA/{wildcards.line}.fa -t 40 
        # snakemake -s /home/jmhan/pansdall/SD1.smk  --config base_dir={base_dir}03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/filt/split/ fa={DATA_dir}/FA/{wildcards.line}.fa -j 40 ##记得路径加/
        touch  {output.outbed}
        """



rule expandt2trmtrf:
    input:
        ha=rules.cluster_expand.output,
    output:
        outbed="{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/expand.out",
    shell:
        """
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/filt/split
        paadd="{base_dir}03_nonref/02_ITER/1_SD_CANPATH/output/{wildcards.line}can_sd.corregion.bed.add"
        paste <( cat $paadd  |cut -f2 |tr ":" "\\t" |tr "-" "\\t" ) <(cat $paadd |awk 'OFS="\\t"{{print $8,$9+$3,$9+$4}}') >winsim_t2t.bed
        bash {script_dir}/nonref_SD/scripts/expandall.sh "{base_dir}01_intersim/rm_trf_cut1M.bed"  {DATA_dir}/FA/{wildcards.line}.fa.fai  /home/jmhan/pangenome/WGSfea/sdseq/expand_pansd.r ###bash /home/jmhan/pansdall/expandall.sh /home/jmhan/pangenome/WGSfea/allW/APG3allWnew/scripts/t2trmtrf.cut1M.bed  /home/jmhan/allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/C001-CHA-E01-Mat.fa.fai  /home/jmhan/pangenome/WGSfea/sdseq/expand_pansd.r
        rm rigin|| true
        rm lefin|| true
        rm *pair|| true
        rm *end.txt|| true
        cp {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/filt/split/expand.out {output.outbed}
        """

rule expandt2trmtrf1:
    input:
        f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/{{line}}/expand.out",
    output:
        parts = expand(f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/{{{{line}}}}/minimapsplit/part_{{i}}",i=["{:02d}".format(i) for i in range(1, 51)])
    shell:
        """
        mkdir -p {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/minimapsplit
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/minimapsplit
        split -n l/50 {input} --numeric-suffixes=1 --suffix-length=0 part_
        """


rule gene_expan_end:
    input:
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/expand.out",line=lines,base_dir=base_dir)
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/gene_expand.OK"
    shell:
        '''
        touch {output}
        '''


rule minimap_reg:
    input:
        last=rules.gene_expan_end.output,
        allmini="{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapsplit/part_{i}"
    output:
        paf="{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/part_{i}.paf"
    shell:
        '''
        cd "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/minimapsplit"
        if [ -s {input.allmini} ]; then
        prefix=$(basename {input.allmini})
        minimap_dir="{base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/minimapspnow"
        less -S {input.allmini} |awk '$1==$4 && $2<=$6 && $3>=$5' >$minimap_dir/$prefix.inter ||true
        if [ -s $minimap_dir/$prefix.inter ];then
        cat $minimap_dir/$prefix.inter |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' > $minimap_dir/$prefix.inter.bed
        less -S $minimap_dir/$prefix.inter.bed |cut -f 1-3 |sort |uniq >$minimap_dir/$prefix.left
        while IFS=$'\\t' read -r col1 col2 col3
        do
        awk -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' $minimap_dir/$prefix.inter.bed > $minimap_dir/$prefix.right
        minimap2 -c --eqx <(samtools faidx {DATA_dir}/FA/{wildcards.line}.fa  "$col1:$col2-$col3")  <(samtools faidx {DATA_dir}/FA/{wildcards.line}.fa  $(cat $minimap_dir/$prefix.right |awk '{{print $4":"$5"-"$6}}'))  -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -r 1000,1000 -P >> {output.paf}
        done < $minimap_dir/$prefix.left
        fi
        



        cat {input.allmini} |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' > $minimap_dir/$prefix.bed
        
        less -S $minimap_dir/$prefix.bed |cut -f 1-3 |sort |uniq >$minimap_dir/$prefix.left
        while IFS=$'\\t' read -r col1 col2 col3
        do
            awk -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' $minimap_dir/$prefix.bed > $minimap_dir/$prefix.right
            minimap2 -c --eqx <(samtools faidx {DATA_dir}/FA/{wildcards.line}.fa $(cat $minimap_dir/$prefix.right |awk '{{print $4":"$5"-"$6}}')) <(samtools faidx {DATA_dir}/FA/{wildcards.line}.fa "$col1:$col2-$col3")  -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >> {output.paf}
        done < $minimap_dir/$prefix.left
        else
        touch {output.paf}
        fi
        '''



rule mini_comple:
    input:
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/part_{i}.paf", i=["{:02d}".format(i) for i in range(1, 51)],line=lines,base_dir=base_dir)
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/gene_paf.OK"
    shell:
        """
        touch {output}
        """




rule combine_minimap:
    input:
         rules.mini_comple.output
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/all1.end.statistics"
    threads: 1  # 限制线程数为1
    shell:
        """
        minimap_dir="{base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/minimapspnow"
        find $minimap_dir -name "*.paf" -exec cat {{}} + > $minimap_dir/repea_mini
        python {script_dir}/nonref_SD/scripts/paf2bed.py  --paf  $minimap_dir/repea_mini  --o  $minimap_dir/all1.end.statistics
        """


rule filter_process: ###算啦好像本来生成的新的就很少，就先这样不cluster了
    input:
        rules.combine_minimap.output
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/samuniq.pair",
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/other.end"
    threads: 1  
    shell:
        """
        mkdir -p {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/processend ##modify
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/processend
        cat {input} | awk '$10>=0.85 && $11>=850 && $12>=0.5' >all1.end.statistics.filt
        less -S all1.end.statistics.filt  |tail -n +2 |cut -f 1-3  |bedtools sort |bedtools merge >s1.bed
        awk 'OFS="\\t"{{if ($2 == 0) $2 = 1}}1' s1.bed > s1.bedo
        less -S all1.end.statistics.filt  |tail -n +2 |cut -f 4-6  |bedtools sort |bedtools merge >s2.bed
        awk 'OFS="\\t"{{if ($2 == 0) $2 = 1}}1' s2.bed > s2.bedo
        s1=$(cat s1.bed|awk '{{sum+=($3-$2)}}END{{print sum}}')
        s2=$(cat s2.bed|awk '{{sum+=($3-$2)}}END{{print sum}}')
        samfa={DATA_dir}/FA/{wildcards.line}.fa
        if (( $(echo "$s1 < $s2" | bc -l) )); then ##前三列长度更短
            echo "1"
            bash {script_dir}nonref_SD/scripts/run_rm_trf.sh s1.bedo $samfa
            cat all1.end.statistics.filt |cut -f 1-12 |tail -n +2  >all1.end.statistics.filt.in
        else  ##后三列长度更短
            echo "2"
            bash {script_dir}nonref_SD/scripts/run_rm_trf.sh s2.bedo $samfa
            cat all1.end.statistics.filt |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' |tail -n +2  >all1.end.statistics.filt.in
        fi

        bash {script_dir}nonref_SD/scripts/endfilt.sh all1.end.statistics.filt.in rm.out rm_trf.s.bed other.end
        ### 删除掉现在T2T(或者是说现在已经有的SDpair)有的SDpair____这里看一下怎么从一个一直迭代添加的数据库中获取
        ## modify
        cat  <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/T2Thave.bed |cut -f 1-6)  <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/T2Thave.bed |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') >haved.pair
        bedtools intersect -a <(cat other.end |cut -f 1-12) -b <(cat haved.pair |cut -f 1-6) -f 0.8 -r -wa -wb|awk '$4==$16' >a
        bedtools intersect -a <(less -S a |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' ) -b <(less -S a |awk 'OFS="\\t"{{print $16,$17,$18,$13,$14,$15}}') -wa -wb  -f 0.8 -r|less -S |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' >lef.del.pair
        bedtools intersect -a <(cat other.end |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}') -b <(cat haved.pair |cut -f 1-6) -f 0.8 -r -wa -wb|awk '$4==$16' >a
        bedtools intersect -a <(less -S a |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' ) -b <(less -S a |awk 'OFS="\\t"{{print $16,$17,$18,$13,$14,$15}}') -wa -wb  -f 0.8 -r|less -S  |cut -f 1-12>right.del.pair
        grep -v -f <(cat lef.del.pair right.del.pair  |sort |uniq) <(cat  other.end |cut -f 1-12) >samuniq.pair
        cd ..
        cp ./processend/samuniq.pair {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/samuniq.pair
        """


rule gene_uniqpair:
    input:
        rules.filter_process.output
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/otherdfaddsam.txt"
    shell:
        '''
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/processend
        cat other.end |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}' |less -S >uniqpair.bed
        less -S {base_dir}03_refSD/sam_cor/Wpath/{wildcards.line}.W |cut -f 4,5,6 |less -S >x
        bedtools intersect -a <(less uniqpair.bed  |bedtools sort |bedtools merge) -b x -wa -wb |cut -f 4-6 |less |sort |uniq |less >process.W
        less -S uniqpair.bed  |bedtools sort  |uniq >../uniqpair.sort
        rm -f *.tsv
        rm -f otherinter.txt
        python {script_dir}nonref_SD/scripts/find_W.py --W {base_dir}03_refSD/sam_cor/Wpath/{wildcards.line}.W --d {base_dir}/03_refSD/sam_cor/seg_db
        cat *.tsv >{wildcards.line}.chooW
        rm -f *.tsv
        less -S {wildcards.line}.chooW |grep -v "node" |tr "@" "\\t" |awk 'OFS="\\t"{{print $7,$4+$8,$5+$8,$1,$2}}' |less -S > {wildcards.line}.Wallposend.tsv
        bedtools intersect -a uniqpair.bed -b {wildcards.line}.Wallposend.tsv -wa -wb|awk 'OFS="\\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' >otherinter.txt
        Rscript {script_dir}nonref_SD/scripts/SDgenepath.r otherinter.txt otherdf.txt
        awk -v d={wildcards.line} '{{print d "\\t" $0}}' otherdf.txt >otherdfaddsam.txt
        '''


rule merge_pair:
    input:
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/otherdfaddsam.txt", line=lines,base_dir=base_dir)
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/df.pair"
    shell:
        '''
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/
        cat {base_dir}/03_nonref/02_ITER/2_SD_GENE/*/processend/otherdfaddsam.txt >allpair
        less -S allpair   |sort -k5,5 >allpair.sort
        python {script_dir}nonref_SD/scripts/flter.py --p {base_dir}

        '''

rule merge_split_pair:
    input:
        rules.merge_pair.output
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/y.temp"
    shell:
        '''
        cd  {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/processend
        if [ -s filt_split.path ]; then
            cat filt_split.path  |awk 'OFS="\\t"{{print $12,$10,$11}}' >x.temp
            bedtools intersect -a x.temp -b x.temp -wa -wb -f 0.95 -r  |awk '$2!=$5 && $3!=$6' >y.temp
        else
            touch filt_split.path
            touch y.temp
        fi
        
        '''


rule filt_pair:
    input:
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/y.temp", line=lines,base_dir=base_dir)
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/allpair.sort.END.e"
    shell:
        '''
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/
        python {script_dir}nonref_SD/scripts/filter2.py --b {base_dir} --s {sam_list}
        less -S allpair.sort.END |cut -f 2-5>allpair.sort.END.e
        '''


rule gene_simapa2:
    input:
        W="{base_dir}/03_refSD/sam_cor/Wpath/{sam}.W",  # WCOM=rules.gene_Wpath.output,
        pain="{base_dir}/03_nonref/02_ITER/2_SD_GENE/allpair.sort.END.e"
    output:
        output_bed="{base_dir}03_nonref/02_ITER/3_PAIR_PATH/output/{sam}can_sd.corregion.bed",  # 目标输出文件
        time_output="{base_dir}03_nonref/02_ITER/3_PAIR_PATH/time/{sam}.output.txt",  # 时间输出
        time_log="{base_dir}03_nonref/02_ITER/3_PAIR_PATH/time/{sam}.txt",  # 日志输出
    params:
        dbpath="{base_dir}/03_refSD/sam_cor/"  # 数据库路径
    shell:
        '''
        cd {base_dir}/03_nonref/02_ITER/3_PAIR_PATH
        /usr/bin/time -v python {script_dir}/ref_SD/sam_sd/scripts/iterpansd.py --W {input.W} --S {wildcards.sam} --i {input.pain} --o {output.output_bed} --p {params.dbpath} > {output.time_output} 2> {output.time_log} ##这里面有个bedtools路径需要改 iterpansd.py 也可以改为两步骤 iterpansd1.py iterpansd2.py  
        '''

rule pair_trans:
    input:
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/other.end",base_dir=base_dir,line=lines),
        pain="{base_dir}/03_nonref/02_ITER/2_SD_GENE/allpair.sort.END.e"
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/{line}.SD.trans" 
    shell:
        '''
        cd {base_dir}/03_nonref/02_ITER/2_SD_GENE/{wildcards.line}/processend
        cat other.end  |awk '{{print $1":"$2"-"$3"\\t"$4":"$5"-"$6}}' >{wildcards.line}.SD.pair
        python {script_dir}nonref_SD/scripts/trans_pair.py --i {wildcards.line}.SD.pair --tran {base_dir}/03_nonref/02_ITER/2_SD_GENE/df.pair.END   --o {wildcards.line}.SD.trans

        
        '''

rule gene_uniqtran:
    input:
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/{line}.SD.trans",base_dir=base_dir,line=lines),
    output:
        "{base_dir}/03_nonref/02_ITER/2_SD_GENE/tranall.uniq"
    run:
        import subprocess
        import pandas as pd
        commo="cat *.trans |sort |uniq >"+base_dir+"/03_nonref/02_ITER/2_SD_GENE/tranall"
        T2Tpos = pd.read_csv(base_dir+"/03_nonref/02_ITER/2_SD_GENE/tranall",sep="\t",header=None)
        T2Tpos['sorted_pair'] = T2Tpos.apply(lambda row: tuple(sorted([row[0], row[1]])), axis=1)
        T2Tpos_unique = T2Tpos.drop_duplicates(subset='sorted_pair').drop(columns='sorted_pair')
        T2Tpos_unique.to_csv(base_dir+"/03_nonref/02_ITER/2_SD_GENE/tranall.uniq", index=False, header=False,sep="\t")



# rule gene_align:
#     input:
#         W=expand("{base_dir}nonreference/Wpath/{sam}.W",base_dir=base_dir,sam=ALLSAM)  # WCOM=rules.gene_Wpath.output,
#         corpath=expand("{base_dir}nonreference/01_ITER/3_PAIR_PATH/output/{sam}can_sd.corregion.bed",base_dir=base_dir,sam=ALLSAM)
#     output:
#     shell:
#         '''
#         mkdir -p {base_dir}/nonreference/01_ITER/2_SD_GENE/{wildcards.sam}/ALIGN
#         cat {input.W} |cut -f 2-6 |awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' |less -S >{wildcards.sam}.init.W
#         cat "{base_dir}nonreference/01_ITER/3_PAIR_PATH/output/{wildcards.sam}can_sd.corregion.bed" |tail -n +2| sed \"s/_rows_/\t/\" | less -S >{wildcards.sam}.pos
#         python /share/home/zhanglab/user/maoyafei/project/pansd/trans_pair.py --pos {wildcards.sam}.pos --W {wildcards.sam}.init.W --o {wildcards.sam}.uniqpair
#         less -S {wildcards.sam}.uniqpair |tr \":\" \"\t\" |awk 'BEGIN{OFS=\"\t\"} {gsub(\"-\", \"\t\", \$2); gsub(\"-\", \"\t\", \$4); print}' -  |cut -f 1-6 >{wildcards.sam}.minimapin
#         less -S {wildcards.sam}.minimapin |sort |uniq |awk '!(\$1==\$4 && ((\$5<=\$2 && \$6>=\$2)||(\$5<=\$3 && \$6>=\$3)||(\$2<=\$5 && \$3>=\$5)||(\$2<=\$6 && \$3>=\$6)))'  >{wildcards.sam}.minimapin1
#         fa=\"/share/home/zhanglab/user/maoyafei/allfa/$sam.fa\"
#         less -S $sam.minimapin1 |cut -f 1-3 |sort |uniq >$sam.left
#         rm -f $sam.paf
#         while IFS=$'\t' read -r col1 col2 col3
#         do
#             awk -F\"\t\" -v col1=\"\$col1\" -v col2=\"\$col2\" -v col3=\"\$col3\" '\$1==col1 && \$2==col2 && \$3==col3' $sam.minimapin1 > $sam.right
#             cat $sam.right
#             minimap2 -c --eqx <(samtools faidx \$fa \"\$col1:\$col2-\$col3\") <(samtools faidx \$fa \$(cat $sam.right |awk '{{print \$4\":\"\$5\"-\"\$6}}'))   -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >> $sam.paf
#         done <$sam.left
#         python /share/home/zhanglab/user/maoyafei/project/pansd/zdcopy/paf2bed.py  --paf  $sam.paf  --o  $sam.all1.end.statistics
#         '''





