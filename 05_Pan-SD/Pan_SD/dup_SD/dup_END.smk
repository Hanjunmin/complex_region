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

# script_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/Pan_SD/"
# base_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/run/"
# DATA_dir="/home/jmhan/project/APG/github_final/PanSD_endtest/example/"
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
SAMRMTRF=base_dir+"/02_refSD/sam_cor/2_ALIGN/" ###每个样本对应的RM和TRF的区域

# rule all:
#     input:
#         f"{base_dir}/03_nonref/non_ref.OK"


rule filter_process: ##yes ###算啦好像本来生成的新的就很少，就先这样不cluster了
    input:
        rules.merge_end_samSD.output
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/samuniq.pair",
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/other.end"
    threads: 1  
    shell:
        """
        mkdir -p {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/processend ##modify
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/processend
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
        cat  <(cat {base_dir}/02_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/T2Thave.bed |cut -f 1-6)  <(cat {base_dir}/02_refSD/sam_cor/2_ALIGN/{wildcards.line}/t2tsdminimap/T2Thave.bed |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') >haved.pair
        bedtools intersect -a <(cat other.end |cut -f 1-12) -b <(cat haved.pair |cut -f 1-6) -f 0.8 -r -wa -wb|awk '$4==$16' >a
        bedtools intersect -a <(less -S a |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' ) -b <(less -S a |awk 'OFS="\\t"{{print $16,$17,$18,$13,$14,$15}}') -wa -wb  -f 0.8 -r|less -S |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' >lef.del.pair
        bedtools intersect -a <(cat other.end |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}') -b <(cat haved.pair |cut -f 1-6) -f 0.8 -r -wa -wb|awk '$4==$16' >a
        bedtools intersect -a <(less -S a |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' ) -b <(less -S a |awk 'OFS="\\t"{{print $16,$17,$18,$13,$14,$15}}') -wa -wb  -f 0.8 -r|less -S  |cut -f 1-12>right.del.pair
        grep -v -f <(cat lef.del.pair right.del.pair  |sort |uniq) <(cat  other.end |cut -f 1-12) >samuniq.pair
        cd ..
        cp ./processend/samuniq.pair {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/samuniq.pair
        """

rule filter_otherend:
    input:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/other.end"
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/other.end.filt"
    shell:
        '''
        cd "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/processend"
    
        cat other.end |cut -f 1-6 >x
        cat x |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}' >y

        bedtools intersect -a x -b x -wa -wb |awk '$4==$10' |awk '$5<$12 && $6>$11' |awk '!($2==$8 && $3==$9 && $5==$11 && $6==$12)' >inter.x
        bedtools intersect -a y -b x -wa -wb |awk '$4==$10' |awk '$5<$12 && $6>$11' |awk '!($2==$8 && $3==$9 && $5==$11 && $6==$12)' >inter.y

        cat inter.x |awk 'OFS="\\t"{{print $1":"$2"-"$3"@"$4":"$5"-"$6,$7":"$8"-"$9"@"$10":"$11"-"$12}}' |sort |uniq >inter.reg.x
        cat inter.y |awk 'OFS="\\t"{{print $4":"$5"-"$6"@"$1":"$2"-"$3,$7":"$8"-"$9"@"$10":"$11"-"$12}}' |sort |uniq >inter.reg.y
        cat inter.reg.x inter.reg.y |sort |uniq >inter.reg.z
        if [ -s inter.reg.z ]; then
            Rscript {script_dir}/nonref_SD/scripts/filter_SD.r
        else
            cat other.end >other.end.filt
        fi
        

        '''
    




rule gene_uniqpair:  ##yes
    input:
        rules.filter_process.output,
        rules.filter_otherend.output
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/otherdfaddsam.txt"
    shell:
        '''
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/processend
        cat other.end.filt |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}' |less -S >uniqpair.bed
        ### filter


        less -S {base_dir}02_refSD/sam_cor/Wpath/{wildcards.line}.W |cut -f 4,5,6 |less -S >x
        bedtools intersect -a <(less uniqpair.bed  |bedtools sort |bedtools merge) -b x -wa -wb |cut -f 4-6 |less |sort |uniq |less >process.W
        less -S uniqpair.bed  |bedtools sort  |uniq >../uniqpair.sort
        rm -f *.tsv
        rm -f otherinter.txt
        python {script_dir}nonref_SD/scripts/find_W.py --W {base_dir}02_refSD/sam_cor/Wpath/{wildcards.line}.W --d {base_dir}/02_refSD/sam_cor/seg_db
        cat *.tsv >{wildcards.line}.chooW
        rm -f *.tsv
        less -S {wildcards.line}.chooW |grep -v "node" |tr "@" "\\t" |awk 'OFS="\\t"{{print $7,$4+$8,$5+$8,$1,$2}}' |less -S > {wildcards.line}.Wallposend.tsv
        bedtools intersect -a uniqpair.bed -b {wildcards.line}.Wallposend.tsv -wa -wb|awk 'OFS="\\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' >otherinter.txt
        Rscript {script_dir}nonref_SD/scripts/SDgenepath.r otherinter.txt otherdf.txt
        awk -v d={wildcards.line} '{{print d "\\t" $0}}' otherdf.txt >otherdfaddsam.txt
        '''


rule merge_pair: ##yes
    input:
        expand("{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/otherdfaddsam.txt", line=lines,base_dir=base_dir)
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/df.pair"
    shell:
        '''
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/
        cat {base_dir}/03_nonref/03_DUP/2_SD_GENE/*/processend/otherdfaddsam.txt >allpair
        less -S allpair   |sort -k5,5 >allpair.sort
        python {script_dir}nonref_SD/scripts/flter.py --p {base_dir}

        '''

rule merge_split_pair: ##yes
    input:
        rules.merge_pair.output
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/y.temp"
    shell:
        '''
        cd  {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/processend
        if [ -s filt_split.path ]; then
            cat filt_split.path  |awk 'OFS="\\t"{{print $12,$10,$11}}' >x.temp
            bedtools intersect -a x.temp -b x.temp -wa -wb -f 0.95 -r  |awk '$2!=$5 && $3!=$6' >y.temp
        else
            touch filt_split.path
            touch y.temp
        fi
        
        '''


rule filt_pair: ##yes
    input:
        expand("{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/y.temp", line=lines,base_dir=base_dir)
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/allpair.sort.END.e"
    shell:
        '''
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/
        python {script_dir}nonref_SD/scripts/filter2.py --b {base_dir} --s {sam_list}
        less -S allpair.sort.END |cut -f 2-5 |sort |uniq >allpair.sort.END.e
        '''


rule gene_simapa2: 
    input:
        W="{base_dir}/02_refSD/sam_cor/Wpath/{sam}.W",  # WCOM=rules.gene_Wpath.output,
        pain="{base_dir}/03_nonref/03_DUP/2_SD_GENE/allpair.sort.END.e"
    output:
        output_bed="{base_dir}03_nonref/03_DUP/3_PAIR_PATH/output/{sam}can_sd.corregion.bed",  # 目标输出文件
        time_output="{base_dir}03_nonref/03_DUP/3_PAIR_PATH/time/{sam}.output.txt",  # 时间输出
        time_log="{base_dir}03_nonref/03_DUP/3_PAIR_PATH/time/{sam}.txt",  # 日志输出
    params:
        dbpath="{base_dir}/02_refSD/sam_cor/"  # 数据库路径
    shell:
        '''
        cd {base_dir}/03_nonref/03_DUP/3_PAIR_PATH
        python {script_dir}/ref_SD/sam_sd/scripts/iterpansd.py --W {input.W} --S {wildcards.sam} --i {input.pain} --o {output.output_bed} --p {params.dbpath} > {output.time_output} 2> {output.time_log} ##这里面有个bedtools路径需要改 iterpansd.py 也可以改为两步骤 iterpansd1.py iterpansd2.py  
        '''

rule pair_trans: ##yes
    input:
        expand("{base_dir}/03_nonref/03_DUP/2_SD_GENE/{line}/processend/other.end.filt",base_dir=base_dir,line=lines),
        pain="{base_dir}/03_nonref/03_DUP/2_SD_GENE/allpair.sort.END.e"
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/ALLPAIR/{line}.SD.trans" 
    shell:
        '''
        mkdir -p {base_dir}/03_nonref/03_DUP/2_SD_GENE/ALLPAIR/
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/ALLPAIR/
        cat {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.line}/processend/other.end.filt  |awk '{{print $1":"$2"-"$3"\\t"$4":"$5"-"$6}}' >{wildcards.line}.SD.pair
        python {script_dir}nonref_SD/scripts/trans_pair.py --i {wildcards.line}.SD.pair --tran {base_dir}/03_nonref/03_DUP/2_SD_GENE/df.pair.END   --o {wildcards.line}.SD.trans

        
        '''

rule gene_uniqtran:
    input:
        expand("{base_dir}/03_nonref/03_DUP/2_SD_GENE/ALLPAIR/{line}.SD.trans",base_dir=base_dir,line=lines),
    output:
        f"{base_dir}/03_nonref/03_DUP/2_SD_GENE/tranall.uniq"
    run:
        import subprocess
        import pandas as pd
        commo="cat "+base_dir+"/03_nonref/03_DUP/2_SD_GENE/ALLPAIR/*.trans |sort |uniq >"+base_dir+"/03_nonref/03_DUP/2_SD_GENE/tranall"
        subprocess.run(commo, shell=True)
        T2Tpos = pd.read_csv(base_dir+"/03_nonref/03_DUP/2_SD_GENE/tranall",sep="\t",header=None)
        T2Tpos['sorted_pair'] = T2Tpos.apply(lambda row: tuple(sorted([row[0], row[1]])), axis=1)
        T2Tpos_unique = T2Tpos.drop_duplicates(subset='sorted_pair').drop(columns='sorted_pair')
        T2Tpos_unique.to_csv(base_dir+"/03_nonref/03_DUP/2_SD_GENE/tranall.uniq", index=False, header=False,sep="\t")



rule gene_align:
    input:
        W="{base_dir}/02_refSD/sam_cor/Wpath/{sam}.W",
        corpath="{base_dir}03_nonref/03_DUP/3_PAIR_PATH/output/{sam}can_sd.corregion.bed",
        staal=f"{base_dir}/03_nonref/03_DUP/2_SD_GENE/tranall.uniq"
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{sam}/ALIGN/{sam}.all1.end.statistics"
    shell:
        '''
        mkdir -p {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.sam}/ALIGN
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.sam}/ALIGN
        cat {input.W} |cut -f 2-6 |awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' |less -S >{wildcards.sam}.init.W
        cat "{base_dir}03_nonref/03_DUP/3_PAIR_PATH/output/{wildcards.sam}can_sd.corregion.bed" |tail -n +2| sed "s/_rows_/\\t/" | less -S >{wildcards.sam}.pos
        python {script_dir}nonref_SD/scripts/trans_pair_w.py --pos {wildcards.sam}.pos --W {wildcards.sam}.init.W --o {wildcards.sam}.uniqpair --trans "{base_dir}/03_nonref/03_DUP/2_SD_GENE/tranall.uniq"
        less -S {wildcards.sam}.uniqpair |tr ":" "\\t" |awk 'BEGIN{{OFS="\\t"}} {{gsub("-", "\\t", $2); gsub("-", "\\t", $4); print}}' -  |cut -f 1-6 >{wildcards.sam}.minimapin
        less -S {wildcards.sam}.minimapin |sort |uniq |awk '!($1==$4 && (($5<=$2 && $6>=$2)||($5<=$3 && $6>=$3)||($2<=$5 && $3>=$5)||($2<=$6 && $3>=$6)))'  >{wildcards.sam}.minimapin1
        fa={DATA_dir}/FA/{wildcards.sam}.fa
        less -S {wildcards.sam}.minimapin1 |cut -f 1-3 |sort |uniq >{wildcards.sam}.left
        rm -f {wildcards.sam}.paf
        while IFS=$'\\t' read -r col1 col2 col3
        do
            awk -F "\\t" -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' {wildcards.sam}.minimapin1 > {wildcards.sam}.right
            cat {wildcards.sam}.right
            minimap2 -c --eqx <(samtools faidx $fa "$col1:$col2-$col3") <(samtools faidx $fa $(cat {wildcards.sam}.right |awk '{{print $4":"$5"-"$6}}'))   -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >> {wildcards.sam}.paf
        done <{wildcards.sam}.left
        python {script_dir}/nonref_SD/scripts/paf2bed.py  --paf  {wildcards.sam}.paf  --o {wildcards.sam}.all1.end.statistics
        '''


rule run_rm_sam:
    input:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{sam}/ALIGN/{sam}.all1.end.statistics"
    output:
        "{base_dir}/03_nonref/03_DUP/2_SD_GENE/{sam}/ALIGN/{sam}.eebed.rev.out",
        "{base_dir}03_nonref/03_DUP/2_SD_GENE/{sam}/ALIGN/{sam}.ITER.SD.al"
    shell:
        '''
        cd {base_dir}/03_nonref/03_DUP/2_SD_GENE/{wildcards.sam}/ALIGN
        [ -d "lefdic" ] && rm -r lefdic
        mkdir lefdic
        [ -d "rigdic" ] && rm -r rigdic
        mkdir rigdic
        samfa={DATA_dir}/FA/{wildcards.sam}.fa
        cat {wildcards.sam}.all1.end.statistics |tail -n +2 |cut -f 1-3 >./lefdic/lef.bed
        cat {wildcards.sam}.all1.end.statistics |tail -n +2|cut -f 4-6 >./rigdic/rig.bed
        cd lefdic
        bash {script_dir}nonref_SD/scripts/run_rm_trf.sh lef.bed $samfa
        cd ..
        cd rigdic
        bash {script_dir}nonref_SD/scripts/run_rm_trf.sh rig.bed $samfa
        cd ..
        cat lefdic/rm_trf.s.bed rigdic/rm_trf.s.bed  |bedtools sort |bedtools merge >rm_trf.ee.bed
        cat lefdic/rm.out rigdic/rm.out >rm.ee.out
        cat {wildcards.sam}.all1.end.statistics |cut -f 1-12 >{wildcards.sam}.eebed 
        bash {script_dir}nonref_SD/scripts/endfilt.sh {wildcards.sam}.eebed rm.ee.out rm_trf.ee.bed {wildcards.sam}.eebed.rev.out
        cd {base_dir}03_nonref/03_DUP/2_SD_GENE/{wildcards.sam}/ALIGN/
        python {script_dir}/SD_SUM/scripts/trans_pair2.py --pos "{base_dir}03_nonref/03_DUP/2_SD_GENE/{wildcards.sam}/ALIGN/{wildcards.sam}.pos" --W {wildcards.sam}.init.W --o {wildcards.sam} --pair "{base_dir}03_nonref/03_DUP/2_SD_GENE/tranall.uniq"
        python {script_dir}/SD_SUM/scripts/construct_dict.py --pair "{base_dir}03_nonref/03_DUP/2_SD_GENE/tranall.uniq" --sam {wildcards.sam} --scp {script_dir}/ref_SD/sam_sd/scripts/filter_sam_SD.r
        
        '''


rule sta_end:
    input:
        expand("{base_dir}/03_nonref/03_DUP/2_SD_GENE/{sam}/ALIGN/{sam}.eebed.rev.out",base_dir=base_dir, sam=ALLSAM)
    output:
        "{base_dir}/03_nonref/non_ref.OK"
    shell:
        '''
        touch {output}
        '''



