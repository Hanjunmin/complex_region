import os

# Snakemake 命令示例
# snakemake -s /dssg/home/acct-clsmyf/clsmyf-user1/pangenome/Pan-SD/pansam/datafrommao/run2.smk --keep-going --keep-incomplete -j 10

script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"
base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"
DATA_dir="/home/jmhan/project/APG/github_final/PanSD_end/example/"


# 定义输入和输出文件路径
samcanpath=base_dir+"/03_refSD/sam_cor/1_SD_CANPATH/"
cordipath = samcanpath+"/output/"
sdregion = samcanpath+"/03_refSD/sam_cor/1_SD_CANPATH/alldfe.txt"
canregion = samcanpath+"/alldfe1.txt"
T2Tedge=samcanpath+"sd.edge"
genomefa=DATA_dir+"FA/"
rm="/home/jmhan/HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/asm_repeatmasker.out.bed"
RMsoft="/home/Public/software/RepeatMasker-4.1.4/"
T2Tsate="/home/jmhan/PANSDEND/DATA/CHM13.sat"
script=script_dir+"/ref_SD/sam_sd/scripts/",
nowdir=base_dir+"03_refSD/sam_cor/2_ALIGN"


import os
## cd /home/jmhan/iterall/
#snakemake -s /home/jmhan/iterall/interdb.smk --keep-going --keep-incomplete -j 10
SAMLIS = DATA_dir+"samlis.txt"
import pandas as pd
samdf = pd.read_csv(SAMLIS,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
lines = samname

rule all:
    input:
        #expand("{li}/all_values.txt", li=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/RMTRFSAT.region",base_dir=base_dir, li=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/{li}.sd", li=lines,base_dir=base_dir),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/can/filin.bedend", li=lines,base_dir=base_dir),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/SDend", li=lines,base_dir=base_dir),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/1.OK", li=lines,base_dir=base_dir),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/2.OK", li=lines,base_dir=base_dir),
        f"{base_dir}/03_refSD/sam_cor/2_ALIGN/align.OK"
        


rule init_sd_can:
    input:
        samnow="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{li}can_sd.corregion.bed",
        samsim="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{li}can_sd.corregion.bed.sim",
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/{li}.sd"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        less -S {base_dir}/03_refSD/sam_cor/Wpath/{wildcards.li}.W | cut -f 2-6 | awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' > init.W
        cat {input.samnow} | sed \"s/_rows_/\\t/\" |tail -n +2 >{wildcards.li}.sd.pos
        Rscript {script}splidtcansd.r {sdregion} {canregion} {wildcards.li} {wildcards.li}.sd.pos init.W
        """


rule geneRMTRF:
    input:
        rules.init_sd_can.output,
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/RMTRFSAT.region"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        [ -d RMTRF ] && rm -r  RMTRF
        [ ! -d RMTRF ] && mkdir -p RMTRF
      cd RMTRF
      paste <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}/{wildcards.li}.sd  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}/{wildcards.li}.sd  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}') <(cat ../{wildcards.li}.sd|cut -f5)  |tail -n +2 >t2tsdrm.simn
      paste <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}/{wildcards.li}.can  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}/{wildcards.li}.can  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}') <(cat ../{wildcards.li}.can|cut -f5) |tail -n +2 >t2trm.simn
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



rule process_line1:
    input:
        samnow="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{li}can_sd.corregion.bed",
        samsim="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{li}can_sd.corregion.bed.sim",
        trfrm="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/RMTRFSAT.region"
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/filin.bedend",
        updaten="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/t2tsdminimap/rm_trf.s.bed"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        less -S {wildcards.li}.sd | tail -n +2 | awk 'OFS="\\t"{{print $7"#"$9,$10+$3,$10+$4}}' > allt2tsdcoorend.tsv
        less -S {wildcards.li}.sd |tail -n +2 |awk 'OFS="\\t"{{print $9,$10+$3,$10+$4}}' |less -S > allt2tsdcoorend_true.tsv
        less -S {wildcards.li}.can |tail -n +2 |awk 'OFS="\\t"{{print $7"#"$9,$10+$3,$10+$4}}' |less -S >candidateall.tsv
        less -S {wildcards.li}.can |tail -n +2 |awk 'OFS="\\t"{{print $9,$10+$3,$10+$4}}' |less -S >candidateall_true.tsv
        Rscript {script}T2Tmatch.r {T2Tedge} {wildcards.li}.sd
        cat <(less -S process.txt |grep -v "nosam" |tail -n +2 |cut -f1,2 |tr "\\t" "\\n" |less -S |sort |uniq) <(cat end.txt |tail -n +2|cut -f1,2 |tr "\\t" "\\n") |less -S >reprocess.pair
        cat <(less -S process.txt|grep -v "nosam" |tail -n +2 |cut -f1,2) <(cat end.txt |tail -n +2 |cut -f1,2) |less -S >reprocess.pairall
        grep -Ff <(cut -f1 reprocess.pair) {wildcards.li}.sd > {wildcards.li}.multicoordi.inte.tsv
        less -S {wildcards.li}.multicoordi.inte.tsv |awk 'OFS="\\t"{{print $2,$7"@"$9"@"$10"@"$11"@"$3"@"$4}}' |less -S  >cor.id
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
        fa=$(find {genomefa}{wildcards.li}*.fa)
        cat $fa >{wildcards.li}.fa
        samtools faidx {wildcards.li}.fa
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
        Rscript {script}expand_rmtrf.r "groupin" "{wildcards.li}.fa.fai" "rm_in" 2 ##这里的2表示左右扩的时候多扩两倍
        cat "rm_in" groupin.ana >rm_in.end
        bedtools coverage -a rm_in.end -b ../RMTRF/RMTRFSAT.region |awk '$7!=1' >rmin.in
        bash {script}run_rm_trf.sh rmin.in {wildcards.li} {RMsoft} {script} ##这步为跑trf和rm  
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
      cp ./{wildcards.li}/t2tsdminimap/filin.bedend {output.output_bed}
       
       """



#####很重要的三个文件！！！RMTRFSAT.region RMSAT.sat RMTRF.rmtrf 
rule process_line2:
    input:
        rules.process_line1.output,
        samsim="{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/output/{li}can_sd.corregion.bed.sim"
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/can/filin.bedend",
        updaten="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/can/rm_trf.s.bed"
    shell:
        """
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        [ -d can ] && rm -r  can
        [ ! -d can ] && mkdir -p can
         
      cd can
      paste <(cat ../{wildcards.li}.sd  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.li}.sd  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}')  |less -S |tail -n +2 >t2tsdrm
      bedtools intersect -a t2tsdrm -b <(cat {base_dir}/01_intersim/rm_trf.bed |cut -f 1-3) -wa -wb >groupinsd
      echo "A"
      paste <(cat ../{wildcards.li}.can  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.li}.can  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}')  |less -S |tail -n +2 >t2trm
      bedtools intersect -a t2trm -b <(cat {base_dir}/01_intersim/rm_trf.bed |cut -f 1-3) -wa -wb >groupin
      bedtools subtract -a t2trm -b <(cat  {base_dir}/01_intersim/rm_trf.bed |cut -f 1-3 |bedtools sort |bedtools merge) -A  |cut -f 4-6>groupin.ana
      Rscript {script}expand_rmtrf.r "groupin" "../t2tsdminimap/{wildcards.li}.fa.fai" "rm_in" 1
      bash {script}transe.sh rm_in rm_in.transe
      mv rm_in.transe rm_in
      cat ../t2tsdminimap/T2Thave.bed   ../t2tsdminimap/T2Tpolymor1.bed ../t2tsdminimap/filin.bedend |cut -f 1-12 >t2tcoresbon.end ##*****************
      bedtools subtract -a <(cat rm_in groupin.ana) -b <(cat t2tcoresbon.end |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}') >candidat.sub
      bedtools subtract -a <(cat  candidat.sub |bedtools sort |bedtools merge) -b  ../allt2tsdcoorend_true.tsv >candidat.sub2
      [ -f allSDsample.fa ] && rm allSDsample.fa
      fa="../t2tsdminimap/{wildcards.li}.fa"
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
    cp ../t2tsdminimap/{wildcards.li}.fa ./
     echo "C"

    if [ -s "rmin.in" ]; then
        bash {script}run_rm_trf.sh rmin.in {wildcards.li} {RMsoft} {script} ##这步为跑trf和rm
    else
        touch rm.out
        touch rm_trf.s.bed
    fi

    echo "D"
    cat <(cat ../RMTRF/RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat rm.out |cut -f 1-3,5 )  >allsat
    cat rm_trf.s.bed ../RMTRF/RMTRF.rmtrf  >allrmtrf
    bash {script}transe.sh allrmtrf allrmtrf.trane
    bash {script}endfilt.sh allsat allrmtrf.trane all.end.statistics.filt filin.bedend
    # echo "E1"
    # if [ -s "filin.bedend" ]; then
    # echo "AAA"
    # mkdir -p endfilt1 && cd endfilt1
    # touch SDend
    # else
    # snakemake -s {script}cluste_end.smk --config fasta_file1={nowdir}/{wildcards.li}/t2tsdminimap/{wildcards.li}.fa script={script}  -j 1 ###filin.bedend 输入为这个但是有点冗余所以需要去cluster一下
    # fi
    
    # echo "E"
    

    
    # rm {nowdir}/{wildcards.li}/t2tsdminimap/{wildcards.li}.fa
"""




rule merge_end:
    input:
        rules.process_line2.output,
    output:
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/SDend"
    threads: 1
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}
        mkdir -p endfilt && cd endfilt
        rm -f *
        if [ -s ../can/filin.bedend ]; then
            Rscript {script}ednfilt.r ../can/filin.bedend  {nowdir}/{wildcards.li}/t2tsdminimap/{wildcards.li}.fa {script}"paf2bed.py"
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
            Rscript {script}endfilt2.r panSD_T2T.bed  {nowdir}/{wildcards.li}/t2tsdminimap/{wildcards.li}.fa {script}"paf2bed.py"
            ls |grep "csv"|less -S |tr "_" "\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
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



#### updata_rmtrf

rule update_rmtrf_ref:
    input:
        rmnow="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/t2tsdminimap/rm_trf.s.bed"
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/1.OK"
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}/RMTRF
        cat <(cat {nowdir}/{wildcards.li}/t2tsdminimap/rmin.in |cut -f 1-3)  RMTRFSAT.region  >RMTRFSAT.region.temp
        cat {nowdir}/{wildcards.li}/t2tsdminimap/rm_trf.s.bed  RMTRF.rmtrf >RMTRF.rmtrf.temp
        cat <(cat RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat {nowdir}/{wildcards.li}/t2tsdminimap/rm.out |cut -f 1-3,5 )  >RMSAT.sat.temp
        mv RMTRFSAT.region.temp  RMTRFSAT.region
        mv RMTRF.rmtrf.temp  RMTRF.rmtrf
        mv RMSAT.sat.temp  RMSAT.sat
        touch 1.OK
        """


rule update_rmtrf_can:
    input:
        rmnow="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/can/rm_trf.s.bed"
    output:
        output_bed="{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/RMTRF/2.OK"
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/{wildcards.li}/RMTRF
        cat <(cat {nowdir}/{wildcards.li}/can/rmin.in |cut -f 1-3)  RMTRFSAT.region  >RMTRFSAT.region.temp
        cat {nowdir}/{wildcards.li}/can/rm_trf.s.bed  RMTRF.rmtrf >RMTRF.rmtrf.temp
        cat <(cat RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat {nowdir}/{wildcards.li}/can/rm.out |cut -f 1-3,5 )  >RMSAT.sat.temp
        mv RMTRFSAT.region.temp  RMTRFSAT.region
        mv RMTRF.rmtrf.temp  RMTRF.rmtrf
        mv RMSAT.sat.temp  RMSAT.sat
        touch 2.OK
        """


rule update_can_reg:
    input:
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/can/rm_in",
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/can/groupin.ana",
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/t2tsdminimap/T2Thave.bed",
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/{li}/t2tsdminimap/filin.bedend"
    output:
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/other/{li}.other",
        "{base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD/{li}.t2tsd"
    shell:
        """
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/other
        mkdir -p {base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD
        bedtools subtract -a <(cat {wildcards.li}/can/rm_in {wildcards.li}/can/groupin.ana) -b <(cat {wildcards.li}/t2tsdminimap/T2Thave.bed   {wildcards.li}/t2tsdminimap/T2Tpolymor1.bed {wildcards.li}/t2tsdminimap/filin.bedend {wildcards.li}/SDend |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}') |bedtools sort |bedtools merge >./other/{wildcards.li}".other"
        cat {wildcards.li}/t2tsdminimap/T2Thave.bed   {wildcards.li}/t2tsdminimap/T2Tpolymor1.bed {wildcards.li}/t2tsdminimap/filin.bedend {wildcards.li}/SDend >./T2TPANSD/{wildcards.li}".t2tsd"
        """


rule check_align_status:
    input:
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/other/{li}.other", li=lines,base_dir=base_dir),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD/{li}.t2tsd", li=lines,base_dir=base_dir)
    output:
        f"{base_dir}/03_refSD/sam_cor/2_ALIGN/align.OK"
    shell:
        '''
        cd {base_dir}/03_refSD/sam_cor/2_ALIGN/
        touch "align.OK"
        '''
