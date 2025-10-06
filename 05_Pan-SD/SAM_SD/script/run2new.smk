import os

# Snakemake 命令示例
# snakemake -s .pangenome/Pan-SD/pansam/datafrommao/run2.smk --keep-going --keep-incomplete -j 10

# 定义输入和输出文件路径
cordipath = ".PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/output/"
sdregion = ".PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/alldfe.txt"
canregion = ".PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/alldfe1.txt"
sam_list =".PANSDEND/DATA/allsam"
T2Tedge=".PANSDEND/01_REFSD/edge.bed"
genomefa=".allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/"
T2Trmtrf=".HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/trf_rm.bed"
Wfold=".PANSDEND/DATA/APG_hgsvc_hprcW/"
rm=".HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/asm_repeatmasker.out.bed"
rmtrf=".HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/trf_rm.bed"
RMsoft="/home/Public/software/RepeatMasker-4.1.4/"

T2Tsate=".PANSDEND/DATA/CHM13.sat"
script=".PANSDEND/script_02_SAM_REF/",

nowdir=".PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN"


with open(sam_list) as f:
    lines = [line.strip() for line in f.readlines()]


rule all:
    input:
        #expand("{li}/all_values.txt", li=lines),
        expand("{li}/RMTRF/RMTRFSAT.region", li=lines),
        expand("{li}/{li}.sd", li=lines),
        expand("{li}/can/filin.bedend", li=lines)
        


rule init_sd_can:
    input:
        samnow=cordipath + "{li}can_sd.corregion.bed",
        samsim=cordipath + "{li}can_sd.corregion.bed.sim",
    output:
        output_bed="{li}/{li}.sd"
    params:
        sam_list=sam_list  # 添加 sam_list
    shell:
        """
        cd {nowdir}
        [ ! -d {wildcards.li} ] && mkdir -p {wildcards.li}
        echo {wildcards.li}
        cd {wildcards.li}
        less -S {Wfold}{wildcards.li}.W | cut -f 2-6 | awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' > init.W
        cat {input.samnow} | sed \"s/_rows_/\t/\" |tail -n +2 >{wildcards.li}.sd.pos
        Rscript {script}splidtcansd.r {sdregion} {canregion} {wildcards.li} {wildcards.li}.sd.pos init.W
        """


rule geneRMTRF:
    input:
        rules.init_sd_can.output,
    output:
        output_bed="{li}/RMTRF/RMTRFSAT.region"
    shell:
        """
        cd {nowdir}
        echo {wildcards.li}
        cd {wildcards.li}
        [ -d RMTRF ] && rm -r  RMTRF
        [ ! -d RMTRF ] && mkdir -p RMTRF
      cd RMTRF
      paste <(cat ../{wildcards.li}.sd  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.li}.sd  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}') <(cat ../{wildcards.li}.sd|cut -f5)  |tail -n +2 >t2tsdrm.simn
      paste <(cat ../{wildcards.li}.can  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.li}.can  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}') <(cat ../{wildcards.li}.can|cut -f5) |tail -n +2 >t2trm.simn
      cat t2tsdrm.simn  t2trm.simn |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7}}'>allrm
     awk '{{
    diff1 = $3 - $2; 
    diff2 = $6 - $5; 
    min = (diff1 < diff2) ? diff1 : diff2;  
    max = (diff1 > diff2) ? diff1 : diff2;  
    print $0"\\t"(max == 0 ? 0 : min / max);   
    }}' allrm >allrm.sim2

    less -S allrm.sim2 |awk '$7>=0.99 && $8>=0.99' >higncovn
    bedtools intersect -a {T2Trmtrf} -b  <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -wa -wb >rm_inter_higncov
    bedtools subtract -a <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -b {T2Trmtrf} -A  |cut -f 4-6|bedtools sort |bedtools merge >rm_inter_higncov.no
    Rscript .PANSDEND/script/cal_RMTRF.r
    mv ALL.HIGHCOV.RMTRF ./RMTRF.END
    cat <(cat RMTRF.END |cut -f 1-3 |bedtools sort |bedtools merge) <(cat rm_inter_higncov.no) >RMTRF.region ###这个大样本区域
    cat RMTRF.END |cut -f 4-6 |bedtools sort |bedtools merge >RMTRF.rmtrf ###这个大样本区域对应的重复序列情况

    bedtools intersect -a <(cat {T2Tsate}|cut -f 1-3) -b  <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -wa -wb >rm_inter_higncov
    bedtools subtract -a <(cat higncovn|awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3}}') -b <(cat {T2Tsate}|cut -f 1-3) -A|cut -f 4-6 |bedtools sort |bedtools merge >rm_inter_higncov.no
    Rscript .PANSDEND/script/cal_RMTRF.r
    mv ALL.HIGHCOV.RMTRF ./RMSAT.END
    cat <(cat RMSAT.END |cut -f 1-3 |bedtools sort |bedtools merge) <(cat rm_inter_higncov.no) >RMSAT.region ###这个大样本区域
    cat RMSAT.END |cut -f 4-6 |bedtools sort |bedtools merge >RMSAT.sat ###这个大样本区域对应的重复序列情况
    rm rm_inter_higncov.no
    rm allrm.sim2
    rm higncovn
    bedtools intersect -a RMSAT.region -b RMTRF.region >RMTRFSAT.region
    cd ..
    """



rule process_line1:
    input:
        samnow=cordipath + "{li}can_sd.corregion.bed",
        samsim=cordipath + "{li}can_sd.corregion.bed.sim",
        trfrm="{li}/RMTRF/RMTRFSAT.region"
    output:
        output_bed="{li}/filin.bedend"
    params:
        sam_list=sam_list  # 添加 sam_list
    shell:
        """
        cd {nowdir}
        cd {wildcards.li}
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
        less -S process.txt |grep "nosam" >processnosam.txt ##(这是第二种情况，是找不到对应的区域，相当于是变化极大)
        [ -d t2tsdminimap ] && rm -r  t2tsdminimap
        [ ! -d t2tsdminimap ] && mkdir -p t2tsdminimap
        cd t2tsdminimap
        fa=$(find {genomefa}{wildcards.li}*.fa)
        cat $fa >{wildcards.li}.fa
        samtools faidx {wildcards.li}.fa
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
        less other_T2T |awk 'OFS="\\t"{{print $2,$1}}' | sort | uniq | tr ":" "\\t" | awk 'BEGIN {{OFS="\\t"}} {{gsub("-", "\\t", $2); gsub("-", "\\t", $4); print}}' | less -S |grep -v "T2T" |less -S >t2trm
        bedtools intersect -a t2trm -b <(cat {rm} |cut -f 1-3) -wa -wb >groupin
        bedtools subtract -a t2trm -b <(cat {rm} |cut -f 1-3) -A |cut -f 4-6>groupin.ana
      Rscript {script}expand_rmtrf.r "groupin" "{wildcards.li}.fa.fai" "rm_in" 2 ##这里的2表示左右扩的时候多扩两倍
      cat "rm_in" groupin.ana >rm_in.end
      bedtools coverage -a rm_in.end -b ../RMTRF/RMTRFSAT.region |awk '$7!=1' >rmin.in
      bash {script}run_rm_trf.sh rmin.in {wildcards.li} {RMsoft} {script} ##这步为跑trf和rm
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
        samsim=cordipath + "{li}can_sd.corregion.bed.sim"
    output:
        output_bed="{li}/can/filin.bedend"
    params:
        sam_list=sam_list  # 添加 sam_list
    shell:
        """
        cd {nowdir}
        echo {wildcards.li}
        cd {wildcards.li}
        [ -d can ] && rm -r  can
        [ ! -d can ] && mkdir -p can
         
      cd can
      paste <(cat ../{wildcards.li}.sd  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.li}.sd  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}')  |less -S |tail -n +2 >t2tsdrm
      bedtools intersect -a t2tsdrm -b <(cat {T2Trmtrf} |cut -f 1-3) -wa -wb >groupinsd
      
      paste <(cat ../{wildcards.li}.can  |cut -f2|tr ":" "\\t" |tr "-" "\\t" ) <(cat ../{wildcards.li}.can  |cut -f 3-4,9-10|awk 'OFS="\\t"{{print $3,$4+$1,$4+$2}}')  |less -S |tail -n +2 >t2trm
      bedtools intersect -a t2trm -b <(cat {T2Trmtrf} |cut -f 1-3) -wa -wb >groupin
      bedtools subtract -a t2trm -b <(cat  {T2Trmtrf} |cut -f 1-3 |bedtools sort |bedtools merge) -A  |cut -f 4-6>groupin.ana
      Rscript {script}expand_rmtrf.r "groupin" "../t2tsdminimap/{wildcards.li}.fa.fai" "rm_in" 1
      bash {script}transe.sh rm_in rm_in.transe
      mv rm_in.transe rm_in
      cat ../t2tsdminimap/T2Thave.bed   ../t2tsdminimap/T2Tpolymor1.bed ../t2tsdminimap/filin.bedend |cut -f 1-12 >t2tcoresbon.end ##*****************
      bedtools subtract -a <(cat rm_in groupin.ana) -b <(cat t2tcoresbon.end |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}') >candidat.sub
      bedtools subtract -a <(cat  candidat.sub |bedtools sort |bedtools merge) -b  ../allt2tsdcoorend_true.tsv >candidat.sub2
      [ -f allSDsample.fa ] && rm allSDsample.fa
      fa="../t2tsdminimap/{wildcards.li}.fa"
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
bash /{script}run_rm_trf.sh rmin.in {wildcards.li} {RMsoft} {script} ##这步为跑trf和rm
cat <(cat ../RMTRF/RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat rm.out |cut -f 1-3,5 )  >allsat
cat rm_trf.s.bed ../RMTRF/RMTRF.rmtrf  >allrmtrf
bash {script}transe.sh allrmtrf allrmtrf.trane
bash {script}endfilt.sh allsat allrmtrf.trane all.end.statistics.filt filin.bedend
snakemake -s {script}cluste_end.smk --config fasta_file1={nowdir}/{wildcards.li}/t2tsdminimap/{wildcards.li}.fa script={script}  -j 1 ###filin.bedend 输入为这个但是有点冗余所以需要去cluster一下
rm {nowdir}/{wildcards.li}/t2tsdminimap/{wildcards.li}.fa
"""


