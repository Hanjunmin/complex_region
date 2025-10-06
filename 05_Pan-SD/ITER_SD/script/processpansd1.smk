import os

# Snakemake 命令示例
# cd ./pansdall/run/
# snakemake -s ./pansdall/processpansd1.smk --keep-going --keep-incomplete -j 10
## find ./pansdall/run/ -type f -path "*/filt/beforegroup.filt.txt" |less -S |tr "/" "\t" |less -S |cut -f6 |less -S >have
#### ./pansdall/processpansd1.smk
# 定义输入和输出文件路径
sam_list ="./PANSDEND/DATA/allsam"
# sam_list = "./pangenome/WGSfea/allW/APG3allWnew/allsam"
windir="./PANSDEND/DATA/APG_hgsvc_hprcW/"
win_simdir="./PANSDEND/02_SAMSD/02_02-ITERSD/00_gene_sim/winsimout/" #"./pangenome/WGSfea/allW/APG3allWnew/scripts/winsimout/"
minimapindir="./PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/other/" #"./pansdall/minimapindir/"
fadir="./allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/"
nowdir="./PANSDEND/02_SAMSD/02_02-ITERSD/01_iter/"
t2trmtrf="./PANSDEND/DATA/t2trmtrf.cut1M.bed"
SAMRMTRF="./PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/" ###每个样本对应的RM和TRF的区域

with open(sam_list) as f:
    lines = [line.strip() for line in f.readlines()]

rule all:
    input:
        expand("{li}/newin", li=lines),
        expand("{li}/filt/beforegroup.txt", li=lines),
        expand("{li}/SPLIT.OK", li=lines),
        expand("{li}/expand.out", li=lines)


rule process_line1:
    input:
        samwin=windir + "{li}.W",
        sammiin=minimapindir + "{li}.other",
    output:
        output_bed="{li}/newin"
    shell:
        """
        [ ! -d {wildcards.li} ] && mkdir -p {wildcards.li}
        echo {wildcards.li}
        cd {wildcards.li}
        less -S {windir}{wildcards.li}.W | cut -f 2-6 | awk 'OFS="\\t"{{print $0,$1"_"$2"_"$3"_"$4"_"$5}}' > init.W
        cat {win_simdir}{wildcards.li}sd.winsim.bed |sed \"s/_rows_/\t/\" | awk '{{OFS="\\t"}}NR>1 {{print $1,$2,$3,$4,$6}} NR==1 {{print $0}}' >{wildcards.li}.sd.pos
        Rscript ./pansdall/winsimprocess.r {wildcards.li}.sd.pos {win_simdir}{wildcards.li}sd.winsim.bed.add
        less -S {win_simdir}{wildcards.li}sd.winsim.bed.add  |awk 'OFS="\\t"{{print $6"#"$8,$9+$3,$9+$4}}' |less -S > wincoorend.tsv
        less -S {win_simdir}{wildcards.li}sd.winsim.bed.add  |awk 'OFS="\\t"{{print $8,$9+$3,$9+$4}}' |less -S > winend_true.tsv
        bedtools intersect -a {minimapindir}{wildcards.li}.other -b winend_true.tsv |less -S |bedtools sort |bedtools merge >winbef.bed
        bedtools makewindows -b winbef.bed -w 1000 -s 500 > ../{output.output_bed}
        """

rule prallel_kmercal:
    input:
        fasta=fadir + "{li}.fa",
        win_bed="{li}/newin"
    output:
        output_bed="{li}/filt/beforegroup.txt" #    threads: 7
    shell:
        """
        cd {wildcards.li}
        rm -r filt || true
        rm *.txt || true
        rm *.pkl || true
        samtools faidx {input.fasta}
        python  ./pansdall/parallel_kmer.py --fa {input.fasta} --i newin --n 20
        python ./pangenome/WGSfea/SDsimnew/writxt.py
        rm *.pkl || true
        find . -maxdepth 1 -name "*.txt" -type f -size 0 -exec rm {{}} \;
        rm -r .snakemake || true
        snakemake -s ./pansdall/SD.smk --config fold={nowdir}{wildcards.li}/ -j 40
        rm *end.a.txt || true
        rm *end.b.txt || true
        """



rule cluster_expand:
    input:
        bef="{li}/filt/beforegroup.txt",
        fasta=fadir + "{li}.fa",
    output:
        outbed="{li}/SPLIT.OK"
    shell:
        """
        cd {wildcards.li}
        cd filt
        cp {nowdir}{input.bef} ./beforegroup.init
        less -S beforegroup.txt  |cut -f 1-3 |less -s |uniq -c |less -S |awk '$1>100' >x
        less -S x |awk 'OFS="\\t"{{print $2,$3,$4}}' >rmin
        bedtools coverage -a rmin -b {SAMRMTRF}/{wildcards.li}/RMTRF/RMTRFSAT.region |awk '$7!=1' >rmin.in
        bash ./pansdall/run_rm_trf.sh rmin.in {input.fasta}  ##这步为跑trf和rm 只是那些低覆盖区域的
        bedtools coverage -a rmin -b <(cat rm_trf.s.bed {SAMRMTRF}/{wildcards.li}/RMTRF/RMTRF.rmtrf)  |awk '$7>=0.9'|cut -f 1-3 >del.pair
        awk 'NR==FNR {{a[$1,$2,$3]; next}} !($1,$2,$3) in a' del.pair beforegroup.txt >beforegroup.filt.txt
        cp beforegroup.filt.txt beforegroup.txt
        cp beforegroup.txt ../
        Rscript ./pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/other/splitgroup.r {nowdir}{wildcards.li}"/filt/split/" {nowdir}{wildcards.li}/
        cd split
        snakemake -s ./pansdall/SD1.smk  --config base_dir={nowdir}{wildcards.li}"/filt/split/" fa={input.fasta} -j 40
        cd ..
        cd ..
        touch  ../{output.outbed}
        """


rule expandt2trmtrf:
    input:
        ha=rules.cluster_expand.output,
        fasta=fadir + "{li}.fa",
    output:
        outbed="{li}/expand.out"
    shell:
        """
        cd {wildcards.li}
        cd filt
        cd split
        paste <( cat {win_simdir}{wildcards.li}sd.winsim.bed.add |cut -f2 |tr ":" "\t" |tr "-" "\t" ) <(cat {win_simdir}{wildcards.li}sd.winsim.bed.add |awk 'OFS="\\t"{{print $8,$9+$3,$9+$4}}') >winsim_t2t.bed
        bash ./pansdall/expandall.sh {t2trmtrf}  {input.fasta}.fai  ./pangenome/WGSfea/sdseq/expand_pansd.r ###bash ./pansdall/expandall.sh ./pangenome/WGSfea/allW/APG3allWnew/scripts/t2trmtrf.cut1M.bed  ./allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/C001-CHA-E01-Mat.fa.fai  ./pangenome/WGSfea/sdseq/expand_pansd.r
        rm rigin|| true
        rm lefin|| true
        cd ..
        rm *pair|| true
        rm *end.txt|| true
        cd ..
        cp ./filt/split/expand.out ../{output.outbed}

        """





###./pansdall/run/no_expand_out_dirs.txt
# input_file="./pansdall/run/no_expand_out_dirs.txt"

# # 遍历每一行
# while IFS= read -r line; do
#     mkdir $line
#     cp ./pansdall/run/$line/expand.out $line/
#     echo "Copied $src_file to current directory."
# done < "$input_file"