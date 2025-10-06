import os
#snakemake -s run.smk --keep-going --keep-incomplete -j 10
# 定义输入和输出文件路径
#inputallres = ".pangenome/WGSfea/allW/reproce.input"
#T2TSD=".PANSDEND/01_REFSD/SD.refine"
T2TSD=".PANSDEND/01_REFSD/para/SD.refine"
CANREGION=".PANSDEND/000_SD_sim/2k/candidate.merge"
sam_list = ".pangenome/WGSfea/allW/APG3allWnew/allsam"
acro=".PANSDEND/DATA/acro.region"
T2Tsegpos=".PANSDEND/DATA/CHM13_APG_HPRC_HGSVC.pos.txt" ##.sim/CHM13_APG_HPRC_HGSVC.pos.txt
script=".PANSDEND/script/"

# 读取 sam.list 文件中的每一行
with open(sam_list) as f:
    lines = [line.strip() for line in f.readlines()]

# 规则定义
rule all:
    input:
        "reproce.input",

rule gene_path:
    input:
        sd={T2TSD},
        can={CANREGION}
    output:
        "reproce.input"
    shell:
        """
        cat {input.sd} |awk 'OFS="\\t"{{print $1,$2,$3"\\n"$4,$5,$6}}' |less -S |sort |uniq |bedtools sort >endendcopy
        cat {input.can} |sort |uniq  |bedtools sort >endendcopy1
        bedtools subtract -a endendcopy -b {acro} -A |grep -v "chrY" >endendcopy.filtacro
        bedtools subtract -a endendcopy1 -b {acro} -A |grep -v "chrY" >endendcopy.filtacro1
        bedtools intersect -a endendcopy.filtacro -b {T2Tsegpos} -wa -wb|awk 'OFS="\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' >sdinter.txt
        bedtools intersect -a endendcopy.filtacro1 -b {T2Tsegpos} -wa -wb|awk 'OFS="\t"{{print $1":"$2"-"$3,$4,$5,$6,$7,$8}}' >sdinter1.txt
        Rscript {script}SDgenepath.r sdinter.txt alldf.txt
        Rscript {script}SDgenepath.r sdinter1.txt alldf1.txt
        cat alldf.txt |tail -n +2 |less -S  >alldfe.txt
        cat alldf1.txt |tail -n +2 |less -S  >alldfe1.txt
        cat alldfe1.txt alldfe.txt >{output}
        """

