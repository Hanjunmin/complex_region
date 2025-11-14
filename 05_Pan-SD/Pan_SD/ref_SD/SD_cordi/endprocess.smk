import os
import glob

# 需要在split下运行 snakemake -s /dssg/home/acct-clsmyf/clsmyf-user1/data/SD/endprocess.smk
base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"
script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"

# repeat_trf="/home/jmhan/project/APG/github_final/PanSD/output/reference/rm_trf.bed"  ##/home/jmhan/HPRC/integrate/sedef/allsedef/{s}/{s1}/Masked/rm_trf.bed
# cigarcall="/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/cigar_call.py"
# endfiltscript="/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/endfilt.sh"
# ednfilt={script_dir}/ref_SD/SD_cordi/scripts/ednfilt.r "/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/ednfilt.r"
# run2filt={script_dir}/ref_SD/SD_cordi/scripts/run2filt.r "/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/run2filt.r"
# ednfilt2={script_dir}/ref_SD/SD_cordi/scripts/endfilt2.r "/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/endfilt2.r"
DATA_dir="/home/jmhan/project/APG/github_final/PanSD_end/example/"


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
repeatmasker="/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/asm_repeatmasker.out.bed" ## /home/jmhan/HPRC/integrate/sedef/allsedef/{s}/{s1}/Masked/asm_repeatmasker.out.bed

minimap_dir1 = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapspnow")
minimap_dir2 = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapspnow_mer")
allclus_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/endfilt")
allclus_dir1= os.path.join(base_dir, "03_refSD/ref_cor/filt/split/endfilt1")
END= os.path.join(base_dir, "03_refSD/ref_cor/filt/split/END")
fasta_file1=FA
fasta_file2=FA


# rule all:
#     input:
#         os.path.join(minimap_dir1, "res/refine.bed"),
#         os.path.join(allclus_dir, "ourend"),
#         os.path.join(allclus_dir1, "SDend"),
#         os.path.join(END, "REF.end"),
        


rule refine_combine:
    input:
        f"{base_dir}/03_refSD/ref_cor/filt/split/check.OK",
        os.path.join(minimap_dir1, "other.sd")
    output:
        os.path.join(minimap_dir1, "res/refine.bed")
    shell:
        """
        cd {minimap_dir1}
        cat  {minimap_dir1}/res/*.txt | awk 'OFS="\\t"{{print $4,$18,$19,$1,$20,$21,$16,$17,$22,$13,$14,$15}}' >{minimap_dir1}/res/refine || touch {minimap_dir1}/res/refine
        bash  {script_dir}/ref_SD/SD_cordi/scripts/endfilt.sh {repeatmasker} {base_dir}/01_intersim/rm_trf.bed {minimap_dir1}/res/refine {output}

        """

rule merge_end:
    input:
        refineres=os.path.join(minimap_dir1, "res/refine.bed"),
        minires=os.path.join(minimap_dir1, "other.sd"),
    output:
        os.path.join(allclus_dir, "ourend")
    threads: 1  
    shell:
        """
        mkdir -p {base_dir}/03_refSD/ref_cor/filt/split/endfilt && cd {base_dir}/03_refSD/ref_cor/filt/split/endfilt
        cat {input.refineres} {input.minires} >other.sd
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/ednfilt.r other.sd  {fasta_file1} {script_dir}/ref_SD/SD_cordi/scripts/cigar_call.py
        ls |grep "csv"|less -S |tr "_" "\\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/run2filt.r
        echo "A"
        shopt -s nullglob
        cat *beforeminus0.9.csv *nomatch.csv *after.csv *nointe.csv  >reg.temp 
        shopt -u nullglob
        echo "D"
        cat  reg.temp | tr "," "\\t" | tail -n  +2 |awk -F'\\t' '$10>=0.9 && $11>=1000 && $12>=0.5{{print $0}}'|less -S |grep -v "V1" |sort |uniq >{output}
        echo "B"
        cd ..
        """


rule merge_end2:
    input:
        os.path.join(allclus_dir, "ourend")
    output:
        os.path.join(allclus_dir1, "SDend")
    threads: 1  
    shell:
        """
        mkdir -p {base_dir}/03_refSD/ref_cor/filt/split/endfilt1 && cd {base_dir}/03_refSD/ref_cor/filt/split/endfilt1
        less -S  {input} |cut -f 1-12 >panSD_T2T.bed
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/endfilt2.r panSD_T2T.bed  {fasta_file1} {script_dir}/ref_SD/SD_cordi/scripts/cigar_call.py
        ls |grep "csv"|less -S |tr "_" "\\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/run2filt.r
        shopt -s nullglob
        cat *beforeminus0.9.csv *nomatch.csv *after.csv *nointe.csv  | tr "," "\\t" | tail -n  +2 |awk -F'\\t' '$10>=0.9 && $11>=1000 && $12>=0.5{{print $0}}'|less -S |grep -v "V1" |sort |uniq >{output}
        shopt -u nullglob
        """



rule GENE_ENDSD:
    input:
        os.path.join(allclus_dir1, "SDend")
    output:
        os.path.join(END, "REF.end")
    threads: 1  
    shell:
        """
        cd {END}
        cat {input} |cut -f 1-12 |grep -v "chrY" >T2Tchm13v1.1SD_delchrY.txt
        cat {repeatmasker} | grep Satellite | cut -f 1-3 | bedtools sort -i - |  bedtools merge -i - | bedtools coverage -header -a <(cat T2Tchm13v1.1SD_delchrY.txt |tail -n +2) -b -  >end1.statistics.ascov.txt
        less -S  end1.statistics.ascov.txt |awk '$16<=0.7{{print $0}}' |cut -f 1-12  >end1.statistics.ascov.end
        bedtools coverage -a <(less -S end1.statistics.ascov.end) -b  {base_dir}/01_intersim/rm_trf.bed|less -S |awk '{{print $0,$15-$14}}' |awk '$17>=500{{print $0}}' |cut -f 1-12 >SDend1.txt
        less -S SDend1.txt |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' |less -S >SDend1.txt.lef
        bedtools coverage -a <(less -S SDend1.txt.lef) -b  {base_dir}/01_intersim/rm_trf.bed|less -S |awk '{{print $0,$15-$14}}' |awk '$17>=500{{print $0}}'   >SDend1.txt
        less -S SDend1.txt | awk '!($1==$4 && (($5<=$2 && $6>=$2)||($5<=$3 && $6>=$3)||($2<=$5 && $3>=$5)||($2<=$6 && $3>=$6)))' |sort |uniq |awk 'OFS="\\t"{{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}}' >{output}
        """
