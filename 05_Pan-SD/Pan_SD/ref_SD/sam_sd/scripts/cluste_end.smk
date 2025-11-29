import os
import glob

# 需要在split下运行 snakemake -s /home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/cluste_end.smk
script_dir =config["script"] #"/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/Pan-SD/pansam/datafrommao/pansdscript/"
#cigarcall="/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/Pan-SD/pansam/datafrommao/pansdscripts/paf2bed.py"
ednfilt=script_dir+"ednfilt.r"
ednfilt2=script_dir+"endfilt2.r"
run2filt=script_dir+"run2filt.r"
#fasta_file1="/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/t2tsdminimap/HG00438_Pat.fa"
fasta_file1=config["fasta_file1"]
rule all:
    input:
        "SDend",
        "ourend"
        



rule merge_end:
    input:
        "filin.bedend"
    output:
        "ourend"
    threads: 1  
    shell:
        """
        mkdir -p endfilt && cd endfilt
        rm -f *
        Rscript {ednfilt} ../{input}  {fasta_file1} {script_dir}"paf2bed.py"
        ls |grep "csv"|less -S |tr "_" "\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
        Rscript {run2filt}
        touch nomatch.csv
        touch beforeminus0.9.csv
        cat *beforeminus0.9.csv *nomatch.csv *after.csv *nointe.csv  | tr "," "\t" | tail -n  +2 |awk -F'\t' '$10>=0.9 && $11>=1000 && $12>=0.5{{print $0}}'|less -S |grep -v "V1" |sort |uniq >ourend
        cd ..
        cp ./endfilt/ourend ./{output}
        """


rule merge_end2:
    input:
        "ourend"
    output:
        "SDend"
    threads: 1  
    shell:
        """
        mkdir -p endfilt1 && cd endfilt1
        rm -f *
        less -S  ../{input} |cut -f 1-12 >panSD_T2T.bed
        Rscript {ednfilt2} panSD_T2T.bed  {fasta_file1} {script_dir}"paf2bed.py"
        ls |grep "csv"|less -S |tr "_" "\t" |less -S |tail -n +3 |awk 'OFS="_"{{print $1,$2}}' |sort |uniq |less -S >chlist
        Rscript {run2filt}
        touch nomatch.csv
        touch beforeminus0.9.csv
        cat *beforeminus0.9.csv *nomatch.csv *after.csv *nointe.csv  | tr "," "\t" | tail -n  +2 |awk -F'\t' '$10>=0.9 && $11>=1000 && $12>=0.5{{print $0}}'|less -S |grep -v "V1"|cut -f 1-12 |sort |uniq >ourend
        cd ..
        cp ./endfilt1/ourend ./{output}
        """



