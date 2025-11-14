import os
import glob

script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"
base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"
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

# 基础路径 snakemake -s /home/jmhan/pangenome/WGSfea/SD1.smk  -j 2
# base_dir = "/home/jmhan/project/APG/github_final/PanSD/output/reference/filt/split/"
# rpfile =base_dir+"/01_intersim/rm_trf_cut1M.bed"     #"/home/jmhan/project/APG/github_final/PanSD/output/reference/rm_trf_cut1M.bed"
# repeat_trf=base_dir+"/01_intersim/rm_trf.bed"    #"/home/jmhan/project/APG/github_final/PanSD/output/reference/rm_trf.bed"  ##/home/jmhan/HPRC/integrate/sedef/allsedef/{s}/{s1}/Masked/rm_trf.bed
# groupscript=script_dir+"/ref_SD/SD_cordi/scripts/group1_parallel.r"      #"/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/group1_parallel.r"
# runcript=script_dir+"/ref_SD/SD_cordi/scripts/run.r"    #"/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/run.r"
# cigarcall=script_dir+"/ref_SD/SD_cordi/scripts/cigar_call.py" #"/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/cigar_call.py"
# endfiltscript=script_dir+"/ref_SD/SD_cordi/scripts/endfilt.sh"   #"/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/endfilt.sh"
# ednfilt=script_dir+"/ref_SD/SD_cordi/scripts/ednfilt.r"   #"/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/ednfilt.r"
# run2filt=script_dir+"/ref_SD/SD_cordi/scripts/run2filt.r"    #"/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/script/run2filt.r"
repeatmasker="/home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/asm_repeatmasker.out.bed" ## /home/jmhan/PANSDEND/000_SD_sim/2k/T2TSD/asm_repeatmasker.out.bed
# FA="/home/jmhan/project/APG/github_final/PanSD/example/DATA/FA/CHM13v2.fa"
# inteRscipt=script_dir+"/ref_SD/SD_cordi/scripts/interprocess.r"

split_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapsplit")
minimap_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapspnow")
allclus_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/endfilt")
fasta_file1=FA
fasta_file2=FA

all_folders = glob.glob(base_dir + "03_refSD/ref_cor/filt/split/*/")
folder_names = [os.path.basename(f.rstrip('/')) for f in all_folders if os.path.basename(f.rstrip('/')) not in ['minimapsplit', 'minimapspnow','endfilt','minimapspnow_mer','minimapsplit_mer']]


# rule all:
#     input:
#         expand("{base_dir}03_refSD/ref_cor/filt/split/{folder_name}/{folder_name}.log", folder_name=folder_names,base_dir=base_dir),
#         os.path.join(base_dir, "03_refSD/ref_cor/filt/split/allclus"),
#         expand(os.path.join(split_dir, "part_{i}"), i=["{:03d}".format(i) for i in range(1, 851)]),
#         expand(os.path.join(minimap_dir, "mini_{j}"), j=["{:03d}".format(i) for i in range(1, 301)]),
#         expand(os.path.join(minimap_dir, "mini_{j}.inter.paf"), j=["{:03d}".format(i) for i in range(1, 301)]),
#         expand(os.path.join(minimap_dir, "mini_{j}.other.paf"), j=["{:03d}".format(i) for i in range(1, 301)]),
#         expand(os.path.join(minimap_dir, "lar_{m}"), m=["{:03d}".format(i) for i in range(1, 101)]),
#         expand(os.path.join(minimap_dir, "lar_{m}.inter.paf"), m=["{:03d}".format(i) for i in range(1, 101)]),
#         expand(os.path.join(minimap_dir, "lar_{m}.other.paf"), m=["{:03d}".format(i) for i in range(1, 101)]),
#         expand(os.path.join(split_dir, "part_{i}.addre.txt"), i=["{:03d}".format(i) for i in range(1, 851)]),
#         expand(os.path.join(split_dir, "part_{i}.minimapbefore"), i=["{:03d}".format(i) for i in range(1, 851)]),
#         os.path.join(split_dir, "minimapbefore"),
#         os.path.join(minimap_dir, "other.sd"),
#         os.path.join(minimap_dir, "1.txt"),
#         os.path.join(minimap_dir, "2.txt"),
#         expand(os.path.join(minimap_dir, "data1/part1_{d1}"), d1=["{:03d}".format(i) for i in range(1, 501)]),
#         expand(os.path.join(minimap_dir, "data2/part2_{d2}"), d2=["{:02d}".format(i) for i in range(1, 51)]),
#         expand(os.path.join(minimap_dir, "res/part1_{d1}.txt"), d1=["{:03d}".format(i) for i in range(1, 501)]),
#         expand(os.path.join(minimap_dir, "res/part2_{d2}.txt"), d2=["{:02d}".format(i) for i in range(1, 51)]),
#         f"{base_dir}/03_refSD/ref_cor/filt/split/check.OK"
#         # os.path.join(allclus_dir, "ourend")
        

#         # expand("{folder_name}/{folder_name}_cluster.bed", folder_name=[os.path.basename(f.rstrip('/')) for f in all_folders]) +
#         # expand("{folder_name}/{folder_name}_largecluster.txt", folder_name=[os.path.basename(f.rstrip('/')) for f in all_folders]) +
#         # expand("{folder_name}/{folder_name}_large.txt", folder_name=[os.path.basename(f.rstrip('/')) for f in all_folders])




rule run_rscript:
    input:
        init=f"{base_dir}prepro.OK",
        bef=f"{base_dir}03_refSD/ref_cor/chrom_split.OK",
        txt="{base_dir}/03_refSD/ref_cor/filt/split/{folder_name}/{folder_name}.txt",
        clu1="{base_dir}/03_refSD/ref_cor/filt/split/{folder_name}/{folder_name}.clu1"
    output:
        log="{base_dir}/03_refSD/ref_cor/filt/split/{folder_name}/{folder_name}.log",
        # clu2="{folder_name}/{folder_name}_large.txt",
        # clu3="{folder_name}/{folder_name}_largecluster.txt",
        # clu4="{folder_name}/{folder_name}_cluster.bed"
    shell:
        """
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/group1_parallel.r {input.txt} {input.clu1} {wildcards.folder_name} {base_dir}/03_refSD/ref_cor/filt/split/ > {output.log}
        """




rule combine_clusters:
    input:
        expand("{base_dir}/03_refSD/ref_cor/filt/split/{folder_name}/{folder_name}.log", folder_name=folder_names, base_dir=base_dir)
    output:
        allclus = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/allclus")
    shell:
        """
        cd {base_dir}/03_refSD/ref_cor/filt/split/
        find {base_dir}/03_refSD/ref_cor/filt/split/ -name "*_largecluster.txt" -exec cat {{}} + > largecluster.bed || touch largecluster.bed
        find {base_dir}/03_refSD/ref_cor/filt/split/ -name "*_cluster.bed" -exec cat {{}} + > cluster.txt || touch cluster.txt
        cat largecluster.bed cluster.txt > {output.allclus}
        """


rule splitcombine_clusters:
    input:
        rules.combine_clusters.output
    output:
        expand(os.path.join(split_dir, "part_{i}"), i=["{:03d}".format(i) for i in range(1, 851)])
    shell:
        """
            split -n l/850 {input} --numeric-suffixes=1 --suffix-length=0 {split_dir}/part_
        """

rule expand_thereshold:
    input:
        allclus = os.path.join(split_dir, "part_{i}")
    output:
        addre = os.path.join(split_dir, "part_{i}.addre.txt"),
        minimapbefore = os.path.join(split_dir, "part_{i}.minimapbefore")
    shell:
        """
        if [ ! -s {input.allclus} ]; then
            echo "输入文件 {input.allclus} 为空，跳过处理"
            touch {output.addre}
            touch {output.minimapbefore}
            exit 0
        fi
        cat <(nl {input.allclus}) | awk 'OFS="\\t"{{print $2,$3,$4,$5,$6,$7,$1}}' > {split_dir}/part_{wildcards.i}.left.end
        cat <(nl {input.allclus}) | awk 'OFS="\\t"{{print $5,$6,$7,$1}}' > {split_dir}/part_{wildcards.i}.right.end

        bedtools intersect -a {split_dir}/part_{wildcards.i}.left.end -b "{base_dir}/01_intersim/rm_trf_cut1M.bed" -v > {split_dir}/part_{wildcards.i}.leftno.inter
        bedtools intersect -a {split_dir}/part_{wildcards.i}.left.end -b "{base_dir}/01_intersim/rm_trf_cut1M.bed" -wa -wb > {split_dir}/part_{wildcards.i}.left_iner_rmtrf.txt
        bedtools intersect -a {split_dir}/part_{wildcards.i}.right.end -b "{base_dir}/01_intersim/rm_trf_cut1M.bed" -v > {split_dir}/part_{wildcards.i}.rightno.inter
        bedtools intersect -a {split_dir}/part_{wildcards.i}.right.end -b "{base_dir}/01_intersim/rm_trf_cut1M.bed" -wa -wb > {split_dir}/part_{wildcards.i}.right_iner_rmtrf.txt

        Rscript {script_dir}/ref_SD/SD_cordi/scripts/run.r {split_dir}/part_{wildcards.i}.left_iner_rmtrf.txt {split_dir}/part_{wildcards.i}.right_iner_rmtrf.txt {output.addre} {split_dir}/part_{wildcards.i}.leftno.inter {split_dir}/part_{wildcards.i}.rightno.inter

        less -S {output.addre} | cut -f2-7 > {output.minimapbefore}
        """

rule combine_thereshold:
    input:
        expand(os.path.join(split_dir, "part_{i}.minimapbefore"), i=["{:03d}".format(i) for i in range(1, 851)])
    output:
        os.path.join(split_dir, "minimapbefore"),
    shell:
        """
        find {split_dir} -name "part_*minimapbefore" -exec cat {{}} + > {output}
        """


rule splitminimap:
    input:
        rules.combine_thereshold.output,
        f"{base_dir}01_intersim/rm_trf.bed"
    output:
        expand(os.path.join(minimap_dir, "lar_{m}"), m=["{:03d}".format(i) for i in range(1, 101)]),
        expand(os.path.join(minimap_dir, "mini_{j}"), j=["{:03d}".format(j) for j in range(1, 301)])
    shell:
        """
            cd {minimap_dir}
            cat {input}  |awk '$3-$2>=500000' >{minimap_dir}/large
            bedtools coverage -a {minimap_dir}/large -b {base_dir}/01_intersim/rm_trf.bed   >y 
            less -S y |sort -k10,10n |less -S |awk '{{print $0,$9-$8}}' |less -S  |less -S |sort -k2,2n |less -S |awk '!($1==$4 && $2<=$6 && $3>=$5 && $10>=0.98)'|cut -f 1-6 >{minimap_dir}/newlarge 
            cat {input}  |awk '$3-$2<500000' >{minimap_dir}/sm
            bedtools coverage -a {minimap_dir}/sm -b {base_dir}/01_intersim/rm_trf.bed  >y 
            less -S y |sort -k10,10n |less -S |awk '{{print $0,$9-$8}}' |less -S  |less -S |sort -k2,2n |less -S |awk '!($1==$4 && $2<=$6 && $3>=$5 && $10>=0.98)'|cut -f 1-6 >{minimap_dir}/newsm
            split -n l/100 {minimap_dir}/newlarge --numeric-suffixes=1 --suffix-length=0 {minimap_dir}/lar_
            split -n l/300 {minimap_dir}/newsm --numeric-suffixes=1 --suffix-length=0 {minimap_dir}/mini_
        """
rule minimap2:
    input:
        allmini = os.path.join(minimap_dir, "lar_{m}"),
        fa1=fasta_file1,
        fa2=fasta_file2
    output:
        paf1 = os.path.join(minimap_dir, "lar_{m}.inter.paf"),
        paf2 = os.path.join(minimap_dir, "lar_{m}.other.paf"),
    shell:
        """

        prefix=$(basename {input.allmini})
        less -S {input.allmini} |awk '$1==$4 && $2<=$6 && $3>=$5' >{minimap_dir}/$prefix.inter ||true
        if [ -s "{minimap_dir}/$prefix.inter" ]; then
            Rscript {script_dir}/ref_SD/SD_cordi/scripts/interprocess.r {minimap_dir}/$prefix.inter
            less -S {minimap_dir}/$prefix.inter |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' |sort |uniq > {minimap_dir}/$prefix.inter.bed
            grep -vFf {minimap_dir}/$prefix.inter {input.allmini} |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' |sort |uniq > {minimap_dir}/$prefix.other.bed
        else
            cat {input.allmini} |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' |sort |uniq > {minimap_dir}/$prefix.other.bed
        fi
        
        
        less -S {minimap_dir}/$prefix.other.bed |cut -f 1-3 |sort |uniq >{minimap_dir}/$prefix.left
        if [ ! -s "{minimap_dir}/$prefix.other.bed" ]; then
            touch {output.paf2}
        else
            while IFS=$'\\t' read -r col1 col2 col3
            do
                awk -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' {minimap_dir}/$prefix.other.bed > {minimap_dir}/$prefix.right
                minimap2 -c --eqx <(samtools faidx {input.fa2} "$col1:$col2-$col3") <(samtools faidx {input.fa1} $(cat {minimap_dir}/$prefix.right |awk '{{print $4":"$5"-"$6}}'))   -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >> {output.paf2}
            done < {minimap_dir}/$prefix.left
        fi

        less -S {minimap_dir}/$prefix.inter.bed |cut -f 1-3 |sort |uniq >{minimap_dir}/$prefix.left
        if [ ! -s "{minimap_dir}/$prefix.inter.bed" ]; then
            touch {output.paf1}
        else
            while IFS=$'\\t' read -r col1 col2 col3
            do
            awk -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' {minimap_dir}/$prefix.inter.bed > {minimap_dir}/$prefix.right
            minimap2 -c --eqx <(samtools faidx {input.fa2} "$col1:$col2-$col3")  <(samtools faidx {input.fa1} $(cat {minimap_dir}/$prefix.right |awk '{{print $4":"$5"-"$6}}'))  -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -r 1000,1000 -P >> {output.paf1}
            done < {minimap_dir}/$prefix.left
        fi

       
        """

rule minimap:
    input:
        allmini = os.path.join(minimap_dir, "mini_{j}"),
        fa1=fasta_file1,
        fa2=fasta_file2
    output:
        paf1 = os.path.join(minimap_dir, "mini_{j}.inter.paf"),
        paf2 = os.path.join(minimap_dir, "mini_{j}.other.paf"),
    shell:
        """
        prefix=$(basename {input.allmini})
        less -S {input.allmini} |awk '$1==$4 && $2<=$6 && $3>=$5' >{minimap_dir}/$prefix.inter ||true
        echo "C"
        if [ -s "{minimap_dir}/$prefix.inter" ]; then
            Rscript {script_dir}/ref_SD/SD_cordi/scripts/interprocess.r {minimap_dir}/$prefix.inter
            echo "X"
            less -S {minimap_dir}/$prefix.inter |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' |sort |uniq > {minimap_dir}/$prefix.inter.bed ||true
            grep -vFf {minimap_dir}/$prefix.inter {input.allmini} |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' |sort |uniq > {minimap_dir}/$prefix.other.bed ||true
            echo "Y"
        else
            cat {input.allmini} |sort -k1,1 -k2,2n -k3,3n  | awk 'OFS="\\t" {{print $1,$2+1,$3,$4,$5+1,$6}}' |sort |uniq > {minimap_dir}/$prefix.other.bed ||true
        fi
        
        less -S {minimap_dir}/$prefix.other.bed |cut -f 1-3 |sort |uniq >{minimap_dir}/$prefix.left
        echo "A"

        if [ ! -s "{minimap_dir}/$prefix.other.bed" ]; then
            touch {output.paf2}
        else
            while IFS=$'\\t' read -r col1 col2 col3
            do
                awk -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' {minimap_dir}/$prefix.other.bed > {minimap_dir}/$prefix.right
                minimap2 -c --eqx  <(samtools faidx {input.fa2} "$col1:$col2-$col3") <(samtools faidx {input.fa1} $(cat {minimap_dir}/$prefix.right |awk '{{print $4":"$5"-"$6}}'))  -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 >> {output.paf2}
            done < {minimap_dir}/$prefix.left
        fi
        echo "B"
        less -S {minimap_dir}/$prefix.inter.bed |cut -f 1-3 |sort |uniq >{minimap_dir}/$prefix.left
        if [ ! -s "{minimap_dir}/$prefix.inter.bed" ]; then
            touch {output.paf1}
        else
        while IFS=$'\\t' read -r col1 col2 col3
        do
            awk -v col1="$col1" -v col2="$col2" -v col3="$col3" '$1==col1 && $2==col2 && $3==col3' {minimap_dir}/$prefix.inter.bed > {minimap_dir}/$prefix.right
            minimap2 -c --eqx  <(samtools faidx {input.fa2} "$col1:$col2-$col3")  <(samtools faidx {input.fa1} $(cat {minimap_dir}/$prefix.right |awk '{{print $4":"$5"-"$6}}')) -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -r 1000,1000 -P >> {output.paf1}
        done < {minimap_dir}/$prefix.left
        fi
        """




rule combine_minimap_REF:
    input:
        expand(os.path.join(minimap_dir, "lar_{m}.inter.paf"), m=["{:03d}".format(i) for i in range(1, 101)]),
        expand(os.path.join(minimap_dir, "mini_{j}.inter.paf"), j=["{:03d}".format(j) for j in range(1, 301)])
    output:
        os.path.join(minimap_dir, "other.sd"),
    threads: 1 
    shell:
        """
        cd {base_dir}/03_refSD/ref_cor/filt/split
        find {minimap_dir} -name "*.paf" -exec cat {{}} + > {minimap_dir}/repea_mini
        python {script_dir}/ref_SD/SD_cordi/scripts/cigar_call.py  --paf  {minimap_dir}/repea_mini  --o  {minimap_dir}/all1.end.statistics
        bash {script_dir}/ref_SD/SD_cordi/scripts/endfilt.sh {repeatmasker} {base_dir}/01_intersim/rm_trf.bed {minimap_dir}/all1.end.statistics {output}

        """

rule process_combine_minimap:
    input:
        expand(os.path.join(minimap_dir, "lar_{m}.inter.paf"), m=["{:03d}".format(i) for i in range(1, 101)]),
        expand(os.path.join(minimap_dir, "mini_{j}.inter.paf"), j=["{:03d}".format(j) for j in range(1, 301)])
    output:
        os.path.join(minimap_dir, "1.txt"),
        os.path.join(minimap_dir, "2.txt")
    threads: 1 
    shell:
        """
        find {minimap_dir} -name "*.other.paf" -exec cat {{}} + > {minimap_dir}/otherrepea_mini
        python {script_dir}/ref_SD/SD_cordi/scripts/cigar_call.py  --paf  {minimap_dir}/otherrepea_mini  --o  {minimap_dir}/all1.end.statisticsn
        awk '$11>1000 && $10>=0.85 && $10<0.9 && $11<100000 {{print $0}}' {minimap_dir}/all1.end.statisticsn > {output[0]}
        awk '$11<1000 && $10>=0.8 && $11>800 {{print $0}}' {minimap_dir}/all1.end.statisticsn  > {output[1]}
       
        """

rule split_files1:
    input:
        os.path.join(minimap_dir, "1.txt")
    output:
        expand("{minimap_dir}/data1/part1_{d1}", d1=["{:03d}".format(i) for i in range(1, 501)],minimap_dir=minimap_dir)
    shell:
        """
        cd {minimap_dir}
        split -n l/500 {input} --numeric-suffixes=1 --suffix-length=0 data1/part1_
        cd ..
        """

rule split_files2:
    input:
        os.path.join(minimap_dir, "2.txt")
    output:
        expand("{minimap_dir}/data2/part2_{d2}", d2=["{:02d}".format(i) for i in range(1, 51)],minimap_dir=minimap_dir)
    shell:
        """
        cd {minimap_dir}
        split -n l/50 {input} --numeric-suffixes=1 --suffix-length=0 data2/part2_
        cd ..
        """

rule run_refine1:
    input:
        "{minimap_dir}/data1/part1_{d1}"
    output:
        "{minimap_dir}/res/part1_{d1}.txt"
    shell:
        "python {script_dir}/ref_SD/SD_cordi/scripts/refine1.py  --i {input} --o {output} --fa {FA}"

rule run_refine2:
    input:
        "{minimap_dir}/data2/part2_{d2}"
    output:
        "{minimap_dir}/res/part2_{d2}.txt"
    shell:
        "python {script_dir}/ref_SD/SD_cordi/scripts/refine1.py   --i {input} --o {output} --fa {FA}"


rule check_al_gene:
    input:
        expand("{minimap_dir}/res/part2_{d2}.txt", d2=["{:02d}".format(i) for i in range(1, 51)],minimap_dir=minimap_dir),
        expand("{minimap_dir}/res/part1_{d1}.txt", d1=["{:03d}".format(i) for i in range(1, 501)],minimap_dir=minimap_dir)
    output:
        f"{base_dir}/03_refSD/ref_cor/filt/split/check.OK"
    shell:
        '''
        cd {base_dir}/03_refSD/ref_cor/filt/split
        touch {output}
        '''




