import pandas as pd


# singularity_image=config["singularity"]["image"]
# singularity_args=config["singularity"]["args"]
import pandas as pd

# script_dir = config["script_dir"]
# base_dir=config["base_dir"]
# data_dir=config['data_dir']

# script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"
# base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"
# DATA_dir="/home/jmhan/project/APG/github_final/PanSD_end/example/"
script_dir = config["script_dir"]
base_dir=config["base_dir"]
data_dir=config['data_dir']
sam_list=DATA_dir+"samlis.txt"
VCF=DATA_dir+"example_mc.vcf"
GFA=DATA_dir+"example_mc.gfa"
sam_list=DATA_dir+"samlis.txt"
# GFABASE="/home/jmhan/pangenome/software/gfabase"
# RM="/home/Public/software/RepeatMasker-4.1.4/"
# TRF="/home/Public/software/trf-4.09.1/trf"
split_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapsplit")
minimap_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapspnow")
allclus_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/endfilt")
minimap_dir1 = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapspnow")
minimap_dir2 = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/minimapspnow_mer")
allclus_dir = os.path.join(base_dir, "03_refSD/ref_cor/filt/split/endfilt")
allclus_dir1= os.path.join(base_dir, "03_refSD/ref_cor/filt/split/endfilt1")
END= os.path.join(base_dir, "03_refSD/ref_cor/filt/split/END")

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



include: os.path.join(script_dir, "inter_sim/inter_sim.smk")
include: os.path.join(script_dir, "ref_SD/ref_SDcor.smk")
include: os.path.join(script_dir, "ref_SD/SD_cordi/SD.smk")
include: os.path.join(script_dir, "ref_SD/SD_cordi/SD1.smk")
include: os.path.join(script_dir, "ref_SD/SD_cordi/endprocess.smk")
include: os.path.join(script_dir, "ref_SD/sam_sd/sam_SD.smk")
include: os.path.join(script_dir, "nonref_SD/iter_END.smk")



rule all:
    input:
        f"{base_dir}prepro.OK",
        f"{base_dir}03_refSD/ref_cor/para_gene.OK",
        f"{base_dir}03_refSD/ref_cor/chrom_split.OK",
        expand("{base_dir}03_refSD/ref_cor/filt/split/{folder_name}/{folder_name}.log", folder_name=folder_names,base_dir=base_dir),
        os.path.join(base_dir, "03_refSD/ref_cor/filt/split/allclus"),
        expand(os.path.join(split_dir, "part_{i}"), i=["{:03d}".format(i) for i in range(1, 851)]),
        expand(os.path.join(minimap_dir, "mini_{j}"), j=["{:03d}".format(i) for i in range(1, 301)]),
        expand(os.path.join(minimap_dir, "mini_{j}.inter.paf"), j=["{:03d}".format(i) for i in range(1, 301)]),
        expand(os.path.join(minimap_dir, "mini_{j}.other.paf"), j=["{:03d}".format(i) for i in range(1, 301)]),
        expand(os.path.join(minimap_dir, "lar_{m}"), m=["{:03d}".format(i) for i in range(1, 101)]),
        expand(os.path.join(minimap_dir, "lar_{m}.inter.paf"), m=["{:03d}".format(i) for i in range(1, 101)]),
        expand(os.path.join(minimap_dir, "lar_{m}.other.paf"), m=["{:03d}".format(i) for i in range(1, 101)]),
        expand(os.path.join(split_dir, "part_{i}.addre.txt"), i=["{:03d}".format(i) for i in range(1, 851)]),
        expand(os.path.join(split_dir, "part_{i}.minimapbefore"), i=["{:03d}".format(i) for i in range(1, 851)]),
        os.path.join(split_dir, "minimapbefore"),
        os.path.join(minimap_dir, "other.sd"),
        os.path.join(minimap_dir, "1.txt"),
        os.path.join(minimap_dir, "2.txt"),
        expand(os.path.join(minimap_dir, "data1/part1_{d1}"), d1=["{:03d}".format(i) for i in range(1, 501)]),
        expand(os.path.join(minimap_dir, "data2/part2_{d2}"), d2=["{:02d}".format(i) for i in range(1, 51)]),
        expand(os.path.join(minimap_dir, "res/part1_{d1}.txt"), d1=["{:03d}".format(i) for i in range(1, 501)]),
        expand(os.path.join(minimap_dir, "res/part2_{d2}.txt"), d2=["{:02d}".format(i) for i in range(1, 51)]),
        f"{base_dir}/03_refSD/ref_cor/filt/split/check.OK",
        os.path.join(minimap_dir1, "res/refine.bed"),
        os.path.join(allclus_dir, "ourend"),
        os.path.join(allclus_dir1, "SDend"),
        os.path.join(END, "REF.end"),
        f"{base_dir}03_refSD/sam_cor/Wpath/seg_db.OK",
        expand("{base_dir}/03_refSD/sam_cor/Wpath/{sam}.W",base_dir=base_dir, sam=ALLSAM),
        expand("{base_dir}/03_refSD/sam_cor/db/time/{sam}.txt",base_dir=base_dir, sam=ALLSAM),
        f"{base_dir}03_refSD/sam_cor/1_SD_CANPATH/reproce.input",
        f"{base_dir}03_refSD/sam_cor/1_SD_CANPATH/sd.edge",
        expand("{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/time/{line}.txt",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/{line}.sd",base_dir=base_dir, line=lines),
        f"{base_dir}/03_refSD/sam_cor/1_SD_CANPATH/check_process_line.OK",
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/other/{line}.other",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/T2TPANSD/{line}.t2tsd",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/1.OK",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_refSD/sam_cor/2_ALIGN/{line}/RMTRF/2.OK",base_dir=base_dir, line=lines),
        expand("{base_dir}03_nonref/02_ITER/2_SD_GENE/{line}/SPLIT.OK",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/expand.out",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapsplit/part_{i}", i=["{:02d}".format(i) for i in range(1, 51)],line=lines,base_dir=base_dir),
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/part_{i}.paf", i=["{:02d}".format(i) for i in range(1, 51)],line=lines,base_dir=base_dir),
        f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/gene_paf.OK",
        f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/gene_expand.OK",
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/minimapspnow/all1.end.statistics",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/samuniq.pair",base_dir=base_dir, line=lines),
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/otherdfaddsam.txt",base_dir=base_dir, line=lines),
        f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/df.pair",
        f"{base_dir}/03_nonref/02_ITER/2_SD_GENE/allpair.sort.END.e",
        expand("{base_dir}03_nonref/02_ITER/3_PAIR_PATH/output/{sam}can_sd.corregion.bed",base_dir=base_dir, sam=ALLSAM),
        expand("{base_dir}/03_nonref/02_ITER/2_SD_GENE/{line}/processend/{line}.SD.trans",base_dir=base_dir, line=lines)
