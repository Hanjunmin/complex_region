import pandas as pd


# singularity_image=config["singularity"]["image"]
# singularity_args=config["singularity"]["args"]

script_dir = config["script_dir"]
print(script_dir)
base_dir=config["base_dir"]
DATA_dir=config['data_dir']
sam_list=DATA_dir+"samlis.txt"
VCF=DATA_dir+"example_mc.vcf"
GFA=DATA_dir+"example_mc.gfa"
sam_list=DATA_dir+"samlis.txt"
split_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/minimapsplit")
minimap_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/minimapspnow")
allclus_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/endfilt")
minimap_dir1 = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/minimapspnow")
minimap_dir2 = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/minimapspnow_mer")
allclus_dir = os.path.join(base_dir, "02_refSD/ref_cor/filt/split/endfilt")
allclus_dir1= os.path.join(base_dir, "02_refSD/ref_cor/filt/split/endfilt1")
END= os.path.join(base_dir, "02_refSD/ref_cor/filt/split/END")

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

print(os.path.join(script_dir, "ref_SD/sam_sd/sam_SD.smk"))

include: os.path.join(script_dir, "inter_sim/inter_sim.smk")
include: os.path.join(script_dir, "ref_SD/ref_SDcor.smk")
include: os.path.join(script_dir, "ref_SD/SD_cordi/endprocess.smk")
include: os.path.join(script_dir, "ref_SD/sam_sd/sam_SD.smk")
include: os.path.join(script_dir, "nonref_SD/iter_END.smk")
include: os.path.join(script_dir, "SD_SUM/ref_SDsum.smk")


rule all:
    input:
        f"{base_dir}prepro.OK",
        f"{base_dir}02_refSD/ref_cor/para_gene.OK",
        f"{base_dir}/02_refSD/ref_cor/filt/split/check.OK",
        os.path.join(END, "REF.end"),
        f"{base_dir}/02_refSD/sam_cor/2_ALIGN/SAM_SD.OK",
        f"{base_dir}/03_nonref/non_ref.OK",
        f"{base_dir}ref.SD.al",
        f"{base_dir}iter.SD.al"
