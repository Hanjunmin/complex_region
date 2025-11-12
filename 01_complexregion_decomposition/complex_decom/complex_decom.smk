import pandas as pd


# singularity_image=config["singularity"]["image"]
# singularity_args=config["singularity"]["args"]
import pandas as pd

script_dir = config["script_dir"]
base_dir=config["base_dir"]
data_dir=config['data_dir']


sam_list=data_dir+"samlis.txt"
samdf = pd.read_csv(sam_list,sep="\t",header=None)
refname=list(samdf[samdf[1]=="reference"][0])[0]
samname=list(samdf[samdf[1]!="reference"][0])
name=list(samdf[0])
lines = name



include: os.path.join(script_dir, "preproc/preprocessing.smk")
include: os.path.join(script_dir, "index_construction/index_construct_anno.smk")

rule all:
    input:
        f"{base_dir}01_preproc/gene_prep.OK",
        f"{base_dir}02_index_construction/gene_indices.OK",
        f"{base_dir}03_model/model.OK"


rule run_model:
    input:
        f"{base_dir}02_index_construction/gene_indices.OK"
    output:
        f"{base_dir}03_model/model.OK"
    shell:
        """
        mkdir -p {base_dir}/03_model/
        cd {base_dir}/03_model/
        python {script_dir}/clustering_model/scripts/model.py --dir {base_dir} --model {script_dir}
        touch {output}
        """
