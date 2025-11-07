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
        f"{base_dir}01_preproc/cor_window.bed",
        f"{base_dir}01_preproc/Wpath/gene_W.OK",
        f"{base_dir}01_preproc/reg.posCHM.txt",
        expand("{base_dir}01_preproc/conti_path/{li}_REFcor.pos",base_dir=base_dir,li=lines),
        f"{base_dir}01_preproc/conti_path/reg.vcf_corr_dot.en",
        f"{base_dir}01_preproc/conti_path/process_vcf.OK",
        f"{base_dir}01_preproc/gene_db.OK",


rule run_model:
    input:
        os.path.join(base_dir, "01_preproc/gene_db.OK"),  
    output:
        f"{base_dir}/03_model/model.OK"
    shell:
        """
        mkdir -p {base_dir}/03_model/
        cd {base_dir}/03_model/
        python {model_script} --dir {base_dir} --model {script_dir}
        touch {output}
        """
