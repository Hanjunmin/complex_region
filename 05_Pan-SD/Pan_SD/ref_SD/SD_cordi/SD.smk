import glob
import os
import pickle
#snakemake -s /home/jmhan/pangenome/WGSfea/SD.smk  -j 70
# Define the directory containing your files

# script_dir="/home/jmhan/project/APG/github_final/PanSD_end/Pan_SD/"
# base_dir="/home/jmhan/project/APG/github_final/PanSD_end/run/"

samples=[os.path.splitext(os.path.basename(f))[0] for f in glob.glob("./03_refSD/ref_cor/*end.pkl")]
a_txt_files = [os.path.join(base_dir+"03_refSD/ref_cor/", f"{sample}.a.txt") for sample in samples]


# Extract sample names (base names without suffix)


# rule all:
#     input:
#         f"{base_dir}03_refSD/ref_cor/filt/beforegroup.txt",
#         f"{base_dir}03_refSD/ref_cor/chrom_split.OK"


rule process_file:
    input:
        bef=f"{base_dir}03_refSD/ref_cor/para_gene.OK",
        inpkl=lambda wildcards: f"{base_dir}/03_refSD/ref_cor/{wildcards.sample}.pkl",
    output:
        a_txt="{base_dir}/03_refSD/ref_cor/{sample}.a.txt",
        b_txt="{base_dir}/03_refSD/ref_cor/{sample}.b.txt",
        c_txt="{base_dir}/03_refSD/ref_cor/{sample}.c.txt"
    run:
        with open(input.inpkl, 'rb') as f:
            loaded_vars = pickle.load(f)
        pairset = loaded_vars['var2']
        jacset = loaded_vars['var3']
        otherid = loaded_vars['var1']
        loaded_vars.clear()
        with open(output.a_txt, 'w') as file:
            for key, value in pairset.items():
                for val in value:
                    file.write(f"{key}\t{val}\n")
        with open(output.b_txt, 'w') as file:
            for key, value in jacset.items():
                for val in value:
                    file.write(f"{key}\t{val}\n")
        with open(output.c_txt, 'w') as file:
            for value in otherid:
                file.write(f"{value}\n")


rule filter_pairs:
    input:
        a_txt=rules.process_file.output.a_txt,
        b_txt=rules.process_file.output.b_txt
    output:
        "{base_dir}03_refSD/ref_cor/filt/{sample}.end.txt",
    shell:
        """
        echo "Processing file: {input.a_txt} and {input.b_txt}"
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/filt_pair.r {input.a_txt} {input.b_txt} {output}
        """

rule bash_pairs:
    input:
         filtered_file=rules.filter_pairs.output
    output:
        "{base_dir}03_refSD/ref_cor/filt/{sample}.filt.pair"
    shell:
        """
        cat {input.filtered_file} | less -S | tr ':' '\t' | tr '-' '\t' | less | awk '((!($1 == $4 && $5 > $2 && $5 < $3) && $10 >= 0.05) || ($1 == $4 && $5 > $2 && $5 < $3 && $10 > 0.5)){{print $0}}' > {output}
        """

rule process_individual_pairs:
    input:
        file=rules.bash_pairs.output
    output:
        "{base_dir}/03_refSD/ref_cor/filt/{sample}.processed.pair"
    group: "processing_group"
    shell:
        """
        cat {input.file} | less -S | tr ':' '\t' | tr '-' '\t' | awk 'OFS="\\t"{{print $1,$2,$3,$4,$5,$6}}' > {output}
        """

rule combine_processed_pairs:
    input:
        f"{base_dir}03_refSD/ref_cor/para_gene.OK",
        expand("{base_dir}/03_refSD/ref_cor/filt/{sample}.processed.pair",base_dir=base_dir,sample=samples)
    output:
        f"{base_dir}03_refSD/ref_cor/filt/beforegroup.txt"
    shell:
        """
        echo {input}
        echo {threads}
        cat {input} > {output}
    """


rule group:
    input:
        rules.combine_processed_pairs.output
    output:
        f"{base_dir}03_refSD/ref_cor/chrom_split.OK"
    shell:
        '''
        cd {base_dir}/03_refSD/ref_cor
        cp ./filt/beforegroup.txt ./
        Rscript {script_dir}/ref_SD/SD_cordi/scripts/1.group.r {base_dir}03_refSD/ref_cor/filt/split/ {base_dir}03_refSD/ref_cor/
        touch chrom_split.OK
        '''

## mkdir split
##Rscript /home/jmhan/pangenome/WGSfea/allW/1.group.r $fold"/filt/split/" $fold
