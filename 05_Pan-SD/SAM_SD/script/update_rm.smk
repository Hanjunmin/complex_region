
import os
#snakemake -s .PANSDEND/02_SAMSD/02_00-REFSD/update_rm.smk --keep-going --keep-incomplete -j 10


###################################################################
sam_list=".PANSDEND/DATA/allsam"
infold=".PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/"
fold=config['fold']
nowf=config['subfold']  # "t2tsdminimap" #nowf="can"
label=str(config['label'])

with open(sam_list) as f:
    lines = [line.strip() for line in f.readlines()]

lines.remove('CHM13v2')

rule all:
    input:
        expand(infold+"{li}/RMTRF/"+label+".OK", li=lines),
        


rule init_sd_can:
    input:
        rmnow=fold+"{li}/"+nowf+"/rm_trf.s.bed"
    output:
        output_bed=infold+"{li}/RMTRF/"+label+".OK"
    shell:
        """
        cd {infold}{wildcards.li}/RMTRF
        cat <(cat {fold}{wildcards.li}/{nowf}/rmin.in |cut -f 1-3)  RMTRFSAT.region  >RMTRFSAT.region.temp
        cat {fold}{wildcards.li}/{nowf}/rm_trf.s.bed  RMTRF.rmtrf >RMTRF.rmtrf.temp
        cat <(cat RMSAT.sat |awk 'OFS="\\t"{{print $0,"Satellite"}}') <(cat {fold}{wildcards.li}/{nowf}/rm.out |cut -f 1-3,5 )  >RMSAT.sat.temp
        mv RMTRFSAT.region.temp  RMTRFSAT.region
        mv RMTRF.rmtrf.temp  RMTRF.rmtrf
        mv RMSAT.sat.temp  RMSAT.sat
        touch {label}.OK
        """


