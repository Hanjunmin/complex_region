fold=".PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/"
subfold="t2tsdminimap"
label=1
snakemake -s .PANSDEND/02_SAMSD/02_00-REFSD/update_rm.smk --config fold=$fold subfold=$subfold label=$label  -j 50


fold=".PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/"
subfold="can"
label=2
snakemake -s .PANSDEND/02_SAMSD/02_00-REFSD/update_rm.smk --config fold=$fold subfold=$subfold label=$label  -j 50
