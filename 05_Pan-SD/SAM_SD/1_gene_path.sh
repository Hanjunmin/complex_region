cd .PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/
#snakemake -s .pangenome/WGSfea/allW/APG3allW/scripts/run.smk ###geneSD and can path 只用前面的gene_path生成一下下一步的reproce.input
snakemake -s .PANSDEND/02_SAMSD/02_00-REFSD/get_p.smk
snakemake -s .iterall/extractpath.smk   --config inputpath=".PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/reproce.input" --keep-going --keep-incomplete -j 30


./project/pansd/zdcopy/allsam
#sample_list="./project/pansd/zdcopy/allsam"
sample_list="./project/PANSDEND/END/1_SD_CANPATH/weno"
cd ./project/PANSDEND/END/1_SD_CANPATH/run
while IFS= read -r line; do
    sample=$(echo "$line" | tr -d '\r\n')
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=${sample}_itermap   # 作业名称
#SBATCH --output=${sample}_itermap.out # 输出文件
#SBATCH --error=${sample}_itermap.err  # 错误文件
#SBATCH --cpus-per-task=15        # 每个任务使用的 CPU 核心数
#SBATCH --partition=cpu64

cd ./project/PANSDEND/END/1_SD_CANPATH/run
/usr/bin/time -v python ./project/PANSDEND/END/1_SD_CANPATH/iterpansd.py --W ./decompose/splinew/${sample}.W --S ${sample} --i ./project/PANSDEND/END/1_SD_CANPATH/reproce.input --o ./output/${sample}can_sd.corregion.bed --p ./project/pansd
EOT
done < "$sample_list"


rsync -av  ./project/PANSDEND/END/1_SD_CANPATH/run/zst clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/zst
cd .PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/
rsync -av   clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/zst ./
cp * .PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/output/

cd ./project/PANSDEND/END/1_SD_CANPATH/run

sample_list=".PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/output/weno"
cd .iterall/transdb/tran
while IFS= read -r line; do
sample=$(echo "$line" | tr -d '\r\n')
cp -r .iterall/transdb/$sample ./ &
done < "$sample_list"

rsync -av  tran clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/tran
rsync -av  clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/tran  ./