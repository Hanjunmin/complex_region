#!/bin/bash


for i in {1..22} X Y; do
    for j in {0..9}; do
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=chr${i}_${j}.rm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --partition=cpu64,cpu128
#SBATCH --output=chr${i}_${j}.rmout
#SBATCH --error=chr${i}_${j}.rmerr
date

mkdir -p ./project/ENDTE/01_INTE/chr${i}
cd  ./project/ENDTE/01_INTE/chr${i}

cat ./project/ENDTE/00_INS/chr${i}/split_${j} |sed "s/seq/\n/g" >split_${j}.fa
samtools faidx split_${j}.fa
RepeatMasker -x -parallel 50 -species human split_${j}.fa -e ncbi
RM2Bed.py split_${j}.fa.out

EOT
done
done

