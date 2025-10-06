#!/bin/bash

# s=59
# e=$(echo "$s+20"|bc)
s=0
e=101
for i in {1..22} X Y; do
    sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=chr$i.${s}.VNTR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=cpu64,cpu128
#SBATCH --output=chr${i}_${s}.out
#SBATCH --error=chr${i}_${s}.err
date
cd ./project/ENDTE/00_INS/
mkdir chr$i
cd chr$i
#mkdir chr$i"_"$s
#cd chr$i"_"$s
#python ./project/ENDTE/graph_VNTRend.py --CHR "chr$i" --T2Tpos "./project/panVNTR/data/chr$i.posCHM.txt" --s $s --e $e
python ./project/ENDTE/graph_TE_addDEL.py --CHR "chr$i" --T2Tpos "./project/panVNTR/data/chr$i.posCHM.txt" --s $s --e $e

EOT
done
