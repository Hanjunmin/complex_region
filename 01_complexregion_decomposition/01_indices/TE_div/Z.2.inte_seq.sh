for i in {1..22} X Y; do
    cd ./project/ENDTE/00_INS/chr${i}
    cat *.csv |grep "@" |grep -v "DEL" >seq.in &
    cat *.csv |grep "@" |grep "DEL" >DEL.in &
done

for i in {1..22} X Y; do
    cd ./project/ENDTE/00_INS/chr${i}
    python ./project/ENDTE/gene_seq.py &
    # samtools faidx output.fasta
done

for i in {1..22} X Y; do
    cd ./project/ENDTE/00_INS/chr${i}
    rm split*
    total_lines=$(wc -l < output.fasta)
    lines_per_file=$(( (total_lines + 9) / 10 ))
    split -l $lines_per_file -d -a 1 output.fasta split_
done


bash ./project/ENDTE/2_run.rm.sh

