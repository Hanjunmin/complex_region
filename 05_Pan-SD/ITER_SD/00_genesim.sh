## sim
# simdir="./sim/"
# cd ./pangenome/WGSfea/allW/APG3allWnew/scripts
# mkdir allsamcansim && cd  allsamcansim 


###这里为根据window的sim获取每个样本更有可能是SD区域的情况
simdir="./PANSDEND/000_SD_sim/2k"    #"./pangenome/APG_seg/allsim/"
cd ./PANSDEND/000_SD_sim/2k/split/
for i in {1..22} X; do
mkdir "chr$i"
cd "chr$i"
less -S ${simdir}/"chr$i"sim.out |cut -f 2-543 >"chr$i".sim
rm *.txt
IFS=$'\t' read -r -a col_names < <(head -n 1 "chr$i".sim)
cols=${#col_names[@]}
(for ((j=1; j<=cols; j++)); do
    cut -f$j "chr$i.sim" > "${col_names[j-1]}.txt"
done) &
cd ..
done


cd ./PANSDEND/000_SD_sim/2k/allsim #./pangenome/APG_seg/allsim
mkdir out
cd out
acro="./PANSDEND/DATA/acro.region"
while read line; do
  echo "$line"
  line=$(echo "$line" | tr -d '\r')
  find ./PANSDEND/000_SD_sim/2k/split/ -name "$line.txt" | sort > mat_files.txt  ##./pangenome/APG_seg/allsim/split/
  find ./PANSDEND/000_SD_sim/2k/split/ -name "win.txt" | sort > win_files.txt
  xargs -a mat_files.txt cat>x
  xargs -a win_files.txt cat >y
  paste y x  >$line.temp
  less -S $line.temp |tail -n +2 |tr ":" "\t" |tr "-" "\t" |less -S  |grep -v "e" |grep -v "win"  |awk '$4>=0.05' |less -S |bedtools sort |bedtools merge  >$line.sim.sd.bed
  bedtools subtract -a $line.sim.sd.bed -b ${acro} >$line.sim.bed.sd
  rm $line.temp
  rm $line.sim.sd.bed
done < <(cat ./DATA/APG/allsam |grep -v "CHM13v2")
#./pangenome/Pan-SD/allsam

cat *sd |less -S |bedtools sort |bedtools merge >cores.pos ##所有的对应的T2TCHM13的pos位置
bedtools intersect -a ./PANSDEND/000_SD_sim/2k/allsim/out/cores.pos -b ./sim/CHM13_APG_HPRC_HGSVC.pos.txt -wa -wb >xx
less -S cores.pos.T2T |cut -f 4-8 >CHM13_APG_HPRC_HGSVC.pos_iter.txt

./pangenome/APG_seg/allsim/out/传到了./sim/split这里了

cd ./PANSDEND/000_SD_sim/2k/split #./sim/split
#T2Tsegpos="./sim/CHM13_APG_HPRC_HGSVC.pos.txt"
mkdir simpath
cd simpath
cat ./pangenome/WGSfea/allW/APG3allW/scripts/allsam |grep -v "CHM13v2" |tail -n +250 | xargs -n 1 -P 30 -I {} bash -c '
  T2Tsegpos="./PANSDEND/000_SD_sim/2k/allsim/out/CHM13_APG_HPRC_HGSVC.pos_iter.txt"
  line=$(echo "{}" | tr -d "\r")
  echo "$line"
  bedtools intersect -a ./PANSDEND/000_SD_sim/2k/allsim/out/"$line".sim.bed.sd -b ${T2Tsegpos} -wa -wb | awk "OFS=\"\t\"{print \$1\":\"\$2\"-\"\$3,\$4,\$5,\$6,\$7,\$8}" > "$line"siminter.txt
  Rscript ./pangenome/WGSfea/allW/SDgenepath.r "$line"siminter.txt "$line"alldf.txt
  less -S "$line"alldf.txt | tail -n +2 | less -S > "$line"alldf.e.txt
  rm "$line"alldf.txt
  rm "$line"siminter.txt
'

cd ./PANSDEND/02_SAMSD/02_02-ITERSD/00_gene_sim/  ## tmux at -t 12
snakemake -s ./pangenome/WGSfea/allW/APG3allWnew/scripts/winsim1.smk --keep-going --keep-incomplete -j 30
snakemake -s ./pangenome/WGSfea/allW/APG3allWnew/scripts/winsim2.smk --keep-going --keep-incomplete -j 30
