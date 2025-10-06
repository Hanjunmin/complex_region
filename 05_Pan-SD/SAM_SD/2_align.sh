cd  .PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/
snakemake -s .PANSDEND/02_SAMSD/02_00-REFSD/run2new.smk --keep-going --keep-incomplete -j 40


mkdir other ###这里存放候选区域除了T2T区域的SD
while IFS= read -r line; do
echo $line
bedtools subtract -a <(cat $line/can/rm_in $line/can/groupin.ana) -b <(cat $line/t2tsdminimap/T2Thave.bed   $line/t2tsdminimap/T2Tpolymor1.bed $line/t2tsdminimap/filin.bedend $line/can/SDend |awk 'OFS="\t"{print $1,$2,$3"\n"$4,$5,$6}') |bedtools sort |bedtools merge >./other/$line".other" &
done < ".PANSDEND/DATA/allsam"

mkdir T2TPANSD ###这里存T2T区域的SD
while IFS= read -r line; do
echo $line
cat $line/t2tsdminimap/T2Thave.bed   $line/t2tsdminimap/T2Tpolymor1.bed $line/t2tsdminimap/filin.bedend $line/can/SDend >./T2TPANSD/$line".t2tsd" &
done < ".PANSDEND/DATA/allsam"
