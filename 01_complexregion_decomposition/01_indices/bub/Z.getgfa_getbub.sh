
#####STEP2.按照2M进行划分 加快之后的1k窗口的计算速度

####获取2M区域的gfa和对应的region左右更改
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/
bedtools makewindows -g <(less -S ./breakpoints/chm13v2.0.fa.fai  |awk 'OFS="\t"{print $1,$2}') -w 2000000 -s 2000000  >WGS.window2M.bed
sed 's/\b0\b/1/g'  WGS.window2M.bed >WGS.window2M1.bed
for i in {1..22} X Y; do
(
mkdir -p ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"
less -S ./COMPLEXFEATURE/09_bubble_only/01_GFA/WGS.window2M1.bed  |awk -v chr="chr$i" '$1==chr{print $0}'  >"chr$i".window
[ -f Sregional.node ] && rm Sregional.node
[ -f region.bed ] && rm region.bed 
befst1=0
befen1=0
while read -r chr start end; do
  reg=$chr":"$start"-"$end
  if [ "$start" -gt "$befst1" ] && [ "$end" -lt "$befen1" ]; then
    echo -e "$chr\t$start\t$end\t$befst1\t$befen1\t"extend"" >> region.bed
    continue
  fi
  ./pangenome/software/gfabase sub ./pangenome/CHM13-APGp1-HPRCp1-HGSVCp3_MC.gfab CHM13v2#0#$chr":"$start"-"$end --range --view --cutpoints 1 --guess-ranges -o $chr":"$start"-"$end.gfa
  st1=$(less -S $chr":"$start"-"$end.gfa  |grep "CHM13v2#" |less -S  |cut -f8 |sed 's/,//g' |tr ":" "\t" |tr "-" "\t" |cut -f4 |sort -k1,1n |head -n 1)
  en1=$(less -S $chr":"$start"-"$end.gfa  |grep "CHM13v2#" |less -S  |cut -f8 |sed 's/,//g' |tr ":" "\t" |tr "-" "\t" |cut -f5 |sort -k1,1n |tail -n 1)
  cat  $chr":"$start"-"$end.gfa  |awk '$1=="S"' |less -S |cut -f 1-2 |awk  -v  reg=$reg 'OFS="\t"{print $1,$2,reg}' |sort -k2,2n |uniq >>Sregional.node
  echo -e "$chr\t$start\t$end\t$st1\t$en1" >> region.bed
  befst1=$st1
  befen1=$en1
done < "chr$i".window
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/
)&
done


####根据2M区域的gfa和对应的region获取相应区域的所有样本的Path
for i in {1..22} X Y;do
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"
cat region.bed >region.bedbef
cat region.bedbef |grep -v "extend" >region.bed
cat ./PANSDEND/02_SAMSD/02_03-VALI/00_PAIR/*cor.pos |awk -v chr="chr$i" '$2==chr' >pair.bed
bash ./COMPLEXFEATURE/09_bubble_only/01_GFA/2M_para.sh #1>log.txt 2>err.txt

done

####从0到29的Path.txt中按照区域进行汇总
for i in {16..22} X Y;do
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"
rm *.path
cat region.bedbef |grep -v "extend" >region.bed
while read -r chr start end x y z; do
  reg=$chr":"$start"-"$end
  (
  for i in {0..29};do
      echo $i
      line_numbers=$(grep -n $reg ${i}path.index | cut -d: -f1)
      if [ -z "$line_numbers" ]; then
      continue
    else
   sed -n "$(echo $line_numbers | sed 's/ /p;/g')p" ${i}path.txt >>${reg}.path
    fi
      
  done)&
done < region.bed
wait
done


####获取原始的>2k的region 然后根据2M的path进行分割计算每个小region(真正region的path)
for i in {1..22} X Y;do
cat ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/region.bed | xargs -I {} -P 10 bash ./COMPLEXFEATURE/09_bubble_only/01_GFA/subpara.sh "{}" "$i" 
done

for i in {1..22} X Y;do
cat ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/region.bed | xargs -I {} -P 300 bash ./COMPLEXFEATURE/09_bubble_only/01_GFA/subpara.inter.sh "{}" "$i" 
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"
find . -name "inter.node" | xargs -I {} awk -v folder="{}" '{OFS="\t"}BEGIN { split(folder, parts, "/"); folder_name=parts[length(parts)-1] } { $4 = folder_name; print }' {} |bedtools sort >inter.nodeal.bef
find . -name "inter.node" | xargs -I {} awk -v folder="{}" '{OFS="\t"}BEGIN { split(folder, parts, "/"); folder_name=parts[length(parts)-1] } { $4 = folder_name; print }' {} |bedtools sort|awk '$3-$2<2000000' >inter.nodeal
less -S inter.nodeal |cut -f 1-3 |sort |uniq -c |awk '$1!=1' |awk 'OFS="\t"{print $2,$3,$4}' >gfa.in
len=$(cat gfa.in |wc -l)
if [ ${len} -eq 0 ]; then
    continue
fi
mkdir -p ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/inter
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/inter
[ -f Sregional.node ] && rm Sregional.node
while read -r chr start end; do
reg=$chr":"$start"-"$end
./pangenome/software/gfabase sub ./pangenome/CHM13-APGp1-HPRCp1-HGSVCp3_MC.gfab CHM13v2#0#$chr":"$start"-"$end --range --view --cutpoints 1 --guess-ranges -o $chr":"$start"-"$end.gfa
cat  $chr":"$start"-"$end.gfa  |awk '$1=="S"' |less -S |cut -f 1-2 |awk  -v  reg=$reg 'OFS="\t"{print $1,$2,reg}' |sort -k2,2n |uniq >>Sregional.node
 done < ../gfa.in
 cp ../inter.nodeal ./
 cp ../region.bed ./
python ./COMPLEXFEATURE/09_bubble_only/01_GFA/subregion_inter.py
done



####!!!!!!!!!!!!!!!!!!!生成一个文件 这个文件为2M区域和subregion的整合文件

in="./COMPLEXFEATURE/09_bubble_only/01_GFA/al.subregional"
[ -f $in ] && rm $in
for i in {1..22} X Y;do
  while read -r chr start end a b; do
  reg=$chr":"$start"-"$end
  echo $reg
  cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/$chr/$reg
  cat $reg".subregion" |awk -v reg=$reg 'OFS="\t"{print reg,$0}' >> $in
done <./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/region.bed
done

[ -f $in2 ] && rm $in2
in2="./COMPLEXFEATURE/09_bubble_only/01_GFA/al.subregional2"

cat  ./COMPLEXFEATURE/09_bubble_only/01_GFA/*/gfa.in >x
less -S x |awk 'OFS="\t"{print "inter",$0}' >$in2

cat $in $in2  |awk '$4-$3<1000000' >./COMPLEXFEATURE/09_bubble_only/01_GFA/al.subregional.1M


#####加入T2TCHM13path
for i in {1..22} X Y;do
  while read -r chr start end a b c; do
    reg=$chr":"$start"-"$end
  echo $reg
  cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/$chr/$reg
  {
    [ -f CHM13.path ] && rm CHM13.path
    while read -r subchr substart subend; do
      x=$subchr":"$substart"-"$subend 
      less -S ${x}.gfa |grep "CHM" |grep "^S" |cut -f2,8 |sed "s/,//g" |sed "s/#/\t/g" |cut -f 1,4 |tr ":" "\t" |tr "-" "\t" |sort -k3,3n |cut -f1|tr "\n" "," |sed "s/,/+,/g" |awk -v reg=$x 'OFS="\t"{print "P",reg"@CHM13",$0,"*"}' >>CHM13.path
    done<$reg".subregion" 
  }&
  done< ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/region.bed
done


for i in {1..22} X Y;do
  cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/chr${i}/inter
  {
    [ -f CHM13.path ] && rm CHM13.path
    while read -r subchr substart subend; do
      x=$subchr":"$substart"-"$subend 
      less -S ${x}.gfa |grep "CHM" |grep "^S" |cut -f2,8 |sed "s/,//g" |sed "s/#/\t/g" |cut -f 1,4 |tr ":" "\t" |tr "-" "\t" |sort -k3,3n |cut -f1|tr "\n" "," |sed "s/,/+,/g" |awk -v reg=$x 'OFS="\t"{print "P",reg"@CHM13",$0,"*"}' >>CHM13.path
    done<../gfa.in
  }&
done




###计算每个小region对应的路径和过滤之后的路径
for i in {1..22} X Y;do
cat ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr$i"/region.bed| xargs -I {} -P 150 bash -c '
  chr=$(echo -e "{}" | cut -f1)
  start=$(echo -e "{}" | cut -f2)
  end=$(echo -e "{}" | cut -f3)
  reg=$chr":"$start"-"$end
  echo $reg
  cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/$chr/$reg
  python ./COMPLEXFEATURE/09_bubble_only/01_GFA/high_freseg.py 2>log${reg}_err.txt
'
done

for i in {1..22} X Y;do
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/chr${i}/inter
python ./COMPLEXFEATURE/09_bubble_only/01_GFA/high_freseg.py
done



####计算具体的segment
##替代下面这俩 inter（覆盖多个2M）和具体的
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA
snakemake -s ./COMPLEXFEATURE/09_bubble_only/01_GFA/para_path.smk --keep-going --keep-incomplete -j 100



####判断结果
file="./COMPLEXFEATURE/09_bubble_only/01_GFA/pathin"
ls */polycomplex* >y
less -S y |sed "s/polycomplex/\t/g" |sed "s/.txt//g" |cut -f2 >have
less -S $file |awk '$1=="chr3"' |awk '{print $1":"$2"-"$3}' >al
grep -Fxv -f have al  |tr ":" "\t" |tr "-" "\t">weno
bedtools coverage -a weno -b ~/allvcf/chr1.bed >xx
less -S xx |awk '$7!=0' |wc -l

####总结结果
for i in {1..22} X Y;do
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/chr$i
ls */polycomplex* >chr${i}.pm.bed &
done


for i in {1..22} X Y;do
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/chr$i
[ -f chr${i}.end.sta ] && rm chr${i}.end.sta
cat chr${i}.pm.bed  |sed "s/polycomplex/\t/g" |sed 's/.txt//g' >chr${i}.pathregion.bed
less -S chr${i}.pathregion.bed |cut -f2 |tr ":" "\t" |tr "-" "\t" >chr${i}.pathregion.se.bed
paste chr${i}.pm.bed chr${i}.pathregion.bed chr${i}.pathregion.se.bed>chr${i}.pm_region.bed
{
  while read -r subfold fold subregion chr start end; do
echo $subregion
echo $fold
echo $subfold
echo $chr
echo $end
poly=$(cat  $subfold)
nodenum=$(cat  ${fold}/${subregion}.sim.gfa|awk '$1=="S"'  |wc -l)
intenode=$(cat ${fold}/simnode${subregion}.txt|cut -f1)
polyratio=$(cat ${fold}/polyratio${subregion}.txt|cut -f1)
polyratio=$(printf "%.20f" "$polyratio")
subpath=$(cat ${fold}/${subregion}.path.dou.filt|wc -l)
allpath=$(cat ${fold}/${subregion}.path.dou.bef|wc -l)
pathdif=$(echo "scale=2; $subpath / $allpath" | bc)
sumnode=$(cat ${fold}/${subregion}.gfa  |grep '^S' |cut -f 2-3|  awk 'length($2) != 1'  |awk '{sum += length($2)} END {print sum}')
truepos=$(echo "$end - $start" | bc)
min_val=$(echo "if ($sumnode < $truepos) $sumnode else $truepos" | bc)
max_val=$(echo "if ($sumnode >= $truepos) $sumnode else $truepos" | bc)
refinelen=$(echo "scale=2; $min_val / $max_val" | bc)
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$chr" "$start" "$end" "$poly" "$polyratio" "$intenode" "$pathdif" "$refinelen" ${fold}${subregion}>> chr${i}.end.sta
done<chr${i}.pm_region.bed
} &
done


cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/
cat  */*.end.sta >al.sta




#####STEP1.生成每个1k窗口对应的gfa序列以及最后的end文件
for i in {1..22} X Y; do
        mkdir "chr$i"
        cd "chr$i"
        rm *
        less -S ./pangenome/WGSfea/WGS.window  |awk -v chr="chr$i" '$1==chr{print $0}'  >"chr$i".window
        split -l 50000 chr${i}.window
        for file in x*; do
            echo "$file"
            nohup bash ./COMPLEXFEATURE/09_bubble_only/01_GFA/para.sh  "$file" &
        done
        cd ..
done

cat */*WGS_bubble_APG.txt |cut -f 1-6>end
#####
library(data.table)
data<-fread("end")
library(dplyr)
library(stringr)

result <- data %>%
  group_by(V4, V5, V1) %>%
  summarise(
    min_col1 = min(V2),
    max_col3 = max(V3),
    combined_V6 = paste(V6, collapse = ", "),  # 将 V6 列的值连接成一个字符串
    .groups = 'drop'
  )

res <- result %>%
  filter(str_detect(combined_V6, "extend"))
res=res[,c(3,4,5)]
fwrite(res,"test",sep="\t")


less -S test |tail -n +2 |bedtools sort |bedtools merge >pathin
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/para_gene_path/
for chr in {1..22} X Y;do
    mkdir -p ./COMPLEXFEATURE/09_bubble_only/01_GFA/para_gene_path/"chr"$chr
    cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr"$chr
    cat ./COMPLEXFEATURE/09_bubble_only/01_GFA/pathin|awk -v chrn="chr"$chr '{OFS="\t"}$1==chrn' >region.bed
done


for chr in {1..22} X Y;do
cd ./COMPLEXFEATURE/09_bubble_only/01_GFA/"chr"$chr
bash ./COMPLEXFEATURE/01_bubble/sub_paraa.sh #1>log.txt 2>err.txt

done



####STEP3.INTEGRATE
cd ./COMPLEXFEATURE/09_bubble_only/02_INTE/
bedtools intersect -b ./COMPLEXFEATURE/09_bubble_only/01_GFA/al.sta -a ./pangenome/WGSfea/WGS.window  -f 1 -wa  -wb >x
bedtools subtract -a ./COMPLEXFEATURE/09_bubble_only/01_GFA/pathin -b ./COMPLEXFEATURE/09_bubble_only/01_GFA/al.sta >weno.regio

library(data.table) ###原先的加权重
data <- read.table("x", header = FALSE, sep = "\t")
sigmoid <- function(x, k, c) {
  1 / (1 + exp(-k * (x - c)))
}

# data$lenrefine=1-sigmoid(data$V6,300,0.85)
data$polyrefine=sigmoid(data$V8,100,0.1) ###polyratio

data$H=0
data[data$V10<0.5,]$H=1 ###path相似性
data[data$V10>0.5,]$H=0



data$V4new=2*sigmoid(data$polyrefine*data$V7,3,3) ###V7 poly
data$V5new=sigmoid(data$V9,3,3) ####simplenode
data$new=data$V4new+data$V5new
data$new=data$new+1*data$H


data <- na.omit(data)
data$newnorm=(data$new-min(data$new))/(max(data$new)-min(data$new))
fwrite(data,"end.txt",sep="\t",col.names=FALSE)
# na_rows <- data[is.na(data$new), ]
# print(na_rows)
# fwrite(nan_rows,"nan.bed",sep="\t")

# end<-rbind(x,y)
# # fwrite(end[,c(6,7,8,4,5,9,10,11,17,18,19,1,2,3)],"end.txt",sep="\t",col.names=FALSE)
# fwrite(end[,c(9,10,11,12,13,14,4,6,22,23,24,25,26,27,1,2,3)],"end.txt",sep="\t",col.names=FALSE)
# fwrite(end[,c(1,2,3,4,6,20)],"test.txt",sep="\t",col.names=FALSE)

# less -S test.txt  |awk '$6>0.2' |cut -f 1-3 |bedtools sort |bedtools merge>bubbletrue.region



# bedtools intersect -a ./COMPLEXFEATURE/06_benchmark/add_buble_coplex/00_HPRC/hprc.end.38.convert.t2t.bed -b end.txt -wa -wb|less -S



# ##########03 model 
# python ./COMPLEXFEATURE/08_model/add_bubble.model.py
# cd ./COMPLEXFEATURE/09_bubble_only/03_model
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py --m KM  --c 5 --b 1 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py --m KM  --c 6 --b 1 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py --m KM  --c 7 --b 1 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py --m KM  --c 5 --b 2 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py --m KM  --c 6 --b 2 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m KM  --c 7 --b 2 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m FCC  --c 5 --b 1 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m FCC  --c 6 --b 1 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m FCC  --c 7 --b 1 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m FCC  --c 5 --b 2 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m FCC  --c 6 --b 2 &
# python ./COMPLEXFEATURE/09_bubble_only/03_model/add_bubble.model.py  --m FCC  --c 7 --b 2 &


