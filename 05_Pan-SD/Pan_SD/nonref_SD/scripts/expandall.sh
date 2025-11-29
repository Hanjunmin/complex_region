### expand the theresehod
## input (allclus)
T2Trmtrfcut=$1
assembyfai=$2 ##/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/t2tsdminimap/HG00438_Pat.fa.fai
processscript=$3 ##/home/jmhan/pangenome/WGSfea/sdseq/expand_pansd.r

#T2Trmtrfcut="/home/jmhan/pangenome/WGSfea/allW/APG3allWnew/scripts/t2trmtrf.cut1M.bed"
less -S allclus |awk '$3-$2>=500 && $6-$5>=500' >allclusnew
cat allclusnew | cut -f 1-3 | parallel --pipe -j 10 --block 50M 'sort -u' > lef
less -S lef |sort -u >lefn
cat allclusnew | cut -f 4-6 | parallel --pipe -j 10 --block 50M 'sort -u' > rig
less -S rig |sort -u >rign
sed -i 's/_complex//g' winsim_t2t.bed

###下面这个的列分别为：想扩展的区域  原始候选区域 候选区域对应的T2T坐标（大概）
bedtools intersect -b <(cat  winsim_t2t.bed  |awk 'OFS="\t"{print $4,$5,$6,$1,$2,$3}') -a lefn -wa -wb >lef1
bedtools intersect -b <(cat  winsim_t2t.bed  |awk 'OFS="\t"{print $4,$5,$6,$1,$2,$3}') -a rign -wa -wb >rig1
#R1
Rscript -e '
library(data.table)
library(dplyr)
data<-fread("lef1",header=FALSE,sep="\t")
result <- data %>%
group_by(V1,V2,V3) %>%  # 将 V4 直接用作列名
summarize(
  can_s =min(V5),
  can_e = max(V6),
  t2t_chr = first(V7),
  t2t_s = min(V8),
  t2t_e = max(V9),
)
fwrite(result,"lef1.out", sep = "\t",col.names=TRUE) '

Rscript -e '
library(data.table)
library(dplyr)
data<-fread("rig1",header=FALSE,sep="\t")
result <- data %>%
group_by(V1,V2,V3) %>%  # 将 V4 直接用作列名
summarize(
  can_s =min(V5),
  can_e = max(V6),
  t2t_chr = first(V7),
  t2t_s = min(V8),
  t2t_e = max(V9),
)
fwrite(result,"rig1.out", sep = "\t",col.names=TRUE) '


bedtools subtract -a <(less -S lef1.out  |awk 'OFS="\t"{print $6,$7,$8,$1,$2,$3,$4,$5}' |tail -n +2 ) -b $T2Trmtrfcut  -A |cut -f 4-6 >lefnointer.out
bedtools intersect -a <(less -S lef1.out  |awk 'OFS="\t"{print $6,$7,$8,$1,$2,$3,$4,$5}' |tail -n +2 )  -b $T2Trmtrfcut  -wa -wb|less -S >lefin

bedtools subtract -a <(less -S rig1.out  |awk 'OFS="\t"{print $6,$7,$8,$1,$2,$3,$4,$5}' |tail -n +2 ) -b $T2Trmtrfcut  -A |cut -f 4-6 >rignointer.out
bedtools intersect -a <(less -S rig1.out  |awk 'OFS="\t"{print $6,$7,$8,$1,$2,$3,$4,$5}' |tail -n +2 )  -b $T2Trmtrfcut  -wa -wb|less -S >rigin

### R2
Rscript $processscript lefin $assembyfai lefin.expand
Rscript $processscript rigin $assembyfai rigin.expand

cat <(cat lefnointer.out |awk 'OFS="\t"{print $1,$2,$3,$2,$3}') lefin.expand > lefin.expand.all
cat <(cat rignointer.out |awk 'OFS="\t"{print $1,$2,$3,$2,$3}') rigin.expand > rigin.expand.all


###R
Rscript -e '
library(data.table)
library(dplyr)
data<-fread("allclusnew",header=FALSE,sep="\t")
colnames(data)<-c("lefchr","lefsinit","lefeinit","rigchr","rigsinit","rigeinit")
expandlef<-fread("lefin.expand.all",header=FALSE,sep="\t")
colnames(expandlef)<-c("lefchr","lefsinit","lefeinit","lefsn","lefen")
expandrigh<-fread("rigin.expand.all",header=FALSE,sep="\t")
colnames(expandrigh)<-c("rigchr","rigsinit","rigeinit","rigsn","rigen")
end1 <- merge(data,expandlef, by = c("lefchr","lefsinit","lefeinit"), all.x = TRUE)
end1 <- merge(end1, expandrigh, by = c("rigchr","rigsinit","rigeinit"), all.x = TRUE)
fwrite(end1[,c("lefchr","lefsn","lefen","rigchr","rigsn","rigen")],"expand.out", sep = "\t",col.names=FALSE) 
'
rm lefn
rm rign 