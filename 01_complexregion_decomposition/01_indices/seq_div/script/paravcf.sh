filein=$1   #./allvcf/chr${m}.vcf
m=$2
cat $filein | awk -F'\t' '{
if($0 ~/^#/){
next
} else {
OFS="\t"
maxlen=0
split($5,spres,",")
same_vals = 1
for(i in spres){
len=length(spres[i])
if(len>maxlen){
maxlen=len
}
if(spres[i] != spres[1]){
same_vals = 0
}
}
if (same_vals == 1){  # query have same varients
if((maxlen-length($4))==0 && maxlen==1){
print $1,$2,$3,$4,$5,length($4),maxlen,"SNV"
}else if(length($4)!=1 && maxlen>1 && maxlen<50 && length($4)<50){
print $1,$2,$3,$4,$5,length($4),maxlen,"smallNM"
}else if(maxlen!=1 && length($4)!=1  && (length($4)>=50 || maxlen>=50)){
print $1,$2,$3,$4,$5,length($4),maxlen,"bigNM"
}else if(length($4)==1 && maxlen-length($4)>0 && maxlen-length($4)<50){
print $1,$2,$3,$4,$5,length($4),maxlen,"indel"
}else if(maxlen==1 && length($4)-maxlen>0 && length($4)-maxlen<50){
print $1,$2,$3,$4,$5,length($4),maxlen,"indel"
}else if(length($4)==1 && maxlen-length($4)>=50){
print $1,$2,$3,$4,$5,length($4),maxlen,"INS"
}else if(maxlen==1 && length($4)-maxlen>=50){
print $1,$2,$3,$4,$5,length($4),maxlen,"DEL"
}

}else{
if (maxlen < 50 && length($4) < 50) {
print $1,$2,$3,$4,$5,length($4),maxlen,"smallcomplex"
}else{
print $1,$2,$3,$4,$5,length($4),maxlen,"bigcomplex"
}

}
}
}' >spe.vcf


cat $filein |awk -F'\t' '{
	if($0 ~/^#/){
		next
	}else{
		OFS="\t"
		maxlen=0
		split($5,spres,",")
		for(i in spres){
			len=length(spres[i])
			if(len>maxlen){
				maxlen=len
			}
		}
		print $1,$2,$2+length($4)
	}
}' >vcf.bed 


awk '{print $8}' spe.vcf| sort | uniq -c
types=("SNV" "indel" "INS" "DEL" "smallNM" "bigNM" "smallcomplex" "bigcomplex")
for type in "${types[@]}"; do
    less -S spe.vcf | awk -v type="$type" '$8==type {OFS="\t"; print $1, $2, $2+length($4)}' > "${type}.bed"
done
types=("SNV" "indel" "INS" "DEL" "smallNM" "bigNM" "smallcomplex" "bigcomplex")
for type in "${types[@]}"; do
	bedtools coverage -a  chr${m}.size -b  "${type}.bed" > "${type}.window.bed"
done

cat spe.vcf |grep -E "indel|smallNM|smallcomplex" |awk 'OFS="\t"{print $1,$2,$2+length($4),$6,$7,$8}' |less -S >yy
bedtools intersect -a chr${m}.size -b yy -wa -wb >y
Rscript -e '
library(data.table)
data<-fread("y")
library(dplyr)
library(data.table)
setDT(data)
data$newval=pmax(data$V7,data$V8)
data$new=data$newval/50
result <- data[, Sum_group := sum(new), by = .(V1, V2, V3)]
result1 <- result[, maxvalue := mean(newval), by = .(V1, V2, V3)]
fwrite(distinct(result[,c("V1","V2","V3","Sum_group","maxvalue")]),"add.txt",sep="\t",col.names = FALSE)
'

## 3.Features:entropy
# cd ./pangenome/vcf/entropy/
cat $filein  |grep -v '^#' |cut -f 10-600 >py_entrin.txt
rm output.txt
python ./allvcf/allfeature/entro_cal.py
cat $filein   |grep -v '^#' >chm13v1.1_mc_allvcf.bed
less -S chm13v1.1_mc_allvcf.bed | awk 'OFS="\t" {
    split($5, arr, ","); 
    max_len = length($4); 
    for (i in arr) {
        if (length(arr[i]) > max_len) {
            max_len = length(arr[i]);
        }
    }
    print $1, $2, $2 + length($4), max_len
}' > lef.bed
#less -S chm13v1.1_mc_allvcf.bed | awk 'OFS="\t" {print $1, $2, $2+length($4), (length($4) > length($5) ? length($4) : length($5))}'>lef.bed
# less -S chm13v1.1_mc_allvcf.bed |awk 'OFS="\t"{print $1,$2,$2+length($4)}' >lef.bed
paste <(less -S lef.bed) <(less -S output.txt |awk '{print $NF}') >entropy_bed.txt
bedtools intersect -a chr${m}.size -b entropy_bed.txt -wa -wb >x
# bedtools map -a  chr${m}.size -b <(cat entropy_bed.txt|bedtools sort) -c 4 -o max >entropyend.txt
Rscript -e '
library(data.table)
data<-fread("x")
library(dplyr)
library(data.table)
setDT(data)
result <- data[, Sum_group := sum(V7), by = .(V1, V2, V3)]
result$weight=result$V7/result$Sum_group
result$weightentro=result$weight*result$V8
end <- result[, endentro := sum(weightentro), by = .(V1, V2, V3)]

fwrite(end,"testAL.txt",sep="\t")
fwrite(end[,c("V1","V2","V3","endentro")],"test.txt",sep="\t",col.names = FALSE)
'


cat test.txt |uniq >entropyendn.txt


# bedtools map -a  chr${m}.size -b <(cat entropy_bed.txt|bedtools sort) -c 4 -o max >entropyend.txt

## 4.allsv
bedtools coverage -a  chr${m}.size -b vcf.bed   >end.txt





less -S spe.vcf  |awk 'OFS="\t"{print $1,$2,$2+$6,$6,$7}' >1
bedtools intersect -a chr${m}.size  -b 1  -wa -wb >2

# for i in {1..22} X Y;do
# cd ./COMPLEXFEATURE/00_vcf/chr${i}
Rscript -e 'library(data.table)
data<-fread("2")
data$newwindos=data$V2
data$newwindoe=data$V3
data[data$V5<data$V2,]$newwindos=data[data$V5<data$V2,]$V5
data[data$V6>data$V3,]$newwindoe=data[data$V6>data$V3,]$V6
data$maxal=pmax(data$V7,data$V8)
setDT(data)
# result <- data[, .(alignlen = sum(maxal), 
# 				   diff = sum(V7), 
#                    ns = min(newwindos), 
#                    ne = max(newwindoe)), 
#                by = .(V1, V2, V3)]

result <- data[, .(
  alignlen   = sum(maxal),
  diff       = sum(V7),
  ns         = min(newwindos),
  ne         = max(newwindoe),
  range_flat = list(unlist(mapply(`:`, V5, V6))),
  range_len  = length(unique(unlist(mapply(`:`, V5, V6))))  # 正确使用 length
), by = .(V1, V2, V3)]


result[, nlen := length(unique(unlist(range_flat))), by = .(V1, V2, V3)]

result$match=(result$ne-result$ns+1)-result$nlen
result$sim=result$match/(result$match+result$alignlen)
fwrite(result[,c("V1","V2","V3","sim")],"simal",sep="\t")

# plot(density(result$sim))'
# done
