



# while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12; do
# minimap2 -c --eqx  <(samtools faidx ~/breakpoints/chm13v2.0.fa ${col1}:${col2}-${col3}) <(samtools faidx ~/breakpoints/chm13v2.0.fa  ${col4}:${col5}-${col6}) -A 5 -B 4 -O 40 -E 1 -s 3000 >>x.paf
# done <  REF.end

awk 'OFS="\t"{diff1 = $3 - $2; diff2 = $6 - $5; min_diff = (diff1 < diff2) ? diff1 : diff2; max_diff = (diff1 > diff2) ? diff1 : diff2; if (max_diff != 0) print $0,min_diff / max_diff}' .PANSDEND/000_SD_sim/2k/T2TSD/filt/split/END/REF.end |awk '$13<0.85' |cut -f 1-12 >p85s.sd
awk 'OFS="\t"{diff1 = $3 - $2; diff2 = $6 - $5; min_diff = (diff1 < diff2) ? diff1 : diff2; max_diff = (diff1 > diff2) ? diff1 : diff2; if (max_diff != 0) print $0,min_diff / max_diff}' .PANSDEND/000_SD_sim/2k/T2TSD/filt/split/END/REF.end |awk '$13>=0.85'  |cut -f 1-12>p85l.sd



mkdir para
cd .PANSDEND/01_REFSD/para
export REF_GENOME=~/breakpoints/chm13v2.0.fa
split -l $(( $(wc -l < ../p85s.sd) / 100 )) ../p85s.sd split_part_
for file in ./split_part_*; do
    filename=$(basename "$file")  # 获取文件名
    suffix="${filename#split_part_}"  # 去掉前缀"split_part_"
    while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12; do
    minimap2 -c --eqx  <(samtools faidx ~/breakpoints/chm13v2.0.fa ${col1}:${col2}-${col3}) <(samtools faidx ~/breakpoints/chm13v2.0.fa  ${col4}:${col5}-${col6}) -A 5 -B 4 -O 40 -E 1 -s 3000 >>$suffix.paf
    done <  $file &
done


cat *.paf  |cut -f 1-9>pafall 


library(data.table)
data<-fread("pafall",header=FALSE)
library(dplyr)
grouped_data <- data %>%
  group_by(V1, V6) %>%
  group_split()

refine=data.frame()
for (group in grouped_data) {
  if(nrow(group)!=1){
    group$lef=group$V4-group$V3
    group$rig=group$V9-group$V8
    leflen=group$V2[1]
    riglen=group$V2[1]
    x1=max(group$lef)/leflen
    x2=max(group$rig)/riglen
    if(x1>0.9 & x2>0.9){
        print("TES")
    }else{
        refine=rbind(refine,group)
    }
  }
}
refdata=fread("../REF.end")
refine=distinct(refine[,c(1,6)])
colnames(refine)=c("REFrig","REFlef")
refine$have=1
refdata$REFlef=paste(refdata$V1,":",refdata$V2,"-",refdata$V3,sep="")
refdata$REFrig=paste(refdata$V4,":",refdata$V5,"-",refdata$V6,sep="")
end<-merge(refdata,refine,by=c("REFlef","REFrig"),all.x=TRUE)
save=end[is.na(end$have),c("V1","V2","V3","V4","V5","V6" , "V7"  ,"V8",  "V9",  "V10" ,"V11" ,"V12")]
fwrite(save,"FIX.sd",sep="\t",col.names=FALSE)
fwrite(refine,"REFINE.sd",sep="\t",col.names=FALSE)

rm refine.paf
while IFS=$'\t' read -r col1 col2 col3; do
minimap2 -c --eqx  <(samtools faidx ~/breakpoints/chm13v2.0.fa ${col1}) <(samtools faidx ~/breakpoints/chm13v2.0.fa  ${col2}) -A 5 -B 4 -O 40 -E 1 -s 3000 >>refine.paf
done <  REFINE.sd

rm=".HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/asm_repeatmasker.out.bed"
rmtrf=".HPRC/integrate/sedef/allsedef/T2Tchm13/Masked/trf_rm.bed"
python .pangenome/WGSfea/allW/paf2bed.py  --paf  refine.paf  --o  all.end.statistics
cat all.end.statistics |cut -f 1-12 >all.end.statistics.12
bash .pangenome/WGSfea/allW/APG3allW/scripts/endfilt.sh $rm $rmtrf all.end.statistics.12 all.end.statistics.end
###小范围过滤-过滤有交集的
Rscript .PANSDEND/script/SR_filt.r 
cat FIX.sd <(cat sd.filt |cut -f 1-12) >SD.refine


.PANSDEND/01_REFSD/para/
less -S SD.refine  |awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6}' >edge.bed
