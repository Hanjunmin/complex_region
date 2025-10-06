### statistic the orient and sample assmbly information
cd ./project/APG/NPHP1/01_minimap/
echo -e 'chr2\t110000000\t111200000' >our.bed
samtools faidx  ~/breakpoints/chm13v2.0.fa chr2:110000000-111200000 >T2Tchmchr2.fa
file_path="./DATA/APG/allsam"
inte_corfold="./PANSDEND/02_SAMSD/02_06-VALI/00_PAIR"
[ -f cov.bed ] && rm cov.bed
while IFS= read -r sam; do
echo $sam
samnew=$(cat ./allfa/alignments/saffire/allsam.link |grep $sam |cut -f1)
saf="./allfa/alignments/saffire/${samnew}.map2chm13.saffire"
COV=$(bedtools coverage -a our.bed -b $saf   |cut -f7)
querynameA=$( bedtools intersect -a our.bed -b $saf -wa -wb  |less -S |sort -k12,12n |less -S   |tail -n 1 |cut -f9)
if [[ ! -n "$querynameA" ]] ;then
    echo -e "$sam\t$COV\tnointer" >>cov.bed
    continue
fi
ques=$(cat $saf |grep $querynameA |awk '$1=="chr2"' |less -S  |sort -k2,2n |head -n 1 |cut -f2)
quee=$(cat $saf |grep $querynameA |awk '$1=="chr2"' |less -S  |sort -k2,2n |tail -n 1 |cut -f3)
if [[ $ques -le 110000000 && $quee -ge 111200000 ]]; then
    add=$( cat $saf |grep $querynameA |awk '$1=="chr2"'  |awk '$5=="+"' |cal)
    minus=$( cat $saf |grep $querynameA |awk '$1=="chr2"'  |awk '$5=="-"' |cal)
    [ -z "$minus" ] && minus=0
    [ -z "$add" ] && add=0
    if [[ $add -ge $minus ]];then
        echo -e "$sam\t$COV\t$querynameA\t$add\t$minus\tadd" >>cov.bed
    else
        echo -e "$sam\t$COV\t$querynameA\t$add\t$minus\tminus" >>cov.bed
    fi
else
    echo -e "$sam\t$COV\tnocontinue..." >>cov.bed
fi

done < "$file_path"




cat cov.bed |grep -v "no" >now.sam
while IFS=$'\t' read -r sam cov name ad mi type; do
{
    echo $sam
    namn=$(echo $name |tr "#" "_")
    if [[ "$type" == "minus" ]]; then
    seqkit seq -r -p <(samtools faidx ./allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/${sam}.fa $namn) > ${sam}.temp.fa
    minimap2 -cx asm5  ${sam}.temp.fa T2Tchmchr2.fa -o "${sam}chr2_T2T.paf"
    else
    minimap2 -cx asm5  <(samtools faidx ./allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/${sam}.fa $namn) T2Tchmchr2.fa -o "${sam}chr2_T2T.paf"
    fi
} &
done < now.sam



while IFS=$'\t' read -r sam cov name ad mi type; do
    echo $sam
    rustybam orient "${sam}chr2_T2T.paf" | rustybam stats --paf > "${sam}chr2_T2T.saffire"
    less -S "${sam}chr2_T2T.saffire" | grep -v "name" | sort -k1 -k2,2n > "${sam}chr2_T2T.txt"
done < now.sam

# cat "$file_path" | xargs -I {} -P 50 bash -c '
# line="$1"
# minimap2 -cx asm5 "./allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/${line}.fa" T2Tchmchr2.fa -o "${line}chr2_T2T.paf"
# rustybam orient "${line}chr2_T2T.paf" | rustybam stats --paf > "${line}chr2_T2T.saffire"
# linum=$(less -S "${line}chr2_T2T.saffire" | cut -f1 | sort | uniq | wc -l)
# echo "####################"
# echo "$line"
# echo "$linum"
# less -S "${line}chr2_T2T.saffire" | grep -v "name" | sort -k1 -k2,2n > "${line}chr2_T2T.txt"
# ' _ {}


library(data.table)
library(GenomicRanges)
alldata<-data.frame()
allrow<-data.frame()
otherdata<-c()
file=fread("now.sam",header=FALSE)
for(sample in file$V1){
    print(sample)
    data<-fread(paste(sample,"chr2_T2T.txt",sep=""))
    filtdata<-data.frame()
    for(id in unique(data$V1)){
        x=data[data$V1==id,]
        gr <- GRanges(seqnames = rep("chrA", nrow(x)), ranges = IRanges(start = x$V7,end = x$V8))
        print(sum(width(reduce(gr))))
        if(sum(width(reduce(gr)))>500000){
            filtdata<-rbind(filtdata,x)
        }
    }
    if(nrow(filtdata)==0){
        print("Oh no! this sample don't have the match contig")
        otherdata<-c(otherdata,sample)
    }else{
        gr <- GRanges(seqnames = rep("chrB", nrow(filtdata)), ranges = IRanges(start = filtdata$V2,end = filtdata$V3))
        if(length(unique(filtdata$V1)) !=1 |sum(width(reduce(gr)))<=500000 | sum(width(reduce(gr)))>1500000 ){
        print("Oh no! this sample may error!!!!!!!!!!!!!")
        
        print(sum(width(reduce(gr))))
        
    }else{
        print("OK!!!")
        value=unique(filtdata$V6)
        ori=substr(value, nchar(value), nchar(value))
        row=c(sample,unique(filtdata$V1),min(filtdata$V2),max(filtdata$V3),min(filtdata$V7),max(filtdata$V8),ori)
        allrow<-rbind(allrow,row)
        alldata<-rbind(alldata,filtdata)
    }
}
}
 colnames(allrow)[1]="V1"
 en=merge(file[,c(1,6)],allrow,by="V1")

fwrite(alldata,"allsample_delo.txt",col.names=FALSE,sep="\t")
fwrite(en,"allprocess_delo.txt",col.names=FALSE,sep="\t")


mkdir fa && cd fa
while IFS=$'\t' read -r col1 seq col2 col3 col4 col5 col6 col7
do
    if [[ "$seq" == "minus" ]]; then
    samtools faidx  ./project/APG/NPHP1/01_minimap/${col1}.temp.fa  "$col2":"$col3"-"$col4"  >${col1}.fa
    else
     samtools faidx  ./allfa/CHM13-APGp1-HPRCp1-HGSVCp3fa/${col1}.fa  "$col2":"$col3"-"$col4" >${col1}.fa
    fi
    minimap2 -c --eqx ${col1}.fa ./pangenome/decompose/NPHP1/sd.fa -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -o ${col1}.paf
done < ../allprocess_delo.txt

samtools faidx  ~/breakpoints/chm13v2.0.fa chr2:110000000-111200000 >CHM13v2.fa
cd ./project/APG/NPHP1/02_PGGB/
file_path="./DATA/APG/allsam"
rm *
while IFS= read -r sam; do
if [ -f "./project/APG/NPHP1/01_minimap/fa/${sam}.fa" ]; then
   # samn=$(echo $sam|tr "-" "_")
    sed "s/^>.*$/>${sam}#1#x/" "./project/APG/NPHP1/01_minimap/fa/${sam}.fa" > "${sam}.fa"
fi
done < "$file_path"
# rm CN1v1.fa
# rm GRCh38.fa
cat *fa >../end.fa
samtools faidx ../end.fa

singularity run ./software/pggb_latest.sif
pggb -i ../end.fa  -n  519  -t 50   -V 'CHM13v2#1#x:1000'     -o output 
##unset -f which

less -S ./pangenome/SD/T2Tchm13v1.1SD.txt |awk '$1=="chr2" && $10=="chr2"' |awk '$2>=110000000 && $3<=112000000' |awk '$11>=110000000 && $12<=112000000' |awk '$3-$2>5000' |less -S
samtools faidx ./pangenome/decompose/NPHP1/sd.fa chr2:110197933-110391878 chr2:110397814-110517542 chr2:110520233-110565558 > now.sd.fa

cd ./project/APG/NPHP1/02_PGGB/output/
gfafile="end.fa.bf3285f.11fba48.85c7139.smooth.final.gfa"
GraphAligner -g ./project/APG/NPHP1/02_PGGB/output/${gfafile}  -f   ./project/APG/NPHP1/02_PGGB/now.sd.fa   -a NPHP1SD.gaf -x vg
samtools faidx  ~/breakpoints/chm13v2.0.fa  chr2:110588729-110670393 >NPHP1.fa
GraphAligner -g ./project/APG/NPHP1/02_PGGB/output/${gfafile}  -f   NPHP1.fa   -a NPHP1gene.gaf -x vg

gfatools  gfa2fa   ${gfafile} >segment.fa
python   ./software/pggb/NPHP1/NPHP1/FINAL_GFA/gen_segemntdic.py   --fa  segment.fa  --o  segemntdic.json.gz
cat NPHP1SD.gaf >NPHP1SD.new.gaf
cat NPHP1SD.gaf |sed 's/</\t</g' |sed 's/>/\t>/g' >NPHP1SDt.new.paf
# Rscript   ./software/pggb/NPHP1/NPHP1/FINAL_GFA/filtSD.r
less -S ${gfafile}|awk '$1=="P"' |less -S >sam.P
less -S sam.P |awk 'OFS="\t"{print $1,$2"_"NR,$3,$4}' >sam.Padid
cat NPHP1SD.new.gaf >NPHP1SD.filt.gaf

# python   ./software/pggb/NPHP1/NPHP1/FINAL_GFA/sdpath.py


# while IFS= read -r sam; do
# python ./software/pggb/NPHP1/NPHP1/FINAL_GFA/sdpath_para.py --sam $x
# done <alsam

cat sam.Padid |cut -f2 | parallel -j 7 "python ./software/pggb/NPHP1/NPHP1/FINAL_GFA/sdpath_para.py --sam {}"
cd ./project/APG/NPHP1/02_PGGB/output/temp
cat *csv*  |grep  -v "SDregion">al.stru




cat NPHP1gene.gaf NPHP1SD.filt.gaf >al.gaf
import pandas as pd
import seaborn as sns
sdpa = pd.read_csv('al.gaf', sep='\t',header=None)
unique_groups = set(sdpa[0])
colors = sns.color_palette("husl", len(unique_groups))
hex_colors = [f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}" for r, g, b in colors]
group_color_dict = dict(zip(unique_groups, hex_colors))
sdpa['color'] = sdpa[0].map(group_color_dict)
df=sdpa.iloc[:,[5,17]]
df.columns=['segments','color']
def process_segments_df_clean(df):
    result_rows = []
    for idx, row in df.iterrows():
        segments_str = row['segments']
        color = row['color']
        segments_list = re.findall(r'[<>](\d+)', segments_str)
        for segment in segments_list:
            result_rows.append({'segment': segment, 'color': color})
    result_df = pd.DataFrame(result_rows)
    return result_df

import re
result_df = process_segments_df_clean(df)
result_df.to_csv("ABC.colorinit.csv",index=False,header=0,sep="\t")cd ./project/APG/NPHP1/03_haplo_sta/

python ./COMPLEXFEATURE/11_figure/fig2.case/gene_ratio.py --chr chr2 --s 110000000 --e 111200000




library(data.table)
library(stringr)
library(dplyr)

data<-fread("./project/APG/NPHP1/02_PGGB/output/temp/al.stru")

mapping <- c("chr2:110197933-110391878" = "A", 
             "chr2:110397814-110517542" = "B", 
             "chr2:110520233-110565558" = "C")


# 替换第一列
data$V1 <- sapply(data$V1 , function(x) {
  if (x %in% names(mapping)) {
    mapping[x]
  } else {
    x
  }
})

data <- data %>%
  mutate(sample = str_split_i(V4, "#", 1))


data$oristru=""
data <- data %>%
  mutate(oristru = ifelse(V5 == "-", tolower(V1), V1))

data$strutype="no"
for(samplename in unique(data$sample)){
  sada<-data[data$sample==samplename,]
  sada <- sada[order(sada$V2), ]
  samstr=paste(sada$oristru, collapse ="")
  data[data$sample==samplename,]$strutype=samstr
}



result <- data[, .(first_sample = first(sample)), by = strutype]
da=data[data$sample %in% result$first_sample,]
draw<-da[,c("strutype","V1","V2","V3","V5")]
colnames(draw)<-c("molecule","gene","start","end","strand")
target_haplotypes <- c("ABCbCa", "AcBCbCa", "AcBcbCa", "AcBcba", "AcbCbCa")
draw=draw[molecule %in% target_haplotypes, molecule := paste0("*", molecule)]



#colnames(draw)<-c("molecule","gene","start","end","strand")

# drawnow$molecule <- factor(draw$molecule, levels = unique(draw$molecule))
library(ggplot2)
library(gggenes)
library(RColorBrewer)
draw$direction <- ifelse(draw$strand == "+", 0,1)
drawnow<-draw
#drawnow<-draw
drawnow<-drawnow[,c("molecule","gene","start","end","direction")]
# pdf("D:/MS/SD/apg/sd/hprc.pdf",width = 10,height = 80)
drawnow$molecule <- sub("\\..*", "", drawnow$molecule)
drawnow$molecule <- factor(drawnow$molecule, levels = c("*ABCbCa", "ABcbCa", "*AcBCbCa", "ABCba", "*AcBcbCa", "*AcBcba","AcBCba", "ABcba", "*AcbCbCa", "AcBcBcba", "ABCa", "ABCbCbCa", "ABcBCa"))




out=distinct(data[,c("sample","strutype")])
SAMPOP<-distinct(fread("~/DATA/APG/sam_pop.csv"))
colnames(SAMPOP)[2]="sample"
end=merge(out,SAMPOP,by="sample")
risk=c("ABCbCa","AcBCbCa","AcBcbCa","AcBcba","AcbCbCa")

end$ifrisk=0
end[end$strutype %in% risk,]$ifrisk=1

nrow(end[end$SP=="EAS" & end$ifrisk==1,])
nrow(end[end$SP=="EAS" & end$ifrisk==0,])
nrow(end[end$SP!="EAS" & end$ifrisk==1,])
nrow(end[end$SP!="EAS" & end$ifrisk==0,])


nrow(end[end$SP=="AFR" & end$ifrisk==1,])
nrow(end[end$SP=="AFR" & end$ifrisk==0,])

nrow(end[end$SP=="AMR" & end$ifrisk==1,])
nrow(end[end$SP=="AMR" & end$ifrisk==0,])

nrow(end[end$SP=="SAS" & end$ifrisk==1,])
nrow(end[end$SP=="SAS" & end$ifrisk==0,])

data <-matrix(c(73,271,56,116), nrow = 2, byrow = TRUE)
colnames(data) <- c("Case", "Control")
rownames(data) <- c("EAS", "non-eas")
library(epitools)
oddsratio(data, rev="c")



target_haplotypes <- c("ABCbCa", "AcBCbCa", "AcBcbCa", "AcBcba", "AcbCbCa")
end=end[strutype %in% target_haplotypes, strutype := paste0("*", strutype)]


haplotype_stats <- end[, .(count = .N), by = .(strutype, SP)]
haplotype_stats[, total := sum(count), by = strutype]
haplotype_stats[, frequency := count / total]

