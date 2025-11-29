### less -S C001-CHA-E01-Mat.sd |tail -n +2 |awk 'OFS="\t"{print $2,$9":"$10+$3"-"$10+$4,$5}' >pos.dic


rev_mer <- function(df) {
  dfrev=df[,c(2,1,3)]
  colnames(dfrev)=colnames(df)
  dfall<-rbind(df,dfrev)
  return(dfall)
}

args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(dplyr)
sam=args[1]
# sam="C001-CHA-E01-Mat"
edgeall=args[2] ## /home/jmhan/project/APG/github_final/PanSD_endtest/run/03_refSD/sam_cor/1_SD_CANPATH/sd.edge
#edgeall="/home/jmhan/PANSDEND/01_REFSD/para/edge.bed"
data<-fread(edgeall,header=FALSE)
colnames(data)<-c("A","B")
data$name=paste(data$A,"@",data$B,sep="")

### 2 link file
mini="./t2tsdminimap/all.end.statistics.filt"
minidf<-fread(mini,header=FALSE)
sam2pos="T2T2SAM"
sam2posdf<-fread(sam2pos,header=FALSE)
if(nrow(sam2posdf)!=0){
colnames(sam2posdf)=c("c1","cs1","ce1","c2","cs2","ce2","chr1_init","chr2_init")
sam2posdf[, samlef := paste(c1, ":", cs1, "-", ce1, sep="")]
sam2posdf[, samrig := paste(c2, ":", cs2, "-", ce2, sep="")]
sam2posdf<-sam2posdf[,c("chr1_init","chr2_init","samlef","samrig")]
colnames(sam2posdf)=c("chr1_init","chr2_init","chr1_source","chr2_source")
sam2posdf=distinct(sam2posdf)
sam2posdf_sim<-merge(sam2posdf,sim,by=c("chr1_init","chr2_init"),all.x=TRUE)
}
###sim.txt
sim<-fread("sim.txt",header=FALSE)
if(nrow(sim)!=0){
sim<-sim[,c(1,2,4,6)]
colnames(sim)<-c("chr1_init","chr2_init","chr2_init.sim","chr1_init.sim")
sim=distinct(sim)

###提取重复的区域（initpair对应多个sim, 可能是这些pair对应的路径被打断为了多个）##有些区域尽管相似性为1，但是由于segment的dup导致对应了多个区域，这部分也作为dot
duplicate_rows <- duplicated(sim[, c("chr1_init", "chr2_init")]) | 
                  duplicated(sim[, c("chr1_init", "chr2_init")], fromLast = TRUE)
sim$chr2_init.sim[duplicate_rows] <- "."
sim$chr1_init.sim[duplicate_rows] <- "."
sim=distinct(sim)}




### 1.processnosam.txt (nosamdf)
nosam="processnosam.txt"
nosamdf<-fread(nosam,header=FALSE)
if(nrow(nosamdf)!=0){
colnames(nosamdf)<-c("chr1_init","chr2_init","chr2_source","chr2_init.sim","chr1_source","chr1_init.sim")
dot_columns <- c("chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel")
nosamdf[, (dot_columns) := "."]
nosamdf <- nosamdf %>%
  mutate(chr2_source = ifelse(chr2_source == "nosam", "not_align", "."))
nosamdf <- nosamdf %>%
  mutate(chr1_source = ifelse(chr1_source == "nosam", "not_align", "."))
}else{
  nosamdf=data.frame()
}


### 2.T2Thave.bed (t2thavedfinitpair)
t2thave="./t2tsdminimap/T2Thave.bed"
t2thavedf<-fread(t2thave,header=FALSE)
if(nrow(t2thavedf)!=0){
t2thavedf<-t2thavedf[,c(1:14)]



colnames(t2thavedf)<-c("chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel","chr2_source","chr1_source")
t2thavedf<-distinct(t2thavedf)
t2thavedfinitpair=merge(sam2posdf_sim,t2thavedf,by=c("chr1_source","chr2_source"))
}else{
  t2thavedfinitpair=data.frame()
}

### 3.T2Tpolymor1.bed
t2tpoly="./t2tsdminimap/T2Tpolymor1.bed"
t2tpolydf<-fread(t2tpoly,header=FALSE)
if(nrow(t2tpolydf)!=0){
t2tpolydf<-distinct(t2tpolydf[,c(1:14)])
colnames(t2tpolydf)<-c("chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel","chr2_source","chr1_source")
t2tpolydf<-distinct(t2tpolydf)
t2tpolydfinitpair=merge(sam2posdf_sim,t2tpolydf,by=c("chr1_source","chr2_source"))
}else{
  t2tpolydfinitpair=data.frame()
}


### 4.T2Tpolymor1.bed

rmrefiltpoly="./t2tsdminimap/filin.bedend"
rmrefiltpolydf<-fread(rmrefiltpoly,header=FALSE)
if(nrow(rmrefiltpolydf)!=0){
rmrefiltpolydfal=merge(rmrefiltpolydf[,c(1,2,3,4,5,6)],minidf,by=c("V1","V2","V3","V4","V5","V6"),all.x=TRUE)
rmrefiltpolydfal<-rmrefiltpolydfal[,c(1:14)]
colnames(rmrefiltpolydfal)<-c("chr1","start1","end1","chr2","start2","end2","match_base","mismatch_base","orient","identity","align_len","matchindel","chr2_source","chr1_source")
rmrefiltpolydfal<-distinct(rmrefiltpolydfal)
rmrefiltpolydfinitpair=merge(sam2posdf_sim,rmrefiltpolydfal,by=c("chr1_source","chr2_source"))}else{
   rmrefiltpolydfinitpair=data.frame()
}



END=rbind(t2thavedfinitpair,t2tpolydfinitpair,rmrefiltpolydfinitpair,nosamdf)
END$name=paste(END$chr1_init,"@",END$chr2_init,sep="")
END$name1=paste(END$chr2_init,"@",END$chr1_init,sep="")

othername=data$name[!(data$name %in% c(END$name,END$name1))]
END$name<-NULL
END$name1<-NULL


if(length(othername)>0){
  library(tidyr)
noSD_df <- data.frame(othername = othername) %>%
  separate(othername, into = c("col1", "col2"), sep = "@")

colnames(noSD_df)<-c("chr1_init","chr2_init")
##
anchor=fread("pos.dic")
duplicate_rows <- duplicated(anchor[, c("V1")]) | 
                  duplicated(anchor[, c("V1")], fromLast = TRUE)
anchor$V2[duplicate_rows] <- "."
anchor$V3[duplicate_rows] <- "."
anchor<-distinct(anchor)

sim=distinct(sim)

colnames(anchor)<-c("chr1_init","chr1_source","chr1_init.sim")
mer1=merge(noSD_df,anchor,by=c("chr1_init"),all.x=TRUE)
colnames(anchor)<-c("chr2_init","chr2_source","chr2_init.sim")
mer2=merge(mer1,anchor,by=c("chr2_init"),all.x=TRUE)
mer2 <- as.data.table(mer2)
dot_columns <- c("chr1","start1","end1","chr2","start2","end2",
                 "match_base","mismatch_base","orient","identity",
                 "align_len","matchindel")
mer2[, (dot_columns) := "."]
mer2[, (names(mer2)) := lapply(.SD, function(x) {
  ifelse(is.na(x), ".", x)
})]
DFAL=rbind(END,mer2)
DFAL$SAM=sam
}else{
DFAL=END
DFAL$SAM=sam
}





fwrite(DFAL,paste(sam,".ref.SD.al",sep=""),sep="\t")
