library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
data<-fread(args[1])
data$label=0
id=1
for(m in 1:nrow(data)){
    rc=data[m,]$V1
    rs=data[m,]$V2
    re=data[m,]$V3
    qc=data[m,]$V4
    qs=data[m,]$V5
    qe=data[m,]$V6
    for(j in 1:nrow(data)){
    if(j==m){next}else{
    nrc=data[j,]$V1
    nrs=data[j,]$V2
    nre=data[j,]$V3
    nqc=data[j,]$V4
    nqs=data[j,]$V5
    nqe=data[j,]$V6
    if(nrc==rc & nqc==qc & nrs>=rs & nre <= re & nqs >= qs & nqe <= qe){
    data[j]$label=id
    data[m]$label=id
    }
    }
    }
    id=id+1
}
processno<-data[data$label==0,c(1,2,3,4,5,6)]
processinte<-data[data$label!=0,]


if(nrow(processinte)!=0){
    
result <- processinte %>%
  group_by(label) %>%
  summarise(
    nrc=V1,
    nrs = min(V2),
    nre = max(V3),
    nqc=V4,
    nqs = min(V5),
    nqe = max(V6),
  )%>% 
  ungroup() %>%    
  select(-label) 
result<-distinct(result)
colnames(processno)=colnames(result)
allres<-rbind(processno,result)

system(paste("mv", args[1], paste(args[1],".bef",sep="")))
fwrite(allres,args[1],sep="\t",col.names=FALSE)
}
