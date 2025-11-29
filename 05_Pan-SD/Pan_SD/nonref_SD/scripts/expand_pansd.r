library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
data<-fread(args[1],header=FALSE,sep="\t") ##"lefin"
fai<-fread(args[2]) ##"/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/process/HG00438_Pat/t2tsdminimap/HG00438_Pat.fa.fai"
fai<-fai[,c(1,2)]
colnames(fai)<-c("pairchr","pairlen")
result <- data %>%
group_by(V1,V2,V3,V4,V5,V6,V7,V8) %>%  # 将 V4 直接用作列名
summarize(
  rmtrf_s =min(V10),
  rmtrf_e = max(V11),
)
result=as.data.frame(result)
colnames(result)<-c("t2tchr","t2ts","t2te","pairchr","pairs","paire","cans","cane","rmtrf_s","rmtrf_e")
result$rmtrfran=result$rmtrf_e-result$t2te
result$rmtrfran2=result$t2ts-result$rmtrf_s
result$rmtrfran[result$rmtrfran< 0] <- 0
result$rmtrfran2[result$rmtrfran2 < 0] <- 0
result$pairsn<-result$pairs-result$rmtrfran2
result$pairen<-result$paire+result$rmtrfran
end <- merge(result, fai, by = "pairchr", all.x = TRUE)
end$pairen[end$pairen > end$pairlen] <- end$pairlen[end$pairen > end$pairlen]
end$pairsn[end$pairsn < 0] <- 0
expand=end[,c("pairchr","pairs","paire","pairsn","pairen")]
fwrite(expand, args[3],sep = "\t",col.names=FALSE)