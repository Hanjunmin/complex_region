library(data.table)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
#edge<-fread("/home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/edge.bed",header=FALSE) ##"/home/jmhan/pangenome/WGSfea/HG00438T2TSDALL/edge.bed"
edge<-fread(args[1],header=FALSE) ##"/home/jmhan/pangenome/WGSfea/HG00438T2TSDALL/edge.bed"

edge<-distinct(edge)
colnames(edge)<-c("A","B")
# edge_sorted <- t(apply(edge, 1, sort))
# edge_unique <- edge_sorted[!duplicated(edge_sorted), ]
# edge_unique <-as.data.frame(edge_unique)
# colnames(edge_unique)<-c("A","B")

edge_sorted <- data.frame(
  A = pmin(edge$A, edge$B),
  B = pmax(edge$A, edge$B)
)

edge_unique <- edge_sorted[!duplicated(edge_sorted), ]


# t2t2sam<-fread("HG00438_Mat.sd")
t2t2sam<-fread(args[2])
t2t2sam<-t2t2sam[,c("SDregion","pc","sim")]
t2tsam2<-t2t2sam
colnames(t2t2sam)<-c("A","sam1","sim1")
colnames(t2tsam2)<-c("B","sam2","sim2")
edge_unique<-merge(edge_unique,t2t2sam,by="A",all.x=TRUE)
edge_unique<-merge(edge_unique,t2tsam2,by="B",all.x=TRUE)

edge_unique$sim1[is.na(edge_unique$sim1)] <- 0
edge_unique$sim2[is.na(edge_unique$sim2)] <- 0
edge_unique$sam1[is.na(edge_unique$sam1)] <- "nosam"
edge_unique$sam2[is.na(edge_unique$sam2)] <- "nosam"

process=edge_unique[edge_unique$sim1<0.97 | edge_unique$sim2<0.97,]
end=edge_unique[edge_unique$sim1>=0.97 & edge_unique$sim2>=0.97,]
fwrite(process,"process.txt",sep="\t")
fwrite(end,"end.txt",sep="\t")



# df_na_rows <- process[is.na(process[, 1]), ]
