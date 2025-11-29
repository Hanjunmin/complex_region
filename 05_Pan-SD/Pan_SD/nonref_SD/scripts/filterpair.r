library(data.table)
library(dplyr)
library(tidyr)
args <- commandArgs(trailingOnly = TRUE)
data1<-fread(args[1],header=FALSE)
data2<-fread(args[2],header=FALSE)
df_combined <- cbind(data1, data2)

colnames(df_combined)<-c("pair1","pair2","pair1.1","jac")
x=df_combined[df_combined$jac>=0.05,]
fwrite(x,args[3],sep="\t",col.names=FALSE)
