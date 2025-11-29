library(data.table)
library(ineq)
args <- commandArgs(trailingOnly=TRUE)
dotpath=args[1]
# data<-fread(paste(chr,".end.out",sep=""))
data<-fread("TEal.txt")

#covdf<-fread("/home/jmhan/project/APG/github_final/01_1p36/index_construction/dot_refine/allsampos")
covdf<-fread(dotpath)

for(sam in colnames(data)){
  print(sam)
  subdf=covdf[covdf$V4==sam,]
  if(nrow(subdf)!=0){
    fwrite(subdf,"x.txt",sep="\t",col.names=FALSE)
    system("bedtools coverage -a TE.head -b x.txt |awk '$7<0.9' |cut -f 1-3 >no")
    noreg=fread("no")
    noreglis=noreg$V2
    data[data$start %in% noreglis,c(sam)]=-1
  }
}



data[data == -1] <- NA
library(dplyr)
data1=data
sub_df <- data1[, -(1:3)]



trimmed_var <- function(x, prop = 0.05) {
  x=x[!is.na(x)]
  lower <- quantile(x, probs = prop, na.rm = TRUE)
  upper <- quantile(x, probs = 1 - prop, na.rm = TRUE)
  x_trimmed <- x[x >= lower & x <= upper]
  if (length(x_trimmed) < 2) {
    NA
  } else {
    sd(x_trimmed)
  }
}



trimmed_gini <- function(x, prop = 0.05) {
  x=x[!is.na(x)]
  lower <- quantile(x, probs = prop, na.rm = TRUE)
  upper <- quantile(x, probs = 1 - prop, na.rm = TRUE)
  x_trimmed <- x[x >= lower & x <= upper]
  if (length(x_trimmed) < 2) {
    NA
  } else {
    Gini(x_trimmed)
  }
}

trimsd<- apply(sub_df, 1, trimmed_var)
trimgini<- apply(sub_df, 1, trimmed_gini)

data =data %>%
  mutate(row_max = pmax(!!!select(., -c(1:3)), na.rm = TRUE))
data$trimsd=trimsd
data$trimgini=trimgini
fwrite(data[,c("chr","start","end","trimsd","trimgini","row_max")],"TEEND",sep="\t")

#x=unname(unlist(sub_df))