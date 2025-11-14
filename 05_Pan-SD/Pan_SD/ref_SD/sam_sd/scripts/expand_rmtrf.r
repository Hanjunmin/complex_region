library(dplyr)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)


df<-fread(args[1],header=FALSE) #"groupin"
df <- df %>% filter_all(all_vars(. != ""))
fai<-fread(args[2],header=FALSE) #"HG00438_Pat.fa.fai"
inte=args[4]
fai<-fai[,c(1,2)]
result <- df %>%
group_by(across(1:6)) %>%
summarize(
  min_col7 = min(V8, na.rm = TRUE),
  max_col8 = max(V9, na.rm = TRUE)
)
result$min_col7 <- pmin(result$min_col7, result$V2)
result$max_col8 <- pmax(result$max_col8, result$V3)
result$mi=result$V2-result$min_col7
result$ma=result$max_col8-result$V3
result$mi[result$mi < 0] <- 0
result$ma[result$ma < 0] <- 0
result$news=result$V5-as.numeric(inte)*result$mi
result$newe=result$V6+as.numeric(inte)*result$ma
end=result[,c("V4","news","newe")]
end$news[end$news < 0] <- 1
end$newe[end$newe < 0] <- 1
colnames(fai)[1]="V4"
end <- merge(end, fai, by = "V4", all.x = TRUE)
end$newe[end$newe > end$V2] <- end$V2[end$newe > end$V2]
fwrite(end[,c(1,2,3)],args[3],sep="\t",col.names=FALSE) ##"rm_in"