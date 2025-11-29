library(data.table)
library(dplyr)
library(tidyr)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
#data<-fread("sdinter.txt",header=FALSE)
data<-fread(args[1],header=FALSE)
alldf<-data.frame()

ii=0
list=unique(data$V1)
# for(id in list){
	
 process_id=function(id){
 	print(id)
	x<-data[V1==id,]
	x<-distinct(x)
	setorder(x, V3)
combined_string <- paste(x$V5, x$V6, sep = "", collapse = "")
x <- x %>% 
  separate(V1, into = c("chromosome", "position"), sep = ":") %>% 
  separate(position, into = c("start", "end"), sep = "-")
x$poss=1
x$pose=as.numeric(x$V4)-as.numeric(x$V3)+1

x$start <- as.numeric(x$start)
x$end <- as.numeric(x$end)
x$V3 <- as.numeric(x$V3)
x$V4 <- as.numeric(x$V4)
setDT(x)
# Condition 1: x[m,2] > x[m,5] & x[m,3] >= x[m,6]
condition1=x[,start>V3 & end >=V4]
#condition1 <- x[, 2] > x[, 5] & x[, 3] >= x[, 6]
x[condition1, "poss"] <- x[condition1, 2] - x[condition1, 5] + 1
x[condition1, "pose"] <- x[condition1, 6] - x[condition1, 5] + 1

# Condition 2: x[m,2] <= x[m,5] & x[m,3] < x[m,6]
condition2=x[,start<=V3 & end <V4]
#condition2 <- x[, 2] <= x[, 5] & x[, 3] < x[, 6]
x[condition2, "poss"] <- 1
x[condition2, "pose"] <- x[condition2, 3] - x[condition2, 5] + 1

# Condition 3: x[m,2] > x[m,5] & x[m,3] < x[m,6]
condition3=x[,start>V3 & end <V4]
#condition3 <- x[, 2] > x[, 5] & x[, 3] < x[, 6]
x[condition3, "poss"] <- x[condition3, 2] - x[condition3, 5] + 1
x[condition3, "pose"] <- x[condition3, 3] - x[condition3, 5] + 1

difference <- x$pose - x$poss+1

# 求和
total_difference <- sum(difference)
as=x[1,]$poss
ae=x[1,]$poss+total_difference-1
new_row<-c(id,as,ae,combined_string)
return(new_row)
}

esults <- mclapply(list, process_id, mc.cores = 5)
for(res in esults){
	alldf<-rbind(alldf,res)
}
#fwrite(alldf,"alldf.txt",sep="\t")
fwrite(alldf,args[2],sep="\t",col.names = FALSE)


