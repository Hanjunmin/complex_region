library(data.table)
args <- commandArgs(trailingOnly = TRUE)
data<-fread(args[1])
data2<-fread(args[2])

#data<-data[!(data$V6 >= data$V2 & data$V7 <= data$V3), ]
#data2<-data2[!(data2$V6 >= data2$V2 & data2$V7 <= data2$V3), ]
library(dplyr)
if(nrow(data)!=0){
  data <- data %>% mutate(row_id = row_number())
  result <- data %>%
    group_by(V1, V2, V3,V7) %>%
    summarise(
  id=V7,
      refchr = V8,
      ref_min= min(V2,V9),
      ref_max = max(V3,V10),
      row_id = first(row_id)  
    ) %>%
    ungroup() %>%
    arrange(row_id) %>%
    select(-row_id)
  left=as.data.frame(distinct(result))
}

if(nrow(data2)!=0){
  data2 <- data2 %>% mutate(row_id = row_number())
  result <- data2 %>%
    group_by(V1, V2, V3,V4) %>%
    summarise(
  id=V4,
     quechr =V5,
     que_min = min(V2,V6),
     que_max = max(V3,V7),
      row_id = first(row_id)
    ) %>%
    ungroup() %>%
    arrange(row_id) %>%
    select(-row_id) 
  right=as.data.frame(distinct(result))
}


datano<-fread(args[4])
data2no<-fread(args[5])

if(nrow(data2no)!=0){
  colnames(data2no)<-c("quechr","que_min" ,"que_max" ,"id" )
}

if(nrow(data2)==0){
  rightall<-data2no

}else{
  rightall<-rbind(right[,c("id","quechr","que_min","que_max")],data2no)

}


if(nrow(data)!=0){
if(nrow(datano)!=0){
  colnames(datano)<-c("refchr" ,"ref_min" ,"ref_max"  ,"quechr" ,"quemin","queend","id" )
  leftall<-rbind(left[,c("id", "refchr" ,"ref_min"  ,"ref_max")],datano[,c("id", "refchr" ,"ref_min"  ,"ref_max")])
}else{
  leftall<-rbind(left[,c("id", "refchr" ,"ref_min"  ,"ref_max")],datano)
}
}else{
  colnames(datano)<-c("refchr" ,"ref_min" ,"ref_max"  ,"quechr" ,"quemin","queend","id" )
  leftall<-datano[,c("id", "refchr" ,"ref_min"  ,"ref_max")]
}

print(leftall)
print("a")
print(rightall)

all<-merge(leftall,rightall,by="id")

fwrite(all,args[3], sep = "\t",col.names=FALSE)
