library(data.table)
library(dplyr)
data<-fread("rm_inter_higncov")
if(nrow(data) == 0){
  df <- data.frame()
  fwrite(df, "ALL.HIGHCOV.RMTRF", sep = "\t", col.names = FALSE)
  quit(save = "no", status = 0)  # 正确的退出命令
}

colnames(data)<-c("rmtrfCHR","rmtrfs","rmtrfe","highcovrefc","highcovrefs","highcovrefe","highcovsamc","highcovsams","highcovsame")
data$sams=0
data$same=0
data$lefs=(data$rmtrfs-data$highcovrefs)/(data$highcovrefe-data$highcovrefs)
data$lefe=(data$highcovrefe-data$rmtrfe)/(data$highcovrefe-data$highcovrefs)
data$lefns=ceiling(data$highcovsams+(data$highcovsame-data$highcovsams)*data$lefs)
data$lefne=ceiling(data$highcovsame-(data$highcovsame-data$highcovsams)*data$lefe)
split_data <- data %>%
  group_by(highcovrefc,highcovrefs,highcovrefe) %>%
  group_split()
alldf<-data.frame()
id<-c()
del<-c()
for (i in 1:length(split_data)) {
  x=split_data[[i]]
  a=min(x$lefs)
  b=min(x$lefe) 
  if((a<0 & abs(a)>0.1)|(b<0 & abs(b)>0.1)){del<-c(del,i)
}else{id<-c(id,i)}
}
selected_groups <- split_data[id]
save_data <- bind_rows(selected_groups)
del_groups <- split_data[del]
save1=data.frame()
for (i in 1:length(del_groups)) {
  if(length(del_groups)==0){break}
  x=del_groups[[i]]
  refi=x[x$lefs>0 & x$lefe>0,]
  if(nrow(refi)!=0){
    inits=x[1,]$highcovrefs
    inite=x[1,]$highcovrefe
    initsams=x[1,]$highcovsams
    initsame=x[1,]$highcovsame
    if(min(x$lefs)<0){news=min(refi$rmtrfs);ratio=(news-inits)/(inite-inits);samnews=ratio*(initsame-initsams)+initsams;refi$highcovsams=ceiling(samnews);refi$highcovrefs=news}
    if(min(x$lefe)<0){newe=max(refi$rmtrfe);ratio=(inite-newe)/(inite-inits);samnewe=initsame-ratio*(initsame-initsams);refi$highcovsame=ceiling(samnewe);refi$highcovrefe=newe}
if(min(refi$highcovsame)-max(refi$highcovsams)<0){print("NO")}else{save1<-rbind(save1,refi)}
if(nrow(refi[refi$lefns-refi$lefne>0,])!=0){print("NO!!");break}
  }	
}
end=rbind(save_data,save1)
df=end[,c("highcovsamc","highcovsams","highcovsame","highcovsamc","lefns","lefne")]
options(scipen = 999)  # 禁用科学计数法
df$lefns[df$lefns < 0] <- 1
fwrite(df,"ALL.HIGHCOV.RMTRF",sep="\t",col.names=FALSE) ##这里左边为样本的区域A 右边为样本该区域下的repeatsB
