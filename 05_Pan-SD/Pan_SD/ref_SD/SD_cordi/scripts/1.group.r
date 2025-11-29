library(data.table)
library(dplyr)
# library(dbscan)
#folder_path="/home/jmhan/pangenome/SD/callSD/panSD/filt/split/"
#data<-fread("/home/jmhan/pangenome/SD/callSD/panSD/filt/beforegroup.txt")
args <- commandArgs(trailingOnly = TRUE)

folder_path=args[1]  ##/home/jmhan/pangenome/WGSfea/sdseq/filt/
data<-fread( paste(args[2],"beforegroup.txt",sep="")) ##"/home/jmhan/pangenome/WGSfea/sdseq/filt/beforegroup.txt"

split_data <- split(data, list(data$V1, data$V4))
### test
for(i in 1:length(split_data)){
  print(i)
  x=split_data[[i]]
  if(nrow(x)==0){
    next
  }
  dir.create(paste(folder_path,"/",names(split_data)[i],sep=""), showWarnings = TRUE, recursive = TRUE)
  system(paste(
    "if [ -d ", folder_path, "/", names(split_data)[i], " ]; then rm -f ", folder_path, "/", names(split_data)[i], "/*; fi", sep=""
  )) 
  fwrite(x,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],".txt",sep=""),sep="\t",col.names=FALSE)
  proname=paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],".txt",sep="")
  comma0<-paste("bedtools subtract -a ",proname," -b ", args[3]," -A >",folder_path,"/",names(split_data)[i],"/a",sep="")
  #system(comma0) 
  comma1<-paste("bedtools subtract -a <(cat ",folder_path,"/",names(split_data)[i],"/a","|awk 'OFS=\"\t\"{print $4,$5,$6,$1,$2,$3}') -b ", args[3]," -A >",folder_path,"/",names(split_data)[i],"/end.txt",sep="")
  #system(paste("bash -c", shQuote(comma1)))
  comma2<-paste("cat " ,folder_path,"/",names(split_data)[i],"/end.txt"," |awk 'OFS=\"\t\"{print $4,$5,$6,$1,$2,$3}' >",folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"end.txt",sep="")
  #system(comma2) 
pronamen<-paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],".txt",sep="")
  command <- paste(
    "nl ",pronamen,"|",
    "awk 'OFS=\"\\t\"{print $2,$3,$4,$5,$6,$7,$1}' |bedtools sort |",
    "bedtools merge -i - -d 10000 -c 7 -o collapse |",
    "less -S >",folder_path,"/",names(split_data)[i],"/",
    names(split_data)[i],".clu1"
  ,sep="")
  system(command)

}

# for(i in 1:length(split_data)){
#   data2<-split_data[[i]]
#   print(i)
#   if(nrow(data2)==0){
#     next
#   }
#   data1<-fread(paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],".clu1",sep=""))
#   nochooseall<-data.frame()
#   nochoosedataall<-data.frame()
#   resultall<-data.frame()
#   for(j in 1:nrow(data1)){
#     print(j)
#     x=data1[j,]
#     cho=strsplit(x$V4,",")
#     cho=as.integer(cho[[1]])
#     new=data2[cho,]
#     dbscan_result <- dbscan(new[,c(2,3,5,6)], eps = 300000, minPts =1)
#     new$cluster=dbscan_result$cluster
#     clustertable=table(dbscan_result$cluster)
#     choose=names(clustertable[clustertable<=1000])
#     nochoose=names(clustertable[clustertable>1000])
#     new_cluster_map <- setNames(seq_along(nochoose),   nochoose)
#     nochoosedata=new[new$cluster %in% nochoose,]
#     nochoosedata[, new_cluster := new_cluster_map[as.character(cluster)]]
    
#   for(clusterid in unique(nochoosedata$cluster)){
#     nocd=nochoosedata[nochoosedata$cluster==clusterid,]
#     dbscan_result <- dbscan(nocd[,c(2,3,5,6)], eps = 100000, minPts =1)
#     nocd$cluster=dbscan_result$cluster
#     result <- nocd %>%
#     group_by(cluster) %>%
#     summarize(
#       ref_chr= first(V1),
#       ref_min = min(V2),
#       ref_max = max(V3),
#       que_chr= first(V4),
#       que_min = min(V5),
#       que_max = max(V6),
#       .groups = 'drop'  # 防止将结果的分组信息保留在结果中
#     )
#     result=as.data.frame(result)
#     result<-result[,c(2,3,4,5,6,7)]
#     nochooseall<-rbind(nochooseall,result)
#     nochoosedataall<-rbind(nochoosedataall,nochoosedata)
#   }
#   alldata<-new[new$cluster %in% choose,]
#   result <- alldata %>%
#   group_by(cluster) %>%
#   summarize(
#     ref_chr= first(V1),
#     ref_min = min(V2),
#     ref_max = max(V3),
#     que_chr= first(V4),
#     que_min = min(V5),
#     que_max = max(V6),
#     .groups = 'drop'  # 防止将结果的分组信息保留在结果中
#   )
#   result=as.data.frame(result)
#   result<-result[,c(2,3,4,5,6,7)]
#   resultall<-rbind(resultall,result)

#   }
#   fwrite(nochoosedataall,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_large.txt",sep=""), sep = "\t",col.names=FALSE)
#   fwrite(nochooseall,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_largecluster.txt",sep=""), sep = "\t",col.names=FALSE)
#   if(nrow(resultall)!=0){fwrite(resultall,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_cluster.bed",sep=""), sep = "\t",col.names=FALSE)}

# }


# #nl test.txt | awk 'OFS="\t"{print $2,$3,$4,$5,$6,$7,$1}' |bedtools merge -i x -d 500000  -c 7 -o collapse |less -S >






# for(i in 1:nrow(data1)){
  
  
#   dbscan_result <- dbscan(new[,c(2,3,5,6)], eps = 300000, minPts =1)

# }



# for(i in 314:length(split_data)){
#   x=split_data[[i]]
#   if(nrow(x)==0){
#     next
#   }
#   print(nrow(x))
#   dir.create(paste(folder_path,"/",names(split_data)[i],sep=""), showWarnings = TRUE, recursive = TRUE)
#   system(paste("rm ",folder_path,"/",names(split_data)[i],"/*",sep=""))
#   dbscan_result <- dbscan(x[,c(2,3,5,6)], eps = 300000, minPts =1)
#   x$cluster=dbscan_result$cluster
#   clustertable=table(dbscan_result$cluster)
#   choose=names(clustertable[clustertable<=1000])
#   nochoose=names(clustertable[clustertable>1000])
#   new_cluster_map <- setNames(seq_along(nochoose),   nochoose)
#   nochoosedata=x[x$cluster %in% nochoose,]
#   nochoosedata[, new_cluster := new_cluster_map[as.character(cluster)]]
#   nochooseall<-data.frame()
#   for(clusterid in unique(nochoosedata$cluster)){
#     nocd=nochoosedata[nochoosedata$cluster==clusterid,]
#     dbscan_result <- dbscan(nocd[,c(2,3,5,6)], eps = 100000, minPts =1)
#     nocd$cluster=dbscan_result$cluster
#     result <- nocd %>%
#     group_by(cluster) %>%
#     summarize(
#       ref_chr= first(V1),
#       ref_min = min(V2),
#       ref_max = max(V3),
#       que_chr= first(V4),
#       que_min = min(V5),
#       que_max = max(V6),
#       .groups = 'drop'  # 防止将结果的分组信息保留在结果中
#     )
#     result=as.data.frame(result)
#     result<-result[,c(2,3,4,5,6,7)]
#     nochooseall<-rbind(nochooseall,result)
#   }
#   fwrite(nochoosedata,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_large.txt",sep=""), sep = "\t",col.names=FALSE)
#   fwrite(nochooseall,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_largecluster.txt",sep=""), sep = "\t",col.names=FALSE)
  

#   alldata<-x[x$cluster %in% choose,]
#   result <- alldata %>%
#   group_by(cluster) %>%
#   summarize(
#     ref_chr= first(V1),
#     ref_min = min(V2),
#     ref_max = max(V3),
#     que_chr= first(V4),
#     que_min = min(V5),
#     que_max = max(V6),
#     .groups = 'drop'  # 防止将结果的分组信息保留在结果中
#   )
#   result=as.data.frame(result)
#   result<-result[,c(2,3,4,5,6,7)]
#   if(nrow(result)!=0){fwrite(result,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_cluster.bed",sep=""), sep = "\t",col.names=FALSE)}
# }


  # for(clusid in choose){
  #   print(clusid)
  #   clusdata<-x[x$cluster==clusid,]
  #   # clusdata=clusdata[(clusdata$V5<clusdata$V3 & clusdata$V5>clusdata$V2 & clusdata$jac>0.5) | (clusdata$V5>=clusdata$V3 |clusdata$V5<=clusdata$V2),]
  #   # clusdata=clusdata[clusdata$jac>=0.05,]
  #   clusdata <- clusdata %>%
  #               arrange(V2, V5)
  #   clusdata$V2=clusdata$V2-1
  #   clusdata$V5=clusdata$V5-1
  #   clusdata$lab=0
  #   df <- data.frame(A = numeric(0), B = numeric(0), C = numeric(0), D = numeric(0))
  #     for (row in 1:nrow(clusdata)){
  #       # print(row)
  #       if(nrow(clusdata)==1 |nrow(clusdata)==0){break}
  #       if(row==nrow(clusdata)){next}
  #       if(clusdata[row,]$lab==1){next}
  #         list=c(row)
  #         rs=clusdata[row,]$V2
  #           re=clusdata[row,]$V3
  #           qs=clusdata[row,]$V5
  #           qe=clusdata[row,]$V6
  #         for (id in (row+1):nrow(clusdata)){
  #           if(clusdata[id,]$lab==1){next}

  #           if((clusdata[id,]$V2>=rs-1000 & clusdata[id,]$V2<=re+1000 & clusdata[id,]$V5>=qs-1000 & clusdata[id,]$V5<=qe+1000)){
  #             rs=min(rs,clusdata[id,]$V2)
  #             re=max(re,clusdata[id,]$V3)
  #             # print(re)
  #             qs=min(qs,clusdata[id,]$V5)
  #             qe=max(qe,clusdata[id,]$V6)
  #             list=append(list,id)
  #             clusdata[id,]$lab=1
  #           }
  #         }
  #         # print(re)
  #       df <- rbind(df, as.data.frame(t(c(rs,re,qs,qe))))
  #   }
    
  #   if(nrow(clusdata)==1){
  #     df=clusdata[,c(1,2,3,4,5,6)]
  #   }else{
  #     df$chr1<-unique(clusdata$V1)
  #   df$chr2<-unique(clusdata$V4)
  #   df<-df[,c(5,1,2,6,3,4)]
  #   }
  #   alldata<-rbind(alldata,df)
  #   #if(nrow(df)!=0){fwrite(df,paste(folder_path,"/",names(split_data)[i],"/",names(split_data)[i],"_",clusid,".bed",sep=""), sep = "\t",col.names=FALSE)}
  # }
    
  # }







# args <- commandArgs(trailingOnly = TRUE)
# chromosome <- args[1]

# print(chromosome)
#  folder_path=paste("/home/jmhan/pangenome/SD/callSD/callSDnew0.05/group/chr",chromosome,sep="")
#  dir.create(folder_path, showWarnings = TRUE, recursive = TRUE)
#  system(paste("rm ",folder_path,"/*",sep=""))
#  system(paste("cd ",folder_path,sep=""))
#  system(paste("cat /home/jmhan/pangenome/SD/callSD/callSDnew0.05/chr",chromosome,".pair.txt |tr '_' '\t' >",folder_path,"/chr",chromosome,".pair.txt",sep=""))
#  system(paste("cat /home/jmhan/pangenome/SD/callSD/callSDnew0.05/chr",chromosome,".jac.txt |tr '_' '\t' >",folder_path,"/chr",chromosome,".jac.txt",sep=""))

# datajac<-fread(paste(folder_path,"/chr",chromosome,".jac.txt",sep=""))
# data$jac=datajac$V4
# data<-data[data$jac>=0.05,]
# split_data <- split(data, list(data$V1, data$V4))


