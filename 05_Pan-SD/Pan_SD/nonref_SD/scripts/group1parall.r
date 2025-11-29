# library(data.table)
# library(dplyr)
# library(dbscan)


# args <- commandArgs(trailingOnly = TRUE)
# folder_path=args[4]
# #"/home/jmhan/pangenome/WGSfea/sdseq/filt/"
# nameall=args[3]
#   data2<-fread(paste(folder_path,args[1],sep=""),header=FALSE,sep="\t")
#   if(nrow(data2)==0){
#     next
#   }
#   data1<-fread(paste(folder_path,args[2],sep=""),header=FALSE,sep="\t")
#   nochooseall<-data.frame()
#   nochoosedataall<-data.frame()
#   resultall<-data.frame()
#   for(j in 1:nrow(data1)){
#     print(j)
#     x=data1[j,]
#     if(nchar(x$V4)==1){
#     cho=as.integer(x$V4)}else{
#     cho <- ifelse(grepl(",", x$V4), as.integer(strsplit(x$V4, ",")[[1]]), as.integer(x$V4))
#   }
#     new=data2[cho,]
#     print(new)
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
#   fwrite(nochoosedataall,paste(folder_path,"/",nameall,"/",nameall,"_large.txt",sep=""), sep = "\t",col.names=FALSE)
#   fwrite(nochooseall,paste(folder_path,"/",nameall,"/",nameall,"_largecluster.txt",sep=""), sep = "\t",col.names=FALSE)
#   if(nrow(resultall)!=0){fwrite(resultall,paste(folder_path,"/",nameall,"/",nameall,"_cluster.bed",sep=""), sep = "\t",col.names=FALSE)}



# #nl test.txt | awk 'OFS="\t"{print $2,$3,$4,$5,$6,$7,$1}' |bedtools merge -i x -d 500000  -c 7 -o collapse |less -S >



library(data.table)
library(dplyr)
library(dbscan)


args <- commandArgs(trailingOnly = TRUE)
folder_path=args[4]
#"/home/jmhan/pangenome/WGSfea/sdseq/filt/"
nameall=args[3]
data2<-fread(paste(folder_path,args[1],sep=""),header=FALSE,sep="\t")
if(nrow(data2)==0){
  next
}
data1<-fread(paste(folder_path,args[2],sep=""),header=FALSE,sep="\t")
resultall_list=list()
nochooseall<-data.frame()
nochoosedataall<-data.frame()
for(j in 1:nrow(data1)){
    print(j)
    x=data1[j,]
    if(grepl(",", x$V4)){
      cho=as.integer(strsplit(x$V4, ",")[[1]])
    }else{
      cho=as.integer(x$V4)
    }
    new=data2[cho,]
    new$B_num <- as.numeric(factor(new$V4)) * 1e6  # 将 A 转换为较大的数字
    dbscan_result <- dbscan(new[,c(2,3,5,6,7)], eps = 300000, minPts =1)
    new$cluster=dbscan_result$cluster
    result <- new %>%
    group_by(cluster, V4) %>%  # 将 V4 直接用作列名
    summarize(
      ref_chr = first(V1),
      ref_min = min(V2),
      ref_max = max(V3),
      que_chr = first(V4),
      que_min = min(V5),
      que_max = max(V6),
      .groups = 'drop'  # 防止将结果的分组信息保留在结果中
    )
    result=as.data.frame(result)
    resultall_list[[paste(j, sep = "_")]] <- result
}

 final_result <- do.call(rbind,resultall_list)
 final_result=as.data.frame(final_result)
 final_result<- final_result[,c("ref_chr","ref_min","ref_max","que_chr","que_min","que_max")]
if(nrow(final_result)!=0){fwrite(final_result,paste(folder_path,"/",nameall,"/",nameall,"_cluster.bed",sep=""), sep = "\t",col.names=FALSE)}



#nl test.txt | awk 'OFS="\t"{print $2,$3,$4,$5,$6,$7,$1}' |bedtools merge -i x -d 500000  -c 7 -o collapse |less -S >









