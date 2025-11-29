library(data.table)
library(dplyr)
library(dbscan)
library(GenomicRanges)
args <- commandArgs(trailingOnly = TRUE)


##生成T2T的标准之后可以用这个整合，目的是把ref和query都有交集的看看符不符合那个条件
process_inte2 <- function(process_bef) {
  ranges_ref <- IRanges(start = process_bef$V2, end = process_bef$V3)
  ranges_query <- IRanges(start = process_bef$V5, end = process_bef$V6)
  process_bef$cluster="NO"
  gr_ref <- GRanges(seqnames = "chrA", ranges = ranges_ref)
  gr_query <- GRanges(seqnames = "chrA", ranges = ranges_query)
  overlaps1 <- DataFrame(findOverlaps(ranges_ref, ranges_ref))
  overlaps1 <-overlaps1[overlaps1$queryHits != overlaps1$subjectHits,]
  overlaps2 <- DataFrame(findOverlaps(ranges_query, ranges_query))
  overlaps2<-overlaps2[overlaps2$queryHits != overlaps2$subjectHits,]
  common_overlaps <- merge(overlaps1, overlaps2, by = c("queryHits", "subjectHits"))
  endall<-cbind(x[common_overlaps$queryHits,],x[common_overlaps$subjectHits,],common_overlaps$queryHits,common_overlaps$subjectHits)
  colnames(endall)[1:6]=c("rc1","rs1","re1","qc1","qs1","qe1")
  colnames(endall)[13:18]=c("rc2","rs2","re2","qc2","qs2","qe2")
  colnames(endall)[(ncol(endall)-1):ncol(endall)]=c("queryHits","subjectHits")
  endall$r1l=endall$re1-endall$rs1
  endall$r2l=endall$re2-endall$rs2
  endall$q1l=endall$qe1-endall$qs1
  endall$q2l=endall$qe2-endall$qs2
  endall$interrs <- apply(endall[, c("rs1", "rs2")], 1, max)
  endall$interre <- apply(endall[, c("re1", "re2")], 1, min)
  endall$interqs <- apply(endall[, c("qs1", "qs2")], 1, max)
  endall$interqe <- apply(endall[, c("qe1", "qe2")], 1, min)
  endall$interrl<-endall$interre-endall$interrs
  endall$interql<-endall$interqe-endall$interqs
  be=endall[,c("queryHits","subjectHits")]
  be$smaller <- pmin(be$queryHits, be$subjectHits)  # 较小的数
  be$larger <- pmax(be$queryHits, be$subjectHits)   # 较大的数
  df_sorted <- be[, c("smaller", "larger")]
  df_unique <- unique(df_sorted)
  library(igraph)
  g <- graph_from_data_frame(df_unique, directed = FALSE)
  components <- components(g)
  df_unique$cluster <- components$membership[as.numeric(V(g)[match(df_unique$smaller, V(g)$name)])]  # 根据 smaller 列查找
  df_unique$cluster[is.na(df_unique$cluster)] <- components$membership[as.numeric(V(g)[match(df_unique$larger, V(g)$name)])]  # 根据 larger 列查找
  for(clus in unique(df_unique$cluster)){
    ch=df_unique[df_unique$cluster==clus,]
    process_bef[unique(c(ch$smaller,ch$larger)),]$cluster=clus
  }
  return(process_bef)
  
}

# chromosome <- args[1]

# print(chromosome)
#  folder_path=paste("/home/jmhan/pangenome/SD/callSD/callSDnew0.05/group/chr",chromosome,sep="")
#  dir.create(folder_path, showWarnings = TRUE, recursive = TRUE)
#  system(paste("rm ",folder_path,"/*",sep=""))
#  system(paste("cd ",folder_path,sep=""))
#  system(paste("cat /home/jmhan/pangenome/SD/callSD/callSDnew0.05/chr",chromosome,".pair.txt |tr '_' '\t' >",folder_path,"/chr",chromosome,".pair.txt",sep=""))
#  system(paste("cat /home/jmhan/pangenome/SD/callSD/callSDnew0.05/chr",chromosome,".jac.txt |tr '_' '\t' >",folder_path,"/chr",chromosome,".jac.txt",sep=""))
# data<-fread(paste(folder_path,"/chr",chromosome,".pair.txt",sep=""))
# datajac<-fread(paste(folder_path,"/chr",chromosome,".jac.txt",sep=""))
# data$jac=datajac$V4
# data<-data[data$jac>=0.05,]
# split_data <- split(data, list(data$V1, data$V4))
library(data.table)
data<-fread(args[1])
#data<-fread("../endend")
#data=data[data$V1=="chr15" & data$V4=="chr22",]
split_data <- split(data, list(data$V1, data$V4))
allname<-names(split_data)
for(i in allname){
  if(nrow(split_data[[i]])==0){
    next
  }
  #dir.create(paste(folder_path,"/",names(split_data)[i],sep=""), showWarnings = TRUE, recursive = TRUE)
  #system(paste("rm ",folder_path,"/",names(split_data)[i],"/*",sep=""))
  x=split_data[[i]]
  x=distinct(x)
  x=process_inte2(x)
  x$combo <- apply(x[, c("V1", "V2", "V3", "V4", "V5", "V6")], 1, function(y) {
  # 将 V1 到 V6 列中的内容排序，确保位置互换的行被认为是重复的
  paste(sort(y), collapse = "-")
})
  x <- x[!duplicated(x$combo), ]
  x <- x %>% select(-combo)
  #fwrite(x, file = "chr15_chr22_cluster.csv")
  clustertable=table(x$cluster)
  # choose=names(clustertable[clustertable<=100 & clustertable>1])
  # nochoose=names(clustertable[clustertable>100])
  nointe="NO"
  chooseall=names(clustertable[clustertable>1])
  chooseall=chooseall[chooseall != "NO"]
  new=data.frame()
  otherlist=c()
  before=x[x$cluster %in% chooseall,]
  nointedata=x[x$cluster == nointe,]
  rc=unique(x$V1)
  qc=unique(x$V4)
  if (file.exists(paste(rc,qc,"before.csv",sep="_") )){
  print("a")
break
}else{
  fwrite(before, file = paste(rc,qc,"before.csv",sep="_") )
}
if (file.exists( paste(rc,qc,"nointe.csv",sep="_") )){
   print("c")
break
}else{
  fwrite(nointedata, file =paste(rc,qc,"nointe.csv",sep="_"))
}  


  if(length(chooseall)==0){
    next
  }
  for(clusid in chooseall){
    clusdata<-x[x$cluster==clusid,]
    rc=unique(clusdata[,1])
    nrs=min(clusdata[,2])
    nre=max(clusdata[,3])
    qc=unique(clusdata[,4])
    nqs=min(clusdata[,5])
    nqe=max(clusdata[,6])
    idd=clusid
    if(clustertable[clusid]>=100){
      large=1
    }else{
      large=0
    }
    str1=paste(rc,":",nrs,"-",nre,sep="")
    str2=paste(qc,":",nqs,"-",nqe,sep="")
    #sample="/home/jmhan/HPRC/HG00438.1.fa"
    sample=args[2]
    path="/home/jmhan/pangenome/SD/callSD/callSDnew0.05/cigar_call.py"
    path=args[3]
    command1=paste("bash -c 'minimap2 -c --eqx <(samtools faidx ",sample," ",str1,") <(samtools faidx ",sample," ",str2,") -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -o a.paf && python ",path," --paf a.paf --o a.txt'",sep="")
    #command1=paste("bash -c 'minimap2 -c --eqx <(samtools faidx ~/breakpoints/chm13v2.0.fa ",str1,") <(samtools faidx ~/breakpoints/chm13v2.0.fa ",str2,") -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -o a.paf && python /home/jmhan/pangenome/SD/callSD/callSDnew0.05/cigar_call.py --paf a.paf --o a.txt'",sep="")
    command <- "bash -c 'minimap2 -c --eqx <(samtools faidx ~/breakpoints/chm13v2.0.fa chr15:166487-222187) <(samtools faidx ~/breakpoints/chm13v2.0.fa chr22:226908-284237) -A 5 -B 4 -O 40 -E 1 -s 3000 -t 12 -o a.paf && python /home/jmhan/pangenome/SD/callSD/callSDnew0.05/cigar_call.py --paf a.paf --o a.txt'"
    system(command1,intern = FALSE)
    output1 <- system("cat a.txt|wc -l", intern = TRUE)
    if(output1=="1"){
      next
    }else{
      output <- system("cat a.txt", intern = TRUE)
      output=output[2:length(output)]
      if(length(output)==1){
        output <- paste0(output, "\n")

      }
      output <- paste(output, collapse = "\n")
    }
    newdata<-fread(input = output, sep = "\t", header = FALSE)
    #newdata<-read.table(text = output, sep = "\t", header = FALSE)
    if(nrow(newdata[newdata$V10<0.9 | newdata$V12<0.5,])!=0){
      otherlist<- c(otherlist,clusid)
    }
    gr <- GRanges(
    seqnames = clusdata$V1,
    ranges = IRanges(start = clusdata$V2, end = clusdata$V3)
    )
    len1=sum(width(reduce(gr)))
    gr <- GRanges(
    seqnames = newdata$V1,
    ranges = IRanges(start = newdata$V2, end = newdata$V3)
    )
    len2=sum(width(reduce(gr)))
    if(len1>len2){
      otherlist<- c(otherlist,clusid)
    }
    newdata$nrs=nrs
    newdata$nre=nre
    newdata$nqs=nqs
    newdata$nqe=nqe
    newdata$cluster=idd
    newdata$large=large
    new<-rbind(new,newdata)
    print(output)
    # newrow=as.data.frame(c(rc,nrs,nre,qc,nqs,nqe,idd))
    # colnames(newrow)=c("rc","nrs","nre","qc","nqs","nqe","cluster")
    # new=rbind(new,newrow)
    # colnames(new)=c("rc","nrs","nre","qc","nqs","nqe","cluster")
    }
after=new
other=setdiff(chooseall,unique(new$cluster))
other=c(otherlist,other)
if(length(other)!=0){
otherdata<-x[x$cluster %in% other,]
fwrite(otherdata, file = paste(rc,qc,"other.csv",sep="_"))
}

if (file.exists(paste(rc,qc,"after.csv",sep="_") )){
  print("b")
break
}else{
  fwrite(after, file = paste(rc,qc,"after.csv",sep="_"))
}


 }











library(dbscan)
cluster_intervals_df <- function(df, fragment_length = 10, eps = 20, minPts = 2) {
  split_intervals_equal_length <- function(start, end, num_segments) {
    seq(from = start, to = end, length.out = num_segments + 1)
  }
  get_length <- function(start, end) {
    return(abs(end - start))
  }
  colnames(df)[1:6]=c("refchr","ref_start","ref_end","quechr","query_start","query_end")
  df$ref_length <- mapply(get_length, df$ref_start, df$ref_end)
  df$query_length <- mapply(get_length, df$query_start, df$query_end)
  df$num_segments <- ceiling(pmin(df$ref_length, df$query_length) / 1000)  # 这里将每个片段的长度设为10，并向上取整
  df$ref_fragments <- mapply(split_intervals_equal_length, df$ref_start, df$ref_end, df$num_segments)
  df$query_fragments <- mapply(split_intervals_equal_length, df$query_start, df$query_end, df$num_segments)
  ref_fragments_flat <- unlist(df$ref_fragments)
  query_fragments_flat <- unlist(df$query_fragments)
  
  
  # 构建用于聚类的数据
  fragment_data <- data.frame(ref_fragments_flat, query_fragments_flat)
  
  # 使用 DBSCAN 进行聚类
  dbscan_result <- dbscan(fragment_data, eps = eps, minPts = minPts)
  dbscan_result <- dbscan(fragment_data, eps = 100000, minPts =1)
  # 将聚类结果添加到片段数据中
  fragment_data$cluster <- dbscan_result$cluster
  
  # 还原聚类结果：找到每个原始区间所属的主要聚类簇
  df$cluster <- sapply(1:nrow(df), function(i) {
    ref_frags <- df$ref_fragments[[i]]
    query_frags <- df$query_fragments[[i]]
    matched_rows <- fragment_data$ref_fragments_flat %in% ref_frags &
                    fragment_data$query_fragments_flat %in% query_frags
    major_cluster <- if (sum(matched_rows) > 0) {
      unique(fragment_data$cluster[matched_rows])[which.max(table(fragment_data$cluster[matched_rows]))]
    } else {
      NA  # 如果没有找到匹配的片段
    }
    return(major_cluster)
  })
  
  # 返回带有聚类结果的数据框
  return(df)
}
