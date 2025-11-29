library(data.table)
library(igraph)


parse_coord <- function(coord_str) {
  parts <- strsplit(coord_str, ":")[[1]]
  ran <- strsplit(parts[2], "-")[[1]]

  list(
    chr = parts[1], ##modify
    start = as.numeric(ran[1]),
    end = as.numeric(ran[2])
  )
}


calculate_length <- function(coord_str) {
  coord <- parse_coord(coord_str)
  return(coord$end - coord$start + 1)
}

calculate_overlap <- function(coord1, coord2) {
  coord1_parsed <- parse_coord(coord1)
  coord2_parsed <- parse_coord(coord2)
  if (coord1_parsed$chr != coord2_parsed$chr) {
    return(0)
  }
  overlap_start <- max(coord1_parsed$start, coord2_parsed$start)
  overlap_end <- min(coord1_parsed$end, coord2_parsed$end)
  if (overlap_start > overlap_end) {
    return(0)
  }
  
  overlap_length <- overlap_end - overlap_start + 1
  length1 <- coord1_parsed$end - coord1_parsed$start + 1
  length2 <- coord2_parsed$end - coord2_parsed$start + 1
  min_length <- min(length1, length2)
  return(overlap_length / min_length)
}



data<-fread("inter.reg.z",header=FALSE)
colnames(data) <- c("from", "to")
g <- graph_from_data_frame(data, directed = FALSE)
comp <- components(g)
subgraphs <- decompose.graph(g)
alld=c()
for(i in 1:length(subgraphs)) {
  print(i)
  subg <- subgraphs[[i]]
  sub_df <- as_data_frame(subg, what = "vertices")
  sub_df$id=i
  coords <- strsplit(sub_df$name, "@")
sub_df$coordA <- sapply(coords, function(x) x[1])
sub_df$coordB <- sapply(coords, function(x) x[2])



sub_df$total_length <- sapply(sub_df$name, function(x) {
  coords <- strsplit(x, "@")[[1]]
  calculate_length(coords[1]) + calculate_length(coords[2])
})
sub_df <- sub_df[order(-sub_df$total_length), ]
sub_df$have=0
for(j in 1:nrow(sub_df)){
  print(j)
  if(sub_df[j,]$have==1){next}
  for(m in 1:nrow(sub_df)){
    if(j==m){next}
    if(sub_df[m,]$have==1){next}
cor1=sub_df[j,]$coordA
cor2=sub_df[j,]$coordB
cor3=sub_df[m,]$coordA
cor4=sub_df[m,]$coordB
A=max(calculate_overlap(cor1,cor3),calculate_overlap(cor1,cor4))
B=max(calculate_overlap(cor2,cor3),calculate_overlap(cor2,cor4))
if(A>=0.9 & B >=0.9){print("x");print(j);print(m);sub_df[m,]$have=1}
}
}
alld<-rbind(alld,sub_df)
}
delreg<-alld[alld$have==1,c("coordA","coordB")]
savreg<-alld[alld$have==0,c("coordA","coordB")]

delreg$name=paste(delreg$coordA,delreg$coordB,sep="@")
delreg$name2=paste(delreg$coordB,delreg$coordA,sep="@")

savreg$name=paste(savreg$coordA,savreg$coordB,sep="@")
savreg$name2=paste(savreg$coordB,savreg$coordA,sep="@")
END<-fread("other.end")
END$name=paste(END$V1,":",END$V2,"-",END$V3,"@",END$V4,":",END$V5,"-",END$V6,sep="")
END <- as.data.frame(END)
ENDdelx <- END[!(END$name %in% c(delreg$name,delreg$name2)), ]
ENDsavx <- END[(END$name %in% c(savreg$name,savreg$name2)), ]
library(dplyr)
out=distinct(rbind(ENDdelx,ENDsavx))

out$normalized_name <- apply(out, 1, function(row) {
  region1 <- paste(row["V1"], row["V2"], row["V3"], sep = ":")
  region2 <- paste(row["V4"], row["V5"], row["V6"], sep = ":")
  regions <- sort(c(region1, region2))
  paste(regions, collapse = "@")
})

# 基于规范化名称去重，保留第一个出现的行
out_unique <- out[!duplicated(out$normalized_name), ]

# 移除临时列
out_unique$normalized_name <- NULL
out_unique$name <- NULL
fwrite(out_unique, "other.end.filt", sep = "\t", col.names = FALSE)