library(data.table)
library(igraph)
args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
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


library(dplyr)
# data<-fread("./t2tsdminimap/T2Thave.bed",header=FALSE)
data<-fread(infile,header=FALSE)



data<-distinct(data)
grouped <- data %>%
  group_by(across(13:14))%>%
  group_split()

alld=c()
for(i in seq_along(grouped)) {
  print(i)
  current_group <- grouped[[i]]
  sub_df <- as.data.frame(current_group)
  if(nrow(sub_df)==1){
      alld<-rbind(alld,sub_df)
  }else{


sub_df$coordA <- paste(sub_df$V1,":",sub_df$V2,"-",sub_df$V3,sep="")
sub_df$coordB <- paste(sub_df$V4,":",sub_df$V5,"-",sub_df$V6,sep="")

sub_df$total_length <- mapply(function(coordA, coordB) {
  calculate_length(coordA) + calculate_length(coordB)
}, sub_df$coordA, sub_df$coordB)

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
ad=sub_df[sub_df$have==0,]
ad$coordA <- NULL
ad$coordB <- NULL
ad$have <- NULL
ad$total_length<- NULL
alld<-rbind(alld,ad)
}}

outfile=args[2]
fwrite(alld, outfile, sep = "\t", col.names = FALSE)