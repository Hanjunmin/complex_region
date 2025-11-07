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


data<-fread("del_dup.in",header=FALSE)
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

sub_df$total_length <- sapply(sub_df$name, function(x) {
  calculate_length(x)
})
sub_df <- sub_df[order(-sub_df$total_length), ]
sub_df$have=0
for(j in 1:nrow(sub_df)){
  print(j)
  if(sub_df[j,]$have==1){next}
  for(m in 1:nrow(sub_df)){
    if(j==m){next}
    if(sub_df[m,]$have==1){next}
cor1=sub_df[j,]$name
cor2=sub_df[m,]$name
A=calculate_overlap(cor1,cor2)
if(A>=0.9){print("x");print(j);print(m);sub_df[m,]$have=1}
}
}
alld<-rbind(alld,sub_df)
}
delreg<-alld[alld$have==1,]$name
writeLines(delreg, "del.reg")


