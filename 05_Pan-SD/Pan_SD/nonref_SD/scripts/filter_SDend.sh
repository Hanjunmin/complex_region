###########################REF DEL
less -S deldup.init.x |tr ":" "\t" |tr "-" "\t" |tail -n +2 >x
less -S x |awk 'OFS="\t"{print $4,$5,$6,$1,$2,$3}' >y
bedtools intersect -a x -b x -wa -wb |awk '$4==$10' |awk '$5<$12 && $6>$11' |awk '!($2==$8 && $3==$9 && $5==$11 && $6==$12)' >inter.x
bedtools intersect -a y -b x -wa -wb |awk '$4==$10' |awk '$5<$12 && $6>$11' |awk '!($2==$8 && $3==$9 && $5==$11 && $6==$12)' >inter.y

cat inter.x |awk 'OFS="\t"{print $1":"$2"-"$3"@"$4":"$5"-"$6,$7":"$8"-"$9"@"$10":"$11"-"$12}' |sort |uniq >inter.reg.x
cat inter.y |awk 'OFS="\t"{print $4":"$5"-"$6"@"$1":"$2"-"$3,$7":"$8"-"$9"@"$10":"$11"-"$12}' |sort |uniq >inter.reg.y
cat inter.reg.x inter.reg.y |sort |uniq >inter.reg.z

library(data.table)
library(igraph)
parse_coord <- function(coord_str) {
  parts <- strsplit(coord_str, ":|-|:")[[1]]
  list(
    chr = parts[1],
    start = as.numeric(parts[2]),
    end = as.numeric(parts[3])
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
fwrite(delreg,"del.reg")


# ###########################other DEL
less -S deldup.init.y |cut -f4 |tr "@" "\t" |cut -f1|tail -n +2 >x
less -S deldup.init.y |cut -f4 |tr "@" "\t" |cut -f2 |tail -n +2>y
paste <(less -S x |tr ":" "\t" |cut -f1) <(less -S x |tr ":" "\t" |cut -f2 |tr "-" "\t") >x1
paste <(less -S y |tr ":" "\t" |cut -f1) <(less -S y |tr ":" "\t" |cut -f2 |tr "-" "\t") >y1
paste x1 y1 >x
less -S x |awk 'OFS="\t"{print $4,$5,$6,$1,$2,$3}' >y
bedtools intersect -a x -b x -wa -wb |awk '$4==$10' |awk '$5<$12 && $6>$11' |awk '!($2==$8 && $3==$9 && $5==$11 && $6==$12)' >inter.x
bedtools intersect -a y -b x -wa -wb |awk '$4==$10' |awk '$5<$12 && $6>$11' |awk '!($2==$8 && $3==$9 && $5==$11 && $6==$12)' >inter.y

cat inter.x |awk 'OFS="\t"{print $1":"$2"-"$3"@"$4":"$5"-"$6,$7":"$8"-"$9"@"$10":"$11"-"$12}' |sort |uniq >inter.reg.x
cat inter.y |awk 'OFS="\t"{print $4":"$5"-"$6"@"$1":"$2"-"$3,$7":"$8"-"$9"@"$10":"$11"-"$12}' |sort |uniq >inter.reg.y
cat inter.reg.x inter.reg.y |sort |uniq >inter.reg.z

library(data.table)
library(igraph)
parse_coord <- function(coord_str) {
  parts <- strsplit(coord_str, ":")[[1]][2]
  list(
    chr = strsplit(coord_str, ":")[[1]][1],
    start = as.numeric(strsplit(parts, "-")[[1]][1]),
    end = as.numeric(strsplit(parts, "-")[[1]][2])
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
fwrite(delreg,"delother.reg")
