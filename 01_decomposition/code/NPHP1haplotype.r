##data<-data[data$oristrutype %in% setdiff(unique(data$oristrutype),strings)]
##data<-data[order(data$oristrutype)]



library(data.table)
library(dplyr)
apgsam<-fread("D:/MS/SD/apgsample.CSV",header=FALSE)
apg<-fread("D:/MS/SD//apg/sd/end",header=FALSE)
hprc<-fread("D:/MS/SD/end_hprc_nphp1",header=FALSE)
hgsvc<-fread("D:/MS/SD/hgsvc.sd.bed",header=FALSE)

hprcsam<-fread("D:/MS/APG/hprc_year1_sample_metadata.txt")
hgsvcsam<-fread("D:/MS/SD/decompose/HGSVC/hgsvcsam.CSV")


sdlist=c("chr2:110197933-110391878","chr2:110397814-110517542","chr2:110520233-110565558","chr2:110520233-110565558")
apg<-apg[apg$V1 %in% sdlist]
hprc<-hprc[hprc$V1 %in% sdlist]
hgsvc<-hgsvc[hgsvc$V1 %in% sdlist]
apg$sample <- apply(apg, 1, function(row) strsplit(row["V6"], ":")[[1]][1])
hprc$sample <- apply(hprc, 1, function(row) strsplit(row["V6"], ":")[[1]][1])
hgsvc$sample <- apply(hgsvc, 1, function(row) strsplit(row["V6"], ":")[[1]][1])
## full srtructure
mapping <- c("chr2:110197933-110391878" = "A", 
             "chr2:110397814-110517542" = "B", 
             "chr2:110520233-110565558" = "C")

# 替换第一列
hprc$V1 <- sapply(hprc$V1 , function(x) {
  if (x %in% names(mapping)) {
    mapping[x]
  } else {
    x
  }
})

apg$V1 <- sapply(apg$V1 , function(x) {
  if (x %in% names(mapping)) {
    mapping[x]
  } else {
    x
  }
})
hgsvc$V1 <- sapply(hgsvc$V1 , function(x) {
  if (x %in% names(mapping)) {
    mapping[x]
  } else {
    x
  }
})
hprc$oristru=""
hprc <- hprc %>%
mutate(oristru = ifelse(V5 == "-", tolower(V1), V1))
samdir=list()
samdirori=list()
hprc$strutype="no"
hprc$oristrutype="no"
for(samplename in unique(hprc$sample)){
  sada<-hprc[hprc$sample==samplename,]
  sada <- sada[order(sada$V8), ]
  samstr=paste(sada$V1, collapse ="")
  samstror=paste(sada$oristru, collapse ="")
  hprc[hprc$sample==samplename,]$strutype=samstr
  hprc[hprc$sample==samplename,]$oristrutype=samstror
  #samdir[[samplename]]=samstr
  #samdirori[[samplename]]=samstror
}
hprc$sampleini <- apply(hprc, 1, function(row) strsplit(row["sample"], "#")[[1]][1])


apg$oristru=""
apg <- apg %>%
  mutate(oristru = ifelse(V5 == "-", tolower(V1), V1))
samdir=list()
samdirori=list()
apg$strutype="no"
apg$oristrutype="no"
for(samplename in unique(apg$sample)){
  sada<-apg[apg$sample==samplename,]
  sada <- sada[order(sada$V8), ]
  samstr=paste(sada$V1, collapse ="")
  samstror=paste(sada$oristru, collapse ="")
  apg[apg$sample==samplename,]$strutype=samstr
  apg[apg$sample==samplename,]$oristrutype=samstror
  #samdir[[samplename]]=samstr
  #samdirori[[samplename]]=samstror
}
apg$sampleini <- apply(apg, 1, function(row) strsplit(row["sample"], "#")[[1]][1])


hgsvc$oristru=""
hgsvc <- hgsvc %>%
  mutate(oristru = ifelse(V5 == "-", tolower(V1), V1))
samdir=list()
samdirori=list()
hgsvc$strutype="no"
hgsvc$oristrutype="no"
for(samplename in unique(hgsvc$sample)){
  sada<-hgsvc[hgsvc$sample==samplename,]
  sada <- sada[order(sada$V8), ]
  samstr=paste(sada$V1, collapse ="")
  samstror=paste(sada$oristru, collapse ="")
  hgsvc[hgsvc$sample==samplename,]$strutype=samstr
  hgsvc[hgsvc$sample==samplename,]$oristrutype=samstror
  #samdir[[samplename]]=samstr
  #samdirori[[samplename]]=samstror
}
hgsvc$sampleini <- apply(hgsvc, 1, function(row) strsplit(row["sample"], "_")[[1]][1])

hprchap<-unique(hprc$oristrutype)
apghap<-unique(apg$oristrutype)
hgsvchap<-unique(hgsvc$oristrutype)




apgsam<-apgsam[,c(2,4)]
colnames(apgsam)<-c("population","sampleini")
hprcsam<-hprcsam[,c(1,6)]
colnames(hprcsam)<-c("sampleini","population")
apgall<-merge(apg,apgsam,by="sampleini")
hprcall<-merge(hprc,hprcsam,by="sampleini")
hgsvcsam<-hgsvcsam[,c(1,6)]
colnames(hgsvcsam)<-c("sampleini","population")
hgsvcall<-merge(hgsvc,hgsvcsam,by="sampleini")


dataall<-apgall[,c(11,14,15)]
dataall<-distinct(dataall)
dataall2<-hprcall[,c(11,14,15)]
dataall2<-distinct(dataall2)
dataall3<-hgsvcall[,c(11,14,15)]
dataall3<-distinct(dataall3)

dataall2 <- subset(dataall2, !is.na(population))
#dataall <- subset(dataall, !isEmpty(population))
dataall$population<-"EAS(APG)"
dataall<-rbind(dataall,dataall2,dataall3)
dataall<-dataall[dataall$population!="",]
group_counts <- prop.table(table(dataall$population, dataall$oristrutype), margin = 1)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
# 将数据转换为长格式
group_counts_df <- as.data.frame.matrix(group_counts)
group_counts_df$Population <- rownames(group_counts_df)

# 将数据转换为长格式
group_counts_melt <- melt(group_counts_df, id.vars = "Population", variable.name = "Group", value.name = "Proportion")
# 使用 'Set3' 调色板，它支持生成 12 种颜色
colors <- brewer.pal(13, "Set3")
colors <-c("#8DD3C7","#F2FFAA","#FFF055","#FFE3C9","#DFCDBE","#CD7100","#FFE26A","#E94400","#EFEBFF","#945B00","#7F1A00","#AA2200","#FFCD94")
colors <-c(colors ,"black")
ggplot(group_counts_melt, aes(x = Population, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Groups in Different Populations", x = "Superpopulation", y = "AF") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = colors) +  # 使用 Set3 调色板
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10, face = "bold"),  # 设置x轴标签的字体大小和粗细
    axis.text.y = element_text(size = 10, face = "bold")
  )


##########dra
hprcallrisky<-c("AcBcbCa","ABCbCa","AcBCbCa","AcBcba")
##
##AFR  AMR
type=1
rs=hprcallrisky[type]
hprcallrs<-distinct(hprcall[hprcall$oristrutype%in%rs &hprcall$population=="AFR" ,c(11,14,15)])
apgallrs<-distinct(apgall[apgall$oristrutype==rs ,c(11,14,15)])
hgsvcallrs<-distinct(hgsvcall[hgsvcall$oristrutype==rs ,c(11,14,15)])

a=nrow(apgallrs)
b=nrow(hprcallrs)
table_data <- matrix(c(a,b , length(unique(apgall$sample))-a, length(unique(hprcall[hprcall$population=="AFR" ,]$sample))-b), nrow = 2)

hprcallrs<-distinct(hprcall[hprcall$oristrutype%in%rs &hprcall$population=="AMR" ,c(11,14,15)])
c=nrow(hprcallrs)
a1=length(unique(apgall$sample))-a
b1=length(unique(hprcall[hprcall$population=="AFR" ,]$sample))-b
c1=length(unique(hprcall[hprcall$population=="AMR" ,]$sample))-c
x=matrix(c(a, a1,b,b1, c,c1),nrow = 2)

## EAS
type=4
rs=hprcallrisky[type]
hprcallrs<-distinct(hprcall[hprcall$oristrutype%in%rs &hprcall$population=="EAS" ,c(11,14,15)])
apgallrs<-distinct(apgall[apgall$oristrutype==rs ,c(11,14,15)])
hgsvcallrs<-distinct(hgsvcall[hgsvcall$oristrutype==rs &hgsvcall$population=="EAS" ,c(11,14,15)])
length(unique(hprcall[hprcall$population=="EAS" ,]$sample))
ength(unique(apgall$sample))
length(unique(hgsvcall[hgsvcall$population=="EAS" ,]$sample))

a=10
b=4
c=35
x=matrix(c(a,97-a,b,44-b,c,328-c),nrow = 2)
chisq.test(x)

### Odd ratio
data <- matrix(c(59, 38, 29, 15), nrow = 2, byrow = TRUE)
data <- matrix(c(29, 15, 234, 70), nrow = 2, byrow = TRUE)
data <- matrix(c(59, 38, 234, 70), nrow = 3, byrow = TRUE)
colnames(data) <- c("Case", "Control")
rownames(data) <- c("AFR", "AMR")
oddsratio(data, rev="c")

a=59
b=29
c=255
library(epitools)
data <- matrix(c(a, 97-a, b, 44-b), nrow = 2, byrow = TRUE)
data <- matrix(c(b, 44-b, c, 328-c), nrow = 2, byrow = TRUE)
data <- matrix(c(a, 97-a, c, 328-c), nrow = 2, byrow = TRUE)
colnames(data) <- c("Case", "Control")
rownames(data) <- c("AFR", "AMR")
oddsratio(data, rev="c")



data <- data.frame(
  Haplotype = rep(c("AcBcbCa", "ABCbCa", "AcBCbCa", "AcBcba", "ALL risky allele"), each = 3),
  Superpopulation = rep(c("AFR", "AMR", "EAS"), times = 5),
  Proportion = c(
    21.6, 4.5, 2.1, # AcBcbCa
    17.5, 45.45, 63.1, # ABCbCa
    10.3, 9.1, 10.6, # AcBCbCa
    11.3, 6.8, 1.8, # AcBcba
    60.8, 65.9, 77.7 # ALL risky allele
  )
)
ggplot(data, aes(x = Superpopulation, y = Proportion, fill = Superpopulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(. ~ Haplotype) +  # 使用 facet_grid 让图表横向排列
  labs(title = "Proportions of Haplotype in Superpopulations",
       x = "Superpopulation",
       y = "AF (%)") +
  scale_fill_brewer(palette = "Set5") +  # 使用颜色调色板
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1) # 为每个子图添加边框
  )





oristrutype_count <- table(hprcall$oristrutype)
oristrutype_prop <- prop.table(oristrutype_count)
oristrutype_df <- as.data.frame(oristrutype_prop)
colnames(oristrutype_df) <- c("oristrutype", "proportion")
#oristrutype_df<-oristrutype_df[oristrutype_df$oristrutype %in%  c("AcBcbCa","ABCbCa","AcBCbCa","AcBcba"),]
ggplot(oristrutype_df, aes(x = oristrutype, y = proportion, fill = oristrutype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.3f", proportion)), vjust = -0.5, size = 2.5) +
  labs(title = "Proportion of Different Strings in 'oristrutype'",
       x = "String",
       y = "Proportion") +
  theme_minimal()

oristrutype_population_count <- hprcall %>%
  group_by(oristrutype, population) %>%
  summarise(count = n()) %>%
  ungroup()

# 绘制堆叠条形图
ggplot(oristrutype_population_count, aes(x = oristrutype, y = count, fill = population)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Count of Different 'oristrutype' from Each 'population'",
       x = "oristrutype",
       y = "Count") +
  theme_minimal()














uniqsam=c()
for(samplename in unique(data$oristrutype)){
  sada<-data[data$oristrutype==samplename,]
  uniqsam<-c(uniqsam,sada$sample[1])
}
data<-data[data$sample %in%  uniqsam,]

# 安装和加载所需的包
install.packages("stringdist")
library(stringdist)

stringdata<-apg[,c(10,13)]
stringdata<-distinct(stringdata)
strings <- c("AcBcbCa", "ABCbCa", "ABCba", "ABCa", "AcBCba", "ABcbCa", "AcBCbCa", "AcBcba", "ABcba")
dist_matrix <- stringdistmatrix(stringdata$oristrutype, method = "lv")
hc <- hclust(as.dist(dist_matrix))
order <- hc$order
plot(hc, labels = stringdata$oristrutype, main = "Hierarchical Clustering of Strings", sub = "", xlab = "")
ordered_strings=stringdata$oristrutype[order]

data<-apg
data$sort_key <- match(data$oristrutype, ordered_strings)
data <- data[order(data$sort_key), ]



# samdf<-unlist(samdir)
# ## structure type
# 
# for (strgroup in unique(samdf)){
#   data[data$sample %in% names(samdf[samdf==strgroup]),]$strutype=strgroup
# }
# 
# for (strgroup in unique(data$strutype)){
#   sada<-data[data$strutype==strgroup,]
#   data[data$sample %in% names(samdf[samdf==strgroup]),]$strutype=num
#   num=num+1
# }
# similarity_matrix <- stringdistmatrix(samdf, method = "lv")
# distance_df <- as.data.frame(similarity_matrix)
# rownames(distance_df) <- colnames(distance_df) <- name
# 
# # 转换数据框以适应 ggplot2
# distance_melted <- melt(distance_df)

# 创建热图
ggplot(distance_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Edit Distance Heatmap", x = "Strings", y = "Strings", fill = "Edit Distance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
names(table(draw$molecule))[which.max(table(draw$molecule))]

library(ggplot2)
library(gggenes)
library(RColorBrewer)
draw<-data[,c(10,1,8,9,5)]
#draw<-fread("D:/MS/sd_trucdbscan.csv")
#draw<-fread("D:/test.txt")
#draw$strand="+"
colnames(draw)<-c("molecule","gene","start","end","strand")
draw<-draw[order(draw$molecul),]


#colnames(draw)<-c("molecule","gene","start","end","strand")

drawnow$molecule <- factor(draw$molecule, levels = unique(draw$molecule))
draw$direction <- ifelse(draw$strand == "+", 0,1)
drawnow<-draw
#drawnow<-draw
drawnow<-drawnow[,c("molecule","gene","start","end","direction")]
pdf("D:/MS/SD/apg/sd/hprc.pdf",width = 10,height = 80)
drawnow$molecule <- sub("\\..*", "", drawnow$molecule)
ggplot(drawnow, aes(xmin = start, xmax = end, y = factor(molecule), fill = gene, forward = direction)) +
  geom_gene_arrow(alpha = 1, arrow_body_height = unit(2.5, "mm"), arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(limits = c(min(drawnow$start), max(drawnow$end))) + # 设置整体x坐标范围
  theme_genes()+
  theme(
    axis.line.x = element_blank(),       # 去掉 x 轴的坐标轴线
    axis.text.x = element_blank(),       # 去掉 x 轴的刻度标签
    axis.ticks.x = element_blank(),      # 去掉 x 轴的刻度线
    axis.title.x = element_blank()       # 去掉 x 轴的标题
  )

dev.off()
library(data.table)
data<-fread("D:/MS/SD/test.csv")
data<-data[data$`chr2:110863736-111057513`!="",]
data<-data[data$end_position-data$start_position>50,]
plot(data$`chr2:110863736-111057513`)
