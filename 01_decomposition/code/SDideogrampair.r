library(data.table)
library(karyoploteR)
genomesize<-fread("D:/chm13 -cp.size")
genomesize2<-fread("D:/cytobondsall")
colnames(genomesize2)<-c("chr","start","end","name","gieStain")
genomesize2=toGRanges(genomesize2)
genomesize=toGRanges(genomesize)
data<-fread("D:/MS/SD/apg/kp/t2t/weno.bed")
data1<-fread("D:/MS/SD/apg/kp/t2t/wemore.bed")
# colnames(data1)<-c("chr","start","end")
# png(width =10000,height =7000,paste("/share/home/zhanglab/user/maoyafei/project/flaggervisulize/",fold,"_",gsub(".bed", "", m),".png",sep = ""),res = 300)
kp <- plotKaryotype(genome=genomesize,cytobands = genomesize2,chromosomes=paste0("chr", c(1:22)))
kpPlotRegions(kp, data,col="red", r0=0, r1=0.5)
kpPlotRegions(kp, data1, col="black",r0=0.5, r1=0.9)
kpPlotRegions(kp, data3, col="#AACCFF", ,r0=0.3, r1=0.5,layer.margin = 0.4)


library(karyoploteR)

sedef<-fread("D:/MS/SD/apg/kp/38/kpsedef.bed")
our<-fread("D:/MS/SD/apg/kp/38/kpgrch38.bed")
wemore<-fread("D:/MS/SD/apg/kp/38/wemoremer.bed")
weno<-fread("D:/MS/SD/apg/kp/38/wenomer.bed")
kp <- plotKaryotype(genome = "hg38")
kpPlotRegions(kp, sedef,col="black", r0=0, r1=0.5)
kpPlotRegions(kp, our, col="red",r0=0.5, r1=0.9)
kpPlotRegions(kp, weno,col="#AACCFF", r0=0, r1=0.5)
kpPlotRegions(kp, wemore, col="#a9d18e",r0=0.5, r1=0.9)


wemore<-fread("D:/MS/SD/apg/kp/38/wemorepair.bed")
wemoreintra<-wemore[wemore$V1==wemore$V4,]
weno<-fread("D:/MS/SD/apg/kp/t2t/whynopair.bed")
wenointra<-weno[weno$V1==weno$V4,]
kp <- plotKaryotype(genome = "hg38")
kpPlotRegions(kp, sedef,col="black", r0=0, r1=0.5)
kpPlotRegions(kp, our, col="red",r0=0.5, r1=0.9)
kpPlotRegions(kp, weno,col="#AACCFF", r0=0, r1=0.5)
kpPlotRegions(kp, wemore, col="#a9d18e",r0=0.5, r1=0.9)
weno<-fread("D:/MS/SD/apg/kp/438/b.bed")
#weno<-weno[weno$V1==weno$V4,]
weno<-weno[weno$V1=="chr1",]
start.regs <- toGRanges(weno[,c(1,2,3)])
end.regs <- toGRanges(weno[,c(4,5,6)])

kp <- plotKaryotype(genome=genomesize,cytobands = genomesize2,chromosomes=paste0("chr", c(1:22)))

kp <- plotKaryotype(chromosomes = c("chr4","chr10"))
TP53.region <- toGRanges("chr4:193,390,000-193,590,")
kp <- plotKaryotype(genome=genomesize,cytobands = genomesize2,chromosomes=c("chr4","chr10"))

kpPlotLinks(kp, data = start.regs, data2 = end.regs, col = "#FF000080") # 使用半透明的红色
kpAddBaseNumbers(kp, tick.dist = 100000, add.units = TRUE, cex = 0.6)

library(circlize)
wenointer<-weno[weno$V1!=weno$V4,]
# 加载circlize包
library(circlize)
wenointer<-wenointer[1:10,]
# 初始化数据
wenointer <- weno[weno$V1 != weno$V4, ]
bed1 <- wenointer[, c(1, 2, 3)]
bed2 <- wenointer[, c(4, 5, 6)]

# 初始化带有染色体图的圈图
circos.initializeWithIdeogram()

# 为每条染色体分配颜色
chromosomes <- unique(c(bed1$V1, bed2$V1))
chromosome_colors <- rand_color(length(chromosomes))
names(chromosome_colors) <- chromosomes

# 为每条连接线分配对应染色体的颜色
col_vec <- chromosome_colors[bed1$V1]

# 绘制连接线并设置线条宽度
circos.genomicLink(
  bed1, bed2,
  col = col_vec,
  lwd = 0.5
)

# 完成图形
circos.clear()


# 清除图像，如果需要的话
# circos.clear()

ggplot(data) +
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2)) +
  labs(x = "chr1", y = "chr3") +
  xlim(min(data$x1)-1000, max(data$x2)+1000) +  # 设置 X 轴范围
  ylim(min(data$y1)-1000, max(data$y2)+1000) +  # 设置 Y 轴范围
  theme_minimal()

