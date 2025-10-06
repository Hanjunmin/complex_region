library(data.table)
library(dplyr)
data<-fread("C:/Users/hjm00000/Desktop/sbnet.csv")
data<-data[,-c(7,8,23,24)]
sd1<-data[,1]
sd2<-data[,15]
data[data > 0.95] <- 1
data[data <= 0.95] <- 0
data[is.na(data)] <- 0
data[,1]=sd1
data[,15]=sd2
data<-distinct(data)
data$sd1len<-rowSums(data[, 3:14])
data$sd1group1<-rowSums(data[, 3:6])
data$sd1group2<-rowSums(data[, 11:14])
data$sd1group3<-rowSums(data[, 7:10])
data$sd2len<-rowSums(data[, 17:28])
data$sd2group1<-rowSums(data[, 17:20])
data$sd2group2<-rowSums(data[, 25:28])
data$sd2group3<-rowSums(data[, 21:24])
node1<-data[,c(1,29,30,31,32)]
node2<-data[,c(15,33,34,35,36)]
node1$sd1group1<-node1$sd1group1/node1$sd1len
node1$sd1group2<-node1$sd1group2/node1$sd1len
node1$sd1group3<-node1$sd1group3/node1$sd1len
node1[is.na(node1)] <- 0


node2$sd2group1<-node2$sd2group1/node2$sd2len
node2$sd2group2<-node2$sd2group2/node2$sd2len
node2$sd2group3<-node2$sd2group3/node2$sd2len
node2[is.na(node2)] <- 0


node1<-distinct(node1)
node2<-distinct(node2)
data<-as.data.frame(data)
alld<-data.frame()
for (i in 1:nrow(data)) {
  for (j in 3:14) {
    if(data[i,j]==1 & data[i,j]==data[i,j+14]){
      alld<-rbind(alld,c(data[i,1],data[i,15],colnames(data)[j]))
    }
  }
}
colnames(alld)<-c("n1","n2","sam")
alld$pop="no"
alld[alld$sam %in% c("HG00438_1","HG00438_2","HG00621_1","HG00621_2"),]$pop="group1"
alld[alld$sam %in% c("HG02622_1","HG02622_2","HG03516_1","HG03516_2"),]$pop="group2"
alld[alld$sam %in% c("HG00735_1","HG00735_2","HG01123_1","HG01123_2"),]$pop="group3"

result <- alld %>%
  group_by(n1, n2, pop) %>%
  summarise(count = n(), .groups = 'drop')
colnames(node2)<-colnames(node1)
node<-rbind(node1,node2)
node<-distinct(node)

aledeg<-data[,c(1,15)]
nedg<-result[,c(1,2)]
colnames(aledeg)<-colnames(nedg)
result <- anti_join(aledeg, nedg)

fwrite(node,"D:/MS/node.txt")
fwrite(result,"D:/MS/edge.txt")
fwrite(result,"D:/MS/othedge.txt",col.names = TRUE)
