rev_mer <- function(df) {
  dfrev=df[,c(2,1,3)]
  colnames(dfrev)=colnames(df)
  dfall<-rbind(df,dfrev)
  return(dfall)
}

args <- commandArgs(trailingOnly = TRUE)
library(data.table)
library(dplyr)
sam=args[1]
edgeall=".PANSDEND/01_REFSD/para/edge.bed"
data<-fread(edgeall,header=FALSE)
colnames(data)<-c("A","B")

### 2 link file
mini="./t2tsdminimap/all.end.statistics.filt"
minidf<-fread(mini,header=FALSE)
sam2pos="T2T2SAM"
sam2posdf<-fread(sam2pos,header=FALSE)
colnames(sam2posdf)=c("c1","cs1","ce1","c2","cs2","ce2","t2tfrom","t2tto")
sam2posdf[, samlef := paste(c1, ":", cs1, "-", ce1, sep="")]
sam2posdf[, samrig := paste(c2, ":", cs2, "-", ce2, sep="")]
sam2posdf<-sam2posdf[,c("t2tfrom","t2tto","samlef","samrig")]

###sim.txt
sim<-fread("sim.txt",header=FALSE)
sim<-sim[,c(1,2,4,6)]
colnames(sim)<-c("t2tfrom","t2tto","t2tto.sim","t2tfrom.sim")
sim=distinct(sim)
sam2posdf=distinct(sam2posdf)
sam2posdf_sim<-merge(sam2posdf,sim,by=c("t2tfrom","t2tto"),all.x=TRUE)
### 1.processnosam.txt (nosamdf)
nosam="processnosam.txt"
nosamdf<-fread(nosam,header=FALSE)
nosamdf$sam1="NOSAM"
nosamdf<-nosamdf[,c(1,2,7)] ###之后可以从这里看究竟是哪个区域没有
nosamdf<-distinct(nosamdf)
colnames(nosamdf)<-c("A","B","sam1")
mer1=rev_mer(nosamdf)



### 2.T2Thave.bed (t2thavedfinitpair)
t2thave="./t2tsdminimap/T2Thave.bed"
t2thavedf<-fread(t2thave,header=FALSE)
t2thavedf<-t2thavedf[,c(13,14)]
colnames(t2thavedf)<-c("samrig","samlef")
t2thavedf<-distinct(t2thavedf)
t2thavedfinitpair=merge(sam2posdf_sim,t2thavedf,by=c("samlef","samrig"))
truehave=t2thavedfinitpair[t2thavedfinitpair$t2tto.sim>=0.8 & t2thavedfinitpair$t2tfrom.sim>=0.8,c("t2tfrom","t2tto")]
otherhave=t2thavedfinitpair[!(t2thavedfinitpair$t2tto.sim>=0.8 & t2thavedfinitpair$t2tfrom.sim>=0.8),c("t2tfrom","t2tto")]
truehave$sam2="HAVE"
otherhave$sam5="POLY"
truehave<-distinct(truehave)
otherhave<-distinct(otherhave)
colnames(truehave)<-c("A","B","sam2")
colnames(otherhave)<-c("A","B","sam5")
mer2=rev_mer(truehave)
mer2_1=rev_mer(otherhave)

### 3.T2Tpolymor1.bed
t2tpoly="./t2tsdminimap/T2Tpolymor1.bed"
t2tpolydf<-fread(t2tpoly,header=FALSE)
t2tpolydf<-t2tpolydf[,c(13,14)]
colnames(t2tpolydf)<-c("samrig","samlef")
t2tpolydf<-distinct(t2tpolydf)
t2tpolydfinitpair=merge(sam2posdf,t2tpolydf,by=c("samlef","samrig"))
t2tpolydfinitpair<-t2tpolydfinitpair[,c("t2tfrom","t2tto")]
t2tpolydfinitpair$sam3="POLY"
t2tpolydfinitpair<-distinct(t2tpolydfinitpair)
colnames(t2tpolydfinitpair)<-c("A","B","sam3")
mer3=rev_mer(t2tpolydfinitpair)


### 4.T2Tpolymor1.bed
rmrefiltpoly="./t2tsdminimap/filin.bedend"
rmrefiltpolydf<-fread(rmrefiltpoly,header=FALSE)
rmrefiltpolydfal=merge(rmrefiltpolydf[,c(1,2,3,4,5,6)],minidf,by=c("V1","V2","V3","V4","V5","V6"),all.x=TRUE)
rmrefiltpolydfal<-rmrefiltpolydfal[,c(13,14)]
colnames(rmrefiltpolydfal)<-c("samrig","samlef")
rmrefiltpolydfal<-distinct(rmrefiltpolydfal)
rmrefiltpolydfinitpair=merge(sam2posdf,rmrefiltpolydfal,by=c("samlef","samrig"))
rmrefiltpolydfinitpair<-rmrefiltpolydfinitpair[,c("t2tfrom","t2tto")]
rmrefiltpolydfinitpair$sam4="POLY2"
rmrefiltpolydfinitpair<-distinct(rmrefiltpolydfinitpair)
colnames(rmrefiltpolydfinitpair)<-c("A","B","sam4")
mer4=rev_mer(rmrefiltpolydfinitpair)

##

allmer=merge(data,mer1,by=c("A","B"),all.x=TRUE)
allmer=merge(allmer,mer2,by=c("A","B"),all.x=TRUE)
allmer=merge(allmer,mer2_1,by=c("A","B"),all.x=TRUE)
allmer=merge(allmer,mer3,by=c("A","B"),all.x=TRUE)
allmer=merge(allmer,mer4,by=c("A","B"),all.x=TRUE)
allmer[] <- lapply(allmer, function(x) ifelse(is.na(x), "", x))
allmer[, al := paste(sam1,sam2,sam3,sam4,sam5,sep="")]
###unique(allmer$al) ???budui
allmer[allmer$al=="",]$al="NO"
allmer<-allmer[,c("A","B","al")]
allmer$al[grepl("HAVE", allmer$al)] <- "HAVE"
allmer$al[grepl("POLY", allmer$al)] <- "POLY"

colnames(allmer)[3]=sam

outfile=paste(sam,".t2tall",sep="")
fwrite(allmer,outfile,sep="\t")
