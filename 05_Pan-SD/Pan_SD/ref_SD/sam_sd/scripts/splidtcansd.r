library(data.table)

args <- commandArgs(trailingOnly = TRUE)
T2Tsd <-args[1] #"/home/jmhan/PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/alldfe.txt" #
T2Tcan<-args[2] #"/home/jmhan/PANSDEND/02_SAMSD/02_00-REFSD/1_SD_CANPATH/alldfe1.txt" #
T2Tcandf <- fread(T2Tcan, header=FALSE)
T2Tsddf <- fread(T2Tsd, header=FALSE)
canid=T2Tcandf$V1
sdid=T2Tsddf$V1
sam<- args[3] #"C022-CHA-S02_2-Mat" #
samfile<-args[4] #"x" #
initW<-args[5] #"init.W" #
initWdf <- fread(initW, header=FALSE)
colnames(initWdf)=c("Sample","lab","Contig","contigs","contige","pc")
T2Tsamfile <- fread(samfile,header=FALSE)
colnames(T2Tsamfile)<-c("SDregion","ps","pe","pc","row","sim")
T2Tsamfile$com <- ifelse(grepl("complex", T2Tsamfile$SDregion), "COMPLEX", "SIMPLE")
T2Tsamfile$SDregion <- gsub("_complex", "", T2Tsamfile$SDregion)
T2Tsamfileadd<-merge(T2Tsamfile,initWdf,by = "pc", all.x = TRUE)
T2Tsamfileadd$row <- NULL
splitcan=T2Tsamfileadd[T2Tsamfileadd$SDregion %in% canid,]
splitsd=T2Tsamfileadd[T2Tsamfileadd$SDregion %in% sdid,]
namesd=paste(sam,".sd",sep="")
namecan=paste(sam,".can",sep="")
fwrite(splitcan,namecan,sep="\t",col.names=TRUE)
fwrite(splitsd,namesd,sep="\t",col.names=TRUE)