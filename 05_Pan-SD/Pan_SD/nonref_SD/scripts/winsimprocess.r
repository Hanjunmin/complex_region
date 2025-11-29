library(data.table)
args <- commandArgs(trailingOnly = TRUE)
input=args[1]
output=args[2]
initW<-"init.W"
initWdf <- fread(initW, header=FALSE)
colnames(initWdf)=c("Sample","lab","Contig","contigs","contige","pc")
T2Tsamfile <- fread(input) ##"{win_simdir}{wildcards.li}sd.winsim.bed"
T2Tsamfileadd<-merge(T2Tsamfile,initWdf,by = "pc", all.x = TRUE)
fwrite(T2Tsamfileadd,output,sep="\t",col.names=FALSE) ##"{win_simdir}{wildcards.li}sd.winsim.bed.add"