library(data.table)
library(IRanges)

data<-fread("chlist",header=FALSE)
chrlist=data$V1
for(chrsring in chrlist){
	chrafterpath=paste(chrsring,"_after.csv",sep="")
	chrotherpath=paste(chrsring,"_other.csv",sep="")
	if (file.exists(chrafterpath) & file.exists(chrotherpath)) {
  	chrafter=fread(chrafterpath)
	chrother=fread(chrotherpath)
	nomatch=setdiff(chrother$cluster,chrafter$cluster)
	match=intersect(chrother$cluster,chrafter$cluster)
	if(length(nomatch)!=0){
		nomatchdata<-chrother[chrother$cluster %in% nomatch,]
		fwrite(nomatchdata, file =paste(chrsring,"_nomatch.csv",sep=""))
	}
	allclus=data.frame()
	for(clus in match){
		matchafter<- chrafter[ chrafter$cluster==clus,]
		matchafter<-matchafter[matchafter$V10>=0.9 & matchafter$V12 >=0.5 & matchafter$V11>=1000,]
		matchother<- chrother[ chrother$cluster==clus,]
		irafter <- IRanges(start = matchafter$V2, end = matchafter$V3)
		irother <- IRanges(start = matchother$V2, end = matchother$V3)
		queirafter <- IRanges(start = matchafter$V5, end = matchafter$V6)
		queirother <- IRanges(start = matchother$V5, end = matchother$V6)
		merged_irafter <- reduce(irafter)
		merged_irother <- reduce(irother)
		quemerged_irafter <- reduce(queirafter)
		quemerged_irother <- reduce(queirother)
		if (!((sum(merged_irafter@width)/sum(merged_irother@width))>=0.95 & (sum(quemerged_irafter@width)/sum(quemerged_irother@width))>=0.95)){allclus<-rbind(allclus,matchother)}
		# if(min(matchafter$V10)<0.9 | min(matchafter$V12)<0.5 | min(matchafter$V11)<1000){
		# 	allclus=c(allclus,clus)
		# }
	}
	#matchdata<-chrother[chrother$cluster %in% allclus,]
	fwrite(allclus, file =paste(chrsring,"_beforeminus0.9.csv",sep=""))
	}	
}