library(data.table)
library(dplyr)
rm(alldf)
allsam<-fread(".DATA/APG/allsam",header=FALSE)
for(sam in allsam$V1){
	print(sam)
	infile=paste(".PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/",sam,"/",sam,".t2tall",sep="")
	samdf<-fread(infile)
	samdf<-distinct(samdf)
	if (!exists("alldf")) {
    alldf <- samdf
  } else {
    alldf <- merge(alldf, samdf, by = c("A", "B"), all.x = TRUE)
  }
}

fwrite(alldf,"allsam.csv",sep=",")


### group
allsam<-fread(".DATA/APG/allsam",header=FALSE)
alldf<-fread("allsam.csv",sep=",")
samsuper=fread(".DATA/APG/sam_pop.csv")
alldfstatisc=alldf

for(sup in unique(na.omit(samsuper$SP))){
	group=samsuper[samsuper$SP==sup,]$V1
	grouplen=length(alldf[, ..group])
	no1 <- apply(alldf[, ..group], 1, function(x) sum(x == "NOSAM"))
	no2 <- apply(alldf[, ..group], 1, function(x) sum(x == "NO"))
	noall=no1+no2
	have=grouplen-noall
	haveall=have/grouplen
	var=paste(sup,"_fre",sep="")
	alldfstatisc[[var]]=haveall
	var=paste(sup,"_have",sep="")
	alldfstatisc[[var]]=have
	var=paste(sup,"_no",sep="")
	alldfstatisc[[var]]=noall
}

###fisher
###  "EAS" "EUR" "AMR" "AFR" "SAS"
### EAS&EUR
### EAS&AMR
### EAS&AFR
### EAS&SAS
### EUR&AMR
### EUR&AFR
### EUR&SAS
### AMR&AFR
### AMR&SAS
### AFR&SAS


allx=c()
for(i in 1:nrow(alldfstatisc)){
print(i)
eas1=alldfstatisc[i,]$EAS_have
eas0=alldfstatisc[i,]$EAS_no
eur1=alldfstatisc[i,]$EUR_have
eur0=alldfstatisc[i,]$EUR_no
amr1=alldfstatisc[i,]$AMR_have
amr0=alldfstatisc[i,]$AMR_no
afr1=alldfstatisc[i,]$AFR_have
afr0=alldfstatisc[i,]$AFR_no
sas1=alldfstatisc[i,]$SAS_have
sas0=alldfstatisc[i,]$SAS_no
table_eas_eur <- matrix(c(eas1,eas0, eur1, eur0), nrow = 2, byrow = TRUE)
table_eas_amr <- matrix(c(eas1,eas0, amr1, amr0), nrow = 2, byrow = TRUE)
table_eas_afr <- matrix(c(eas1,eas0, afr1,afr0), nrow = 2, byrow = TRUE)
table_eas_sas <- matrix(c(eas1,eas0, sas1, sas0), nrow = 2, byrow = TRUE)
table_eur_amr <- matrix(c(eur1, eur0, amr1, amr0), nrow = 2, byrow = TRUE)
table_eur_afr <- matrix(c(eur1, eur0, afr1,afr0), nrow = 2, byrow = TRUE)
table_eur_sas <- matrix(c(eur1, eur0, sas1, sas0), nrow = 2, byrow = TRUE)
table_amr_afr <- matrix(c(amr1, amr0, afr1,afr0), nrow = 2, byrow = TRUE)
table_amr_sas <- matrix(c(amr1, amr0, sas1, sas0), nrow = 2, byrow = TRUE)
table_afr_sas <- matrix(c(afr1,afr0, sas1, sas0), nrow = 2, byrow = TRUE)
tables <- list(table_eas_eur, table_eas_amr, table_eas_afr, table_eas_sas,table_eur_amr,table_eur_afr,table_eur_sas,table_amr_afr,table_amr_sas,table_afr_sas)
p_values <- sapply(tables, function(table) {
  fisher.test(table)$p.value
})
adjusted_p_values_bonferroni <- p.adjust(p_values, method = "bonferroni")
allx=rbind(allx,adjusted_p_values_bonferroni[1:4])
}

allxdf=data.frame(allx)
colnames(allxdf)<-c("eas_eur_p","eas_amr_p","eas_afr_p","eas_sas_p")




freall<-alldfstatisc[,c("A","B","EAS_fre","EUR_fre","AMR_fre","AFR_fre","SAS_fre")]
cols_to_normalize <- freall[, 3:7]
x<- as.data.frame(t(apply(cols_to_normalize, 1, function(x) x / sum(x))))
colnames(x)<-c("EAS_norm","EUR_norm","AMR_norm","AFR_norm","SAS_norm")
freall<-alldfstatisc
alldf<-cbind(freall,x)
alldf<-cbind(alldf,allxdf)
##add identity
t2tsdiden=".PANSDEND/01_REFSD/para/SD.refine"
t2tsdidendf<-fread(t2tsdiden)
t2tadd<-t2tsdidendf[,c(1,2,3,4,5,6,10,11)]
colnames(t2tadd)<-c("Ac","As","Ae","Bc","Bs","Be","iden","len")
t2tadd$A <- paste(t2tadd$Ac, ":", t2tadd$As, "-", t2tadd$Ae, sep = "")
t2tadd$B <- paste(t2tadd$Bc, ":", t2tadd$Bs, "-", t2tadd$Be, sep = "")
t2tadd<-t2tadd[,c("A","B","iden","len")]
t2taddcp<-t2tadd
colnames(t2taddcp)=c("B","A","iden","len")
t2tadddb<-distinct(rbind(t2tadd,t2taddcp))
alldfadlab<-merge(alldf,t2tadddb,by=c("A","B"),all.x=TRUE)
fwrite(alldfadlab,"eas_fisher_all.csv",sep=",")


library(data.table)
library(dplyr)
alldfadlab<-fread("eas_fisher_all.csv")
alldfadlab$count_lt_0.05 <- rowSums(alldfadlab[, c("eas_eur_p", "eas_amr_p", "eas_afr_p", "eas_sas_p")] < 0.05)
fwrite(alldfadlab,"eas_fisher_all.csv",sep=",")



alldf1=alldfadlab[alldfadlab$eas_eur_p<0.05 |alldfadlab$eas_amr_p<0.05|alldfadlab$eas_afr_p<0.05|alldfadlab$eas_sas_p<0.05, ]
alldf2=alldfadlab[alldfadlab$eas_eur_p<0.05  & alldfadlab$eas_amr_p<0.05 & alldfadlab$eas_afr_p<0.05 &alldfadlab$eas_sas_p<0.05, ]

fwrite(alldf1,"eas_fisher.csv",sep=",")
fwrite(alldf2,"eas_fisheralsig.csv",sep=",")



# colnames(alldf)<-c("A","B","EAS_fre","EUR_fre","AMR_fre","AFR_fre","SAS_fre",)

# alldfcopy=alldf
# alldfcopy$max_values=0
# alldfcopy$min_values=0
# alldfcopy$max_values <- apply(alldf[, (ncol(alldf)-5):ncol(alldf)], 1, max)
# alldfcopy$min_values <- apply(alldf[, (ncol(alldf)-5):ncol(alldf)], 1, min)
# eashigh=alldfcopy[alldfcopy$EAS_norm==alldfcopy$max_values,]
# easlow=alldfcopy[alldfcopy$EAS_norm==alldfcopy$min_values,]
# eashigh1=eashigh[eashigh$max_values>0.2,]
# easlow1=easlow[easlow$min_values<0.2,]
# easlow1=easlow1[easlow1$max_values-easlow1$min_values>0.05,]


fwrite(eashigh1,"./eashigh1.csv",sep=",")
fwrite(easlow1,"./pangenome/Pan-SD/pansam/easlow1.csv",sep=",")


###find
eashigh1=fread("./eashigh1.csv")
alldf=fread("allsam.csv")
eashigh1alldf <- merge(eashigh1, alldf, by = c("A", "B"), all.x = TRUE)
fwrite(eashigh1alldf,"eashigh1alldf.csv",sep="\t")
eashigh1alldf=fread("eashigh1alldf.csv",sep="\t")
t2t=fread("./pangenome/APG_seg/T2TSD/pansd.t2t")
t2t=t2t[,c(1,2,3,4,5,6,10)]
colnames(t2t)[7]="iden"
t2t$Amer <- paste(t2t$V1, t2t$V2, sep = ":")
t2t$A <- paste(t2t$Amer, t2t$V3, sep = "-")
t2t$Bmer <- paste(t2t$V4, t2t$V5, sep = ":")
t2t$B <- paste(t2t$Bmer, t2t$V6, sep = "-")
t2t<-t2t[,c("iden","A","B")]
eashigh1alldfien <- merge(eashigh1alldf, t2t, by = c("A", "B"), all.x = TRUE)

fwrite(eashigh1alldfien[eashigh1alldfien$iden>0.92,],"eashigh1alldf1.csv",sep="\t")



easlow1=fread("./pangenome/Pan-SD/pansam/easlow1.csv")
alldf=fread("allsam.csv")
easlow1alldf <- merge(easlow1, alldf, by = c("A", "B"), all.x = TRUE)
t2t=fread("./pangenome/APG_seg/T2TSD/pansd.t2t")
t2t=t2t[,c(1,2,3,4,5,6,10)]
colnames(t2t)[7]="iden"
t2t$Amer <- paste(t2t$V1, t2t$V2, sep = ":")
t2t$A <- paste(t2t$Amer, t2t$V3, sep = "-")
t2t$Bmer <- paste(t2t$V4, t2t$V5, sep = ":")
t2t$B <- paste(t2t$Bmer, t2t$V6, sep = "-")
t2t<-t2t[,c("iden","A","B")]
easlow1aien <- merge(easlow1alldf, t2t, by = c("A", "B"), all.x = TRUE)

fwrite(easlow1aien ,"easlowh1alldf1.csv",sep="\t")


# for(i in 1:nrow(alldf)){
# 	cho=alldf[i,]
# }
# alldf <- subset(alldf, select = -CHM13v2)
# test<-alldf
# alldfstatisc=alldf
# alldfstatisc$HAVE_count <- apply(alldf, 1, function(x) sum(x == "HAVE"))
# alldfstatisc$NO_count <- apply(alldf, 1, function(x) sum(x == "NO"))
# alldfstatisc$NOSAM_count <- apply(alldf, 1, function(x) sum(x == "NOSAM"))
# alldfstatisc$nono<-alldfstatisc$NO_count +alldfstatisc$NOSAM_count