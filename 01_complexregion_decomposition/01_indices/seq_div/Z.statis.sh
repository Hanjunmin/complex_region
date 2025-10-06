

###1. make windows
cd ./COMPLEXFEATURE/
bedtools makewindows -g <(less -S ./breakpoints/chm13v2.0.fa.fai  |awk 'OFS="\t"{print $1,$2}') -w 1000 -s 500  >WGS.window.bed
#bedtools subtract -a WGS.window.bed -b ./COMPLEXFEATURE/data/CHM13v2.full.cent.bed -A >WGS.window.delcen.bed
bedtools subtract -a WGS.window.bed -b ./DATA/REF/chm13v2.0.refine.sat.bed -A >WGS.window.delcen.bed



###2.vcf (./allvcf/)
fold="./COMPLEXFEATURE/00_vcf/"
vcffold="./allvcf/"
for m in {1..22} X Y;do
	mkdir $fold"chr"$m
	cd $fold"chr"$m
	rm *
	cat ./COMPLEXFEATURE/WGS.window.delcen.bed |awk -v chr=chr$m '$1==chr{print $0}'  >chr${m}.size
	bash  ${fold}paravcf.sh $vcffold"chr"${m}.vcf $m &  ##
done 

cd ./COMPLEXFEATURE/00_vcf/


cat */end.txt |less -S |bedtools sort >allend
# cat */entropyend.txt |less -S |bedtools sort >entropyendall
cat */entropyendn.txt |less -S |bedtools sort >entropyendall
cat */entropyendn01.txt |less -S |bedtools sort >entropyendall01 ###entropy/(maxmixentropy)
types=("SNV" "indel" "INS" "DEL" "smallNM" "bigNM" "smallcomplex" "bigcomplex")
for type in "${types[@]}"; do
	cat */"${type}.window.bed" |bedtools sort > "${type}.featurein.txt" &
done

cat */add.txt >small_complex_NM_indel.txt  ###小于50bp区域的情况
cat */simal |grep -v "sim" >sim_al.bed     ### 计算每个窗口的相似性

Rscript -e 'library(data.table)
sigmoid <- function(x, k, mu) {
  return(1 / (1 + exp(-k * (x - mu))))
}
data<-fread("small_complex_NM_indel.txt")
data$new1=sigmoid(data$V4,3,1)
data$lennorm=sigmoid(data$V5,1,10)
data$end=data$lennorm*data$new1+0.05
plot(density(data$V5))
plot(density(data$new1))
plot(density(data[data$V4<1,]$V4))

fwrite(data,"small_complex_NM_indel.norm",sep="\t",col.names=FALSE)'

Rscript ./COMPLEXFEATURE/intevscore.r

###3.SD
cd ./COMPLEXFEATURE/04_SD/
rsync -av clsmyf-user1@sylogin.hpc.sjtu.edu.cn:/dssg/home/acct-clsmyf/clsmyf-user1/pangenome/APG_seg/output.txt ./APGsdsim.bed
cat APGsdsim.bed |tr ":" "\t" |tr "-" "\t" |less -S |tail -n +2 |less -S |cut -f 1-4 >APGsdsim.bedend


###4.VNTR
cd ./COMPLEXFEATURE/03_VNTR/
# cp ./pangenome/WGSfea/VNTR/TRF_VAMOStest/ben/ENDBEN/trf ./  ##之后可以把other的加上
# cat trf |sed 's/\.0//g' |cut -f 1-3 >PANVNTR.REF  ## ./WGSfeatureall/vntrindex.py
# bedtools coverage -a ../WGS.window.delcen.bed -b   PANVNTR.REF  |cut -f 1,2,3,7 >VNTR.cov
cp ./pangenome/WGSfea/VNTR/TRF_VAMOStest/ben/ENDBEN/trf ./  ##之后可以把other的加上
# cat trf |sed 's/\.0//g' |cut -f 1-3 >PANVNTR.REF  ## ./WGSfeatureall/vntrindex.py
# bedtools coverage -a ../WGS.window.delcen.bed -b   PANVNTR.REF  |cut -f 1,2,3,7 >VNTR.cov
cat  ./pangenome/WGSfea/VNTR/TRF_VAMOStest/panvntr/end/STRVNTR120/*.STR_VNTR.out  |awk '$4>2 && $4*$6>=150 && $4*$6<10000 && $6>6'|cut -f 1-3 >PANVNTR.REF
edtools coverage -a ../WGS.window.delcen.bed -b   PANVNTR.REF  |cut -f 1,2,3,7 >VNTR.cov




###5.bubble
cd ./COMPLEXFEATURE/01_bubble
#./pangenome/WGSfea/bubble/111/endbubblerun.sh
#cp ./pangenome/WGSfea/bubble/1210/all.bubblenew ./
#cp ./pangenome/WGSfea/bubble/424/all.bubblenew ./
./pangenome/WGSfea/bubble/424/endbubblerun.sh

###对连续的bubble进行统计
less -S all.bubblenew  |bedtools sort >all.bubblenew.sort


Rscript -e 'library(data.table)
data<-fread("all.bubblenew.sort")
data[, group := rleid(V4)]
group_counts <- data[, .N, by = group]
serial <- group_counts[N >=3]
out=unique(data[data$V4>300 & data$group %in% serial$group,]$group)
end=data[data$group %in%out,]
fwrite(end,"refine_5window_500.region",sep="\t")'  ###连续区域>=5k且simplebubble大于500的重新去给权重确定对应的情况
less -S refine_5window_500.region |cut -f 1-3 |tail -n +2 |bedtools sort |bedtools merge >refine.bubble.region

###refinebubble
cd ./COMPLEXFEATURE/01_bubble/para_gene_path/
for chr in {1..22} X Y;do
	mkdir -p ./COMPLEXFEATURE/01_bubble/para_gene_path/"chr"$chr
	cd ./COMPLEXFEATURE/01_bubble/para_gene_path/"chr"$chr
	cat ./COMPLEXFEATURE/01_bubble/refine.bubble.region |awk -v chrn="chr"$chr '{OFS="\t"}$1==chrn' >region.bed
done


###运行./COMPLEXFEATURE/01_bubble/sub.sh ##./COMPLEXFEATURE/01_bubble/resub.sh
###最终生成了end.txt


###6.总结
cd ./COMPLEXFEATURE/
bedtools coverage -a ./COMPLEXFEATURE/WGS.window.delcen.bed -b ~/allvcf/alvcfbed  |awk '$7!=0' |cut -f 1-3  >WGS.window.delcen.delnvariants.bed

cd ./COMPLEXFEATURE/05_split/
###R

library(data.table)
data<-fread("../WGS.window.delcen.delnvariants.bed")
colnames(data)<-c("win_chr","win_s","win_e")
# sddf=fread("../04_SD/APGsdsim.bedend")
sddf=fread("../04_SD/SDsim.n")
colnames(sddf)<-c("win_chr","win_s","win_e","sd_max",'sd_var')
sddf$win_s=sddf$win_s-1
data=merge(data, sddf, by = c("win_chr","win_s","win_e"),all.x=TRUE)
vntrdf=fread("../03_VNTR/VNTR.cov")
colnames(vntrdf)<-c("win_chr","win_s","win_e","vntr")
data=merge(data, vntrdf, by = c("win_chr","win_s","win_e"),all.x=TRUE)

# data2=fread("../00_vcf/vcf_feature.csv")
data2=fread("../00_vcf/entropyendall01")
colnames(data2)<-c("win_chr","win_s","win_e","entropy")
data=merge(data, data2, by = c("win_chr","win_s","win_e"),all.x=TRUE)

dataTE=fread("./COMPLEXFEATURE/09_bubble_only/03_model/TE.ALLchr")
dataTE[is.na(dataTE$V5) & dataTE$V4 == 0, ]$V5=0
colnames(dataTE)<-c("win_chr","win_s","win_e","TEvar","TEgini","TEmax")
data=merge(data, dataTE, by = c("win_chr","win_s","win_e"),all.x=TRUE)

datasim=fread("./COMPLEXFEATURE/00_vcf/sim_al.bed") ##试一下把Vscore换为序列相似性
colnames(datasim)<-c("win_chr","win_s","win_e","win_sim")
data=merge(data, datasim, by = c("win_chr","win_s","win_e"),all.x=TRUE)




# bubble=fread("../01_bubble/all.bubblenew")
# bubble=bubble[,c(1,2,3,4,5,6)]
# colnames(bubble)<-c("win_chr","win_s","win_e","simple","ins_bubble","super")
# bubble$win_s=bubble$win_s-1
# data=merge(data, bubble, by = c("win_chr","win_s","win_e"),all.x=TRUE)
# data$anno_vscore_na=0
# data$anno_simple_na=0
# data[is.na(data$Vscore), ]$anno_vscore_na=1
# data[is.na(data$simple), ]$anno_simple_na=1

# data[is.na(data$Vscore), ]$Vscore=0
a=data[is.na(data$Vscore), ]
b=data[data$entropy==0, ]
data=data[!(is.na(data$Vscore)),]
data=data[data$entropy!=0,]


fwrite(rbind(a,b), "noalign1or0_vcf.bed", row.names = FALSE,sep="\t")

# data[is.na(data$entropy), ]$entropy=0
fwrite(data[is.na(data$simple), ], "Lar_bubblen.bed", row.names = FALSE,sep="\t")
data[is.na(data$simple), ]$simple=23
data[is.na(data$super), ]$super=4
data[is.na(data$ins_bubble), ]$ins_bubble=4  ###这里的填补是通过复杂区域的中位数计算的 希望能够忽略这些区域的连续性


data=data[data$win_chr!="chrM",]
compleno<-data[!complete.cases(data), ]
cleaned_data <- na.omit(data)
cleaned_data =cleaned_data [,c("win_chr", "win_s", "win_e","entropy","win_sim","TEvar","TEgini","sd_max","sd_var")]
write.csv(cleaned_data, "allfeaturen", row.names = FALSE)



# for m in {1..22} X Y; do
#     bash ./allvcf/allfeature/parell.sh $m &
# done