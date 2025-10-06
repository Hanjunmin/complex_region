rsync -av .pangenome/Pan-SD/pansam/t2tmer.r .PANSDEND/02_SAMSD/02_00-REFSD/
rsync -av .pangenome/Pan-SD/pansam/merall.r .PANSDEND/02_SAMSD/02_00-REFSD/



while read line; do
 sam=$(echo "$line" | tr -d '\r')
 echo $sam
 cd .PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/${sam}
 cat process.txt end.txt |less -S |grep -v "sam" >sim.txt
 Rscript .PANSDEND/02_SAMSD/02_00-REFSD/t2tmer.r $sam
done < .DATA/APG/allsam

cd .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/
Rscript .PANSDEND/02_SAMSD/02_00-REFSD/merall.r 




.PANSDEND/02_SAMSD/02_00-REFSD/2_ALIGN/T2TPANSD/ ##最终的完整的基于REF的每个样本对应的区域


###先把上面的数据按照是否组装完整和是否是SDR等注释一下 下面是看人群stratified的区域的情况
import pandas as pd
data=".PANSDEND/02_SAMSD/02_00-REFSD/3_STA/eas_fisher_all.csv"
df=pd.read_csv(data,sep=",",header=0)
samlis=pd.read_csv("~/DATA/APG/allsam",sep="\t",header=None)
samlis=list(samlis[0])

for sam in samlis:
	indices = df[df[sam].isin(["NOSAM", "NO"])].index
	dfn=df.loc[indices,['A','B',sam]]
	dfn.to_csv('.PANSDEND/02_SAMSD/02_00-REFSD/3_STA/00_NOREGION/'+sam+"_no.sd", sep="\t", index=False, header=True)

cd .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/00_NOREGION/
cat .PANSDEND/DATA/allsam | tr -d '\r' | xargs -I {} -P 100 bash .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/anno_sdregion.para.sh {} ###生成SD区域的注释，是否是ASSEM的问题

import pandas as pd
def count_value(row):
    have_count = row.tolist().count('HAVE')
    poly_count = row.tolist().count('POLY')
    no_count = row.tolist().count('NO')
    nosam_count = row.tolist().count('NOSAM')
    AFh=have_count+poly_count
    AFal=AFh+no_count+nosam_count
    return pd.Series([AFh, AFal])

df=pd.read_csv(".PANSDEND/02_SAMSD/02_00-REFSD/3_STA/eas_fisher_all.csv",header=0,sep=",")
dfnew=df.copy()
samlis=pd.read_csv("~/DATA/APG/allsam",sep="\t",header=None)
samlis=list(samlis[0])
for sam in samlis:
	if sam=="CHM13v2":
			continue
	print(sam)
	dfn=df[['A','B',sam]]
	samsub=pd.read_csv("./00_NOREGION/"+sam+"_sd.assem.anno",header=None,sep="\t")
	anno_dict = pd.Series(samsub[1].values, index=samsub[0]).to_dict()
	dfn['anno1'] = dfn['A'].map(anno_dict)
	dfn['anno2'] = dfn['B'].map(anno_dict)
	dfn[sam]=dfn.apply(lambda row: f"{row['anno1']}_{row['anno2']}" if pd.notna(row['anno1']) or pd.notna(row['anno2']) else row[sam], axis=1)
	dfnew[sam]=dfn[sam]

dfnew[['AFhave', 'AFal']] = dfnew.apply(count_value, axis=1)
dfnew['AF']=dfnew['AFhave']/dfnew['AFal']

dfnew[['A','B','AF']].to_csv("eas_fisher.REF_refine.anno.af", sep="\t", index=False, header=True)
dfnew.to_csv("eas_fisher.REF_refine.anno", sep="\t", index=False, header=True)


less -S eas_fisher.REF_refine.anno.af  |awk '$3>=0.8' |tail -n +2 |awk '{print $1"\n"$2}' |tr ":" "\t" |tr "-" "\t" |less -S  |cal


.PANSDEND/02_SAMSD/02_00-REFSD/merall_add.r
cat eas_fisher_all_n.csv  |less -S  |tr "," "\t" |less -S |awk '$561>=2' |less -S  |grep -v "T2TSDDEL" |wc -l

cat eas_fisher_all_n.csv  |less -S  |tr "," "\t" |less -S |awk '$561>=2' |less -S |tail -n +2 |awk '{print $1"\n"$2}' |tr ":" "\t" |tr "-" "\t" |bedtools sort |bedtools merge  >pop_star.T2T.reg
cat eas_fisher_all_n.csv  |less -S  |tr "," "\t" |less -S |awk '$561<2' |less -S |tail -n +2 |awk '{print $1"\n"$2}' |tr ":" "\t" |tr "-" "\t"|bedtools sort|bedtools merge  >non-pop_star.T2T.reg
bedtools subtract -a pop_star.T2T.reg -b non-pop_star.T2T.reg |cal

cat eas_fisher_all_n.csv    |tr "," "\t" |tail -n +2 |less -S |awk '{print $1"\n"$2}' |tr ":" "\t" |tr "-" "\t"  >we.T2T
bedtools subtract -a  we.T2T -b .pangenome/acrocentric_region -A   >we.T2T.delacro


###根据图是否比对完整去注释 这里注释了组装和图不一致的区域以及inversion
cd .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/00_NOREGION/
cat .PANSDEND/DATA/allsam | tr -d '\r' | xargs -I {} -P 100 bash .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/graph_comple_align.sh {} ###生成SD区域的注释，是否是ASSEM的问题

###重新注释最后的表(根据图是否比对完整去注释)
cd .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/

import pandas as pd
def count_value(row):
    have_count = row.tolist().count('HAVE')
    poly_count = row.tolist().count('POLY')
    no_count = row.tolist().count('NO')
    nosam_count = row.tolist().count('NOSAM')
    AFh=have_count+poly_count
    AFal=AFh+no_count+nosam_count
    return pd.Series([AFh, AFal])

data="eas_fisher_all_n.csv"
df=pd.read_csv(data,header=0,sep=",")
dfnew=df.copy()
samlis=pd.read_csv("~/DATA/APG/allsam",sep="\t",header=None)
samlis=list(samlis[0])
for sam in samlis:
	if sam=="CHM13v2":
		continue
	print(sam)
	dfn=df[['A','B',sam]]
	samsub=pd.read_csv("./00_NOREGION/"+sam+".graph_err.anno",header=None,sep="\t")
	anno_dict = pd.Series(samsub[1].values, index=samsub[0]).to_dict()
	dfn['anno1'] = dfn['A'].map(anno_dict)
	dfn['anno2'] = dfn['B'].map(anno_dict)
	dfn[sam]=dfn.apply(lambda row: f"{row['anno1']}_{row['anno2']}" if pd.notna(row['anno1']) or pd.notna(row['anno2']) else row[sam], axis=1)
	dfnew[sam]=dfn[sam]

dfnew[['AFhave', 'AFal']] = dfnew.apply(count_value, axis=1)
dfnew['AF']=dfnew['AFhave']/dfnew['AFal']

dfnew[['A','B','AF']].to_csv("eas_fisher.REF_refine.graph.anno.af", sep="\t", index=False, header=True)
dfnew.to_csv("eas_fisher.REF_refine.graph.anno", sep="\t", index=False, header=True)



###将上面的结果进行基因注释 
cat .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/eas_fisher.REF_refine.graph.anno |tr "," "\t" |cut -f1 |tail -n +2 |tr ":" "\t" |tr "-" "\t" >A
cat .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/eas_fisher.REF_refine.graph.anno |tr "," "\t" |cut -f2 |tail -n +2 |tr ":" "\t" |tr "-" "\t" >B
bedtools intersect -a A -b .SDR/liftoffgene/protein_coding/T2TCHM13V2.0.gff -wa -wb  >A.inter.al
bedtools intersect -a B -b .SDR/liftoffgene/protein_coding/T2TCHM13V2.0.gff -wa -wb  >B.inter.al
bedtools intersect -a A -b .SDR/liftoffgene/protein_coding/T2TCHM13V2.0.gff>A.inter
bedtools intersect -a B -b .SDR/liftoffgene/protein_coding/T2TCHM13V2.0.gff>B.inter
paste  A.inter.al A.inter |grep -E "pseudogene|protein_coding" |less -S |awk 'OFS="\t"{print $1":"$2"-"$3,$7,$8,($11-$10)/($6-$5),($11-$10)/($3-$2),($6-$5)/($3-$2)}' | grep -v "LOC" >A.gff1
paste  B.inter.al B.inter |grep -E "pseudogene|protein_coding"  |less -S |awk 'OFS="\t"{print $1":"$2"-"$3,$7,$8,($11-$10)/($6-$5),($11-$10)/($3-$2),($6-$5)/($3-$2)}' | grep -v "LOC" >B.gff1

library(dplyr)
library(data.table)
library(stringdist)
data<-fread("eas_fisher.REF_refine.graph.anno")
Agff<-fread("A.gff1")
colnames(Agff)<-c("A","genenameA","genetypeA","ratioA_ge","ratioA_sd","geneA_sd")
Bgff<-fread("B.gff1")
colnames(Bgff)<-c("B","genenameB","genetypeB","ratioB_ge","ratioB_sd","geneB_sd")
data1<-merge(distinct(data),distinct(Agff),by="A",all.x=TRUE,allow.cartesian = TRUE)
data2<-merge(distinct(data1),distinct(Bgff),by="B",all.x=TRUE,allow.cartesian = TRUE)
cleaned_data <- na.omit(data2)
endall<-cleaned_data
endall <- endall %>%
  mutate(lev_similarity = 1-stringdist(genenameA, genenameB, method = "lv") /
           pmin(nchar(genenameA), nchar(genenameB)))

endall=endall[endall$lev_similarity>0,]
endall$first_letter_same <- apply(endall, 1, function(row) {
  first_char_A <- substr(row['genenameA'], 1, 1)  # 提取genenameA的第一个字母
  first_char_B <- substr(row['genenameB'], 1, 1)  # 提取genenameB的第一个字母

  # 判断两个字母是否相同
  return(first_char_A == first_char_B)
})

endall=endall[endall$first_letter_same=="TRUE",]
delpse=endall[endall$genetypeA=="protein_coding" |endall$genetypeB=="protein_coding",]
delpse1=delpse[delpse$eas_eur_p<0.05 |delpse$eas_amr_p<0.05 |delpse$eas_afr_p<0.05 |delpse$eas_sas_p<0.05,]
delpse2=delpse[delpse$eas_eur_p<0.05  &delpse$eas_amr_p<0.05 & delpse$eas_afr_p<0.05 &delpse$eas_sas_p<0.05,]

fwrite(delpse1,"eas_fisher_delinvpse.csv",sep=",")

less -S eas_fisher_delinvpse.csv |head -n 1 >header
less -S eas_fisher_delinvpse.csv |grep -v "ASSEM" |grep -v "GRA"  |grep -v "POLY" |less -S |grep  "NOSAM" |less -S >cho ###这些可以去试一下
cat header cho >cho.al


library(data.table)
data<-fread("cho.al")
data$new_column <- apply(data, 1, function(row) {
  paste(colnames(data)[row == "NOSAM"], collapse = ";")
})
data$new_column1 <- apply(data, 1, function(row) {
  paste(colnames(data)[row == "NO"], collapse = ";")
})
data$nosam <- paste(data$new_column, data$new_column1, sep = ";")

genelis=unique(c(data$genenameA,data$genenameB))
fwrite(as.data.frame(genelis),"genelis.csv")
fwrite(data[,c("A","B","nosam","count_lt_0.05")],"cho.al.sam.csv",sep="\t")


###对这些区域画两个图,一个是样本的svbyeye的图，一个是SD对应的基因的图
##
cd  .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/01_para_gene/
Rscript .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/gene_gene.r ##生成所有pair两个region对应的gene注释，但是看起来好像只有proteincoding的

cd .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/02_para_svbyeye/
# 使用制表符作为分隔符
cat cho.al.sam.csv |tail -n +2 |sort |uniq >cho.al.sam.csv.temp
while IFS=$'\t' read -r col1 col2 col3
do
A=$col1
B=$col2
sam=$(echo $col3|tr ";" "\t" |cut -f1)
file2=.PANSDEND/02_SAMSD/02_00-REFSD/3_STA/02_para_svbyeye/$A"@"$B/C_myplot.png
if [[ -f "$file2" ]]; then
      echo "yes"
       continue
else
  echo "no"
  bash .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/svbyeye.sh $A $B $sam &

fi

done < cho.al.sam.csv.temp

###svbyeye的图和基于基因的图放一起
# cd .PANSDEND/02_SAMSD/02_00-REFSD/3_STA
# while IFS=$'\t' read -r col1 col2 col3
# do
# A=$col1
# B=$col2
# convert $file1 $file2 -resize 2000x -append .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/03_mer_gene_svby/$A"@"$B.png  

# done < cho.al.sam.csv.temp




while IFS=$'\t' read -r col1 col2 col3 col4
do
A=$col1
B=$col2
if [ $col4 -ge 2 ]; then ##只可视化一下大于两个人群的吧
  file1=.PANSDEND/02_SAMSD/02_00-REFSD/3_STA/01_para_gene/$A"@"$B.png
	file2=.PANSDEND/02_SAMSD/02_00-REFSD/3_STA/02_para_svbyeye/$A"@"$B/C_myplot.png

  if [[ -f "$file1" && -f "$file2" ]]; then
          echo "处理: $A@$B - 两个文件都存在"
          convert $file1 $file2 -resize 2000x -append .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/03_mer_gene_svby/$A"@"$B.png
  else
          echo "错误: convert命令执行失败 - $A@$B"
  fi
fi

if [ $col4 -ge 3 ]; then ##只可视化一下大于两个人群的吧
 file1=.PANSDEND/02_SAMSD/02_00-REFSD/3_STA/01_para_gene/$A"@"$B.png
	file2=.PANSDEND/02_SAMSD/02_00-REFSD/3_STA/02_para_svbyeye/$A"@"$B/C_myplot.png
 if [[ -f "$file1" && -f "$file2" ]]; then
          echo "处理: $A@$B - 两个文件都存在"
          convert $file1 $file2 -resize 2000x -append .PANSDEND/02_SAMSD/02_00-REFSD/3_STA/04_mer_3/$A"@"$B.png
  else
          echo "错误: convert命令执行失败 - $A@$B"
  fi
fi

done < cho.al.sam.csv.temp


.PANSDEND/02_SAMSD/02_00-REFSD/4_case/