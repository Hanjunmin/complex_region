.project/ENDVNTR/
###1.预处理trf
.pangenome/WGSfea/VNTR/TRF_VAMOStest/ben/trf_process.vntr.str.sh
###2.1中输出的结果运行 生成INS OTHER REF
bash .project/ENDVNTR/1.run_vntr.sh ##其中的py脚本改了max和min 为only_pkl.py
###3.整理REF的结果 （allCN1bech.csv+VNTRall.fill） 
snakemake -s 2.expand_rm.smk --jobs 10
for i in {1..22} X Y;do
	cd .project/ENDVNTR/newcen/chr$i
	cat *.pkl.csv |grep 'CHM' |awk 'OFS="\t"{print $2,$6,$7,$3,$8}' >REF.specifi.motif &
done
###测试一下基于REF的是否和之前的拷贝数一致
cd .project/ENDVNTR/draft/00_REF_GRAPH-ASSEM_TEST/
cat .project/ENDVNTR/iter1test/*/REF.specifi.motif >GRAPH.REF
cat .project/ENDVNTR/STRVNTR130/* |grep -v "CHR" |awk 'OFS="\t"{print $1"_"$2"_"$3,$4,$5,$6,$7}'>ASSEM.REF
library(data.table)
data1<-fread("GRAPH.REF")
data2<-fread("ASSEM.REF")
data1<-data1[,c(1,4,5)]
data2<-data2[,c(1,2,3)]
colnames(data1)<-c("ID","GRAPH_CN","GRAPH_MOTIF")
colnames(data2)<-c("ID","ASSEM_CN","ASSEM_MOTIF")
data<-merge(data1,data2,by="ID",all.y=TRUE)
nrow(data)
nrow(data[data$GRAPH_CN!="LARRANGE",])
data=data[data$GRAPH_CN!="LARRANGE",]
data$minus=as.numeric(data$GRAPH_CN)/as.numeric(data$ASSEM_CN)
nrow(data[data$minus>2,])

###继续整理,把GRAPH中的T2TCHM13添加进去
cd .project/ENDVNTR/00_REF/
cat .project/ENDVNTR/newcen/*/REF.specifi.motif >T2Tmotifall.csv
cd .project/ENDVNTR/00_REF/allchormo
cp .project/ENDVNTR/newcen/*/*_CNref.csv ./
cd ..
import pandas as pd
import os
folder_path=".project/ENDVNTR/00_REF/allchormo"
intelis=[]
for file_name in os.listdir(folder_path):
    file_path = os.path.join(folder_path, file_name)
    simpair = pd.read_csv(file_path,sep="\t",header=0)
    intelis.append(simpair)

simpair=pd.concat(intelis, axis=0, ignore_index=True)
simpair_filled = simpair.fillna('NAN')
simpair_filled.columns = ['vntr'] + simpair_filled.columns[1:].tolist()

t2t = pd.read_csv("T2Tmotifall.csv",sep="\t",header=None)
t2t['chr'] = t2t[0].str.split('_').str[0]
t2t.columns=['vntr','start','end','cn','motif','chr']
t2t['min_motiflen'] = t2t['motif'].str.split(',').apply(lambda x: min(len(s) for s in x))
simpair_filledal = simpair_filled.merge(t2t, on="vntr", how="left")
filtered_df = simpair_filledal[simpair_filledal['cn']!="LARRANGE"]
filtered_df = filtered_df[filtered_df['cn']!=0]
filtered_df = filtered_df[filtered_df['cn']!='0']
filtered_df.to_csv('VNTRall.fill', index=False,sep="\t")  # 不保存索引列
t2tchm13=filtered_df[['CHM13v2','vntr','start','end','cn','motif','chr']]
t2tchm13.to_csv('allCN1bech.csv', index=False,sep="\t")  # 不保存索引列

###allCN1bech.csv 传到.pangenome/WGSfea/VNTR/TRF_VAMOStest/ben/ENDBEN/中进行benchmark加传回来VNTR

###看一下这些VNTR的相似性
for i in {1..22} X Y;do
cd .project/ENDVNTR/newcen/chr${i}
cat *pkl.csv |less -S  |grep -v "LARRANGE" |less -S  |cut -f 2,4,9|grep "chr" |sort -k1,1  >alsim &
done
#计算每个loci的相似性平均值


for i in {1..22} X Y;do
cd .project/ENDVNTR/newcen/chr${i}
Rscript -e '
library("data.table")
library(dplyr)
data<-fread("alsim")
result <- data %>% 
  group_by(V1) %>% 
  summarise(mean_simi = mean(V3, na.rm = TRUE))
fwrite(result,"sim.out",sep="\t") ' &
done
cd .project/ENDVNTR/newcen
cat */sim.out >allchr.sim.end

###4.整理OTHER的结果 先跑INS吧 因为INS匀出来一些OTHER

for m in {1..22} X Y; do
    cd .project/ENDVNTR/iternew1/"chr"$m/
    cat *other.copynumber.csv |less -S  |grep 'chr' |less -S |awk '$5!=1' |less -S |grep -v "INS" |less -S |cut -f 2-9,13 >other.in & ####这里的INS是那些参考基因组上没有拷贝数（甚至没有这个区域），但是其他样本多了拷贝数（denovo的）
    cat *other.copynumber.csv |less -S  |grep 'chr' |less -S |awk '$5!=1' |less -S |grep  "INS" |less -S |cut -f 2-15  >ins.in & ####这里的INS是那些参考基因组上没有拷贝数（甚至没有这个区域），但是其他样本多了拷贝数（denovo的）
done 

####python脚本处理一下motif motif需要在ref上只有1个motif或者没有motif
for m in {1..22} X Y; do
    cd .project/ENDVNTR/iternew1/"chr"$m
    python .project/panVNTR_10000/refineother.py &
done 

###按照样本进行分割以便计算具体位置坐标
cd .project/ENDVNTR/02_OTHER/
cat .project/ENDVNTR/iternew1/*/other.in.out >other.allchrom
mkdir SAM && cd SAM
library(data.table)
data<-fread("../other.allchrom")
data$sam <- sub("@.*", "", data$V7)
for (name in unique(data$sam)) {
  write.csv(data[data$sam==name,], file = name, row.names = FALSE)
}

cat .project/pansd/zdcopy/allsam | tr -d '\r' | xargs -I {} -P 50 python .project/ENDVNTR/3.find_pos_INSOTHER.py --S {} --type other
cat .project/pansd/zdcopy/allsam | tr -d '\r' | xargs -I {} -P 50 python .project/ENDVNTR/3.1find_pos_spe.py --sam {}

cd .project/ENDVNTR/02_OTHERnew/SAM/
cat *.END |grep 'chr' >ALL.specifi.cor

####INS中重新refine的other的加到这里
cat .project/ENDVNTR/01_INS/OTHER.all ALL.specifi.cor >OTHER.all



mkdir CHR && cd CHR
for i in {1..22} X Y; do
    less -S ../OTHER.all  |awk -v chr="chr$i" '$1==chr'  >"chr$i".other.vntr.end
    cat "chr$i".other.vntr.end |awk 'OFS="\t"{print $1":"$2"-"$3,$7,$10,$5,$6,$4":"$8"-"$9,$4,$8,$9,$11}'  >"chr$i".in
    less -S "chr$i".other.vntr.end |cut -f 1-3 |bedtools sort |uniq >"chr$i".temp
    bedtools intersect -a "chr$i".temp -b "chr$i".temp -f 0.7 -r -wa -wb |awk 'OFS="\t"{print $1":"$2"-"$3,$4":"$5"-"$6}' |less -S  |awk '$1!=$2' |sort |uniq >"chr"$i".pair.id"
done
cat ../OTHER.all |cut -f 1-3 |sort |uniq >other.vntr.all

###把other的对应的位置提出来重新看看有没有对应区域有motif存在 test一下

cat ../OTHER.all  |cut -f 1-3 |bedtools sort |uniq >can.reg

cd .project/ENDVNTR/02_OTHER/SAM/CHR/tes
for i in {1..22} X Y; do
    less -S ../can.reg  |awk -v chr="chr$i" '$1==chr'  >"chr$i".can.reg
    bedtools intersect -a "chr$i".can.reg -b .project/panVNTR/data/chr$i.posCHM.txt -wa -wb  > "chr$i".can.reg.ref.pos
done


cd .project/ENDVNTR/02_OTHER/SAM/CHR/tes
for i in {1..22} X Y; do
echo $i
less -S .project/ENDVNTR/02_OTHER/SAM/OTHER.all  |awk -v chr="chr$i" '$1==chr'  >"chr$i".other.all
done





##根据上面的pos去获取对应的不同样本的segment和对应的位置信息 没有相应pos的segment可以稍微扩一点找对应区域 如果相差大就不管了
bash .project/ENDVNTR/other_pos.sh
snakemake -s .project/ENDVNTR/02_OTHER/other.expand_rm.smk --jobs 10
##根据上述位置寻找对应区域的拷贝数等
##整合并合并
cd .project/ENDVNTR/02_OTHER/SAM/CHR/tes
for i in {1..22} X Y; do
	{
	cat "chr$i"_*pkl.csv |grep "chr" >"chr$i"OTHER.csv
    less -S "chr$i"OTHER.csv |cut -f2 |tr "_" "\t" |sort |uniq >"chr$i".temp
    bedtools intersect -a "chr$i".temp -b "chr$i".temp -f 0.7 -r -wa -wb |awk 'OFS="\t"{print $1"_"$2"_"$3,$4"_"$5"_"$6}' |less -S  |awk '$1!=$2' |sort |uniq >"chr"$i".pair.id"
} &
done




python .project/ENDVNTR/02_OTHER/new.inte_OTHER.py

cp FINAL_OTHERALL.csv  .project/ENDVNTR/02_OTHER/
cp FINAL_OTHERvntr.csv  .project/ENDVNTR/02_OTHER/
cp BENCH_OTHERvntr.csv  .project/ENDVNTR/02_OTHER/





###5.整理INS的结果 （INS不用扩！！！）
cd .project/ENDVNTR/01_INS/
cat .project/ENDVNTR/iternew1/*/ins.in >ins.allchrom
less -S ins.allchrom  |awk '$4>6' >ins.allchrom.filt

mkdir SAM && cd SAM
library(data.table)
data<-fread("../ins.allchrom.filt")
data$sam <- sub("@.*", "", data$V7)
for (name in unique(data$sam)) {
  write.csv(data[data$sam==name,], file = name, row.names = FALSE)
}
cat .project/pansd/zdcopy/allsam | tr -d '\r' | xargs -I {} -P 50 python .project/ENDVNTR/3.find_pos_INSOTHER.py --S {} --type INS ###生成segment对应的pos
cat .project/pansd/zdcopy/allsam | tr -d '\r' | xargs -I {} -P 50 python .project/ENDVNTR/3.2findspeINS_pos.py --sam {} ###上面sam.out中的pos生成具体的结果 smvntr.out  larvntr.refine 
cat .project/pansd/zdcopy/allsam | tr -d '\r' | xargs -I {} -P 50 python .project/ENDVNTR/4.larvntr_refine.py --SAM {} ###lar.refine.out
cat .project/pansd/zdcopy/allsam | tr -d '\r' | xargs -I {} -P 50 python .project/ENDVNTR/5.inte_ins_lars.py --SAM {} ###生成INS OTHER  整合

###在integrate之前需要获取一下INS和OTHER对应区域的信息  here
cd .project/ENDVNTR/01_INS/SAM
cat *smvntr.out |less -S  |cut -f 1,2,3,4,6,11-12 |grep -v "REFCHR" |sort |uniq   >smvntral
cat *larvntr.refine |less -S  |cut -f 1,2,3,4,6,11-12 |grep -v "REFCHR" |sort |uniq   >larvntral
cat smvntral larvntral |less -S  |sort |uniq >smlarvntral

###整合
cd .project/ENDVNTR/01_INS
cat ./SAM/*INS.END |grep -v "REFCHR" >INS.temp
cat ./SAM/*INS.END |grep "REFCHR" |head -n1 >INS.head
cat INS.head INS.temp  >INS.all
cat ./SAM/*OTHER.END |grep -v "REFCHR" >OTHER.all


mkdir CHR
cd CHR
for i in {1..22} X Y; do
    less -S ../INS.all  |awk -v chr="chr$i" '$1==chr'  >"chr$i".INS
done


###在integrate之前需要获取一下INS和OTHER对应区域的信息
for i in {1..22} X Y; do
    python .project/ENDVNTR/01_INS/INSinte.py --chr "chr$i" &
     # python .project/ENDVNTR/01_INS/INSinte2.py --chr "chr$i" &
done


for i in {1..22} X Y; do
    python .project/ENDVNTR/01_INS/INSinte2.py --chr "chr$i" &
done

###.PANSDEND/02_SAMSD/02_06-VALI/00_PAIR/align.al ###所有样本在图上是完整的边界的情况 .project/ENDVNTR/align.al
for i in {1..22} X Y; do
   less -S chr17.INS.inte.refine2.al |awk 'OFS="\t"{print $5"@"$6,$7,$8,$1,$2,$3,$4,$5,$9,$10,$11}' |grep "chr" |less -S  >x.temp
   
done

python .project/ENDVNTR/01_INS/new.inte_INS.py



# ####ins通过染色体进行合并
# cd .project/ENDVNTR/01_INS/
# # cat ./CHR/*INS.inte.out |grep -v "CHR" |bedtools sort |less -S >ins.all
# ###script
# import pandas as pd
# chromosomes = list(range(1, 23)) + ['X', 'Y']
# allsam = pd.read_csv(".project/ENDVNTR/allsam", sep="\t",header=None)
# columns = ["CHR","POS","MOTIF","MOTIFlen","CNmax","vali",'MOTIFsim']+list(allsam[0])
# dfal = pd.DataFrame(columns=columns)

# for chr_num in chromosomes:
# 	print(chr_num)
# 	ins = pd.read_csv("./CHR/chr"+str(chr_num)+".INS.inte.out", sep="\t",header=0)
# 	missing_columns = [col for col in allsam[0] if col not in ins.columns]
# 	for col in missing_columns:
# 		ins[col] = 'Miss'
# 	dfal = pd.concat([dfal, ins], ignore_index=True)

# dfal.to_csv('INS.END', index=False,sep="\t")

# ###只保留HPRC和HGSVC的样本
# import pandas as pd
# allsam = pd.read_csv("INS.END", sep="\t",header=0)
# df = allsam.loc[:, ~allsam.columns.str.startswith(('CHR', 'POS','MOTIF','MOTIFlen','CNmax','vali','C', 'K',))]
# all_no_rows = df.apply(lambda row: (row == 'Miss').all(), axis=1)
# no_rows_index = df[all_no_rows].index
# allsam.drop(no_rows_index).to_csv('BENCH_INSvntr.csv', index=False,sep="\t")                                                    

####################################################END########################################################
###1.REF:.project/ENDVNTR/00_REF/VNTRall.fill(allCN1bech.csv)
###2.INS:.project/ENDVNTR/01_INS/INS.END (BENCH_INSvntr.csv)
###3.OTHER:.project/ENDVNTR/02_OTHER/FINAL_OTHERvntr.csv
###这里的REF没有找VNTR的区域（VNTR+STR），其他两个直接统计的VNTR的区域


####################################################END最新版########################################################
###1.REF:.project/ENDVNTR/00_REF/VNTRall.fill(allCN1bech.csv)
###2.INS:.project/ENDVNTR/01_INS/CHR/FINAL_OTHERALL.csv|FINAL_INSvntr.csv(BENCH_INS_OTHERvntr.csv|BENCH_INS_INSvntr.csv)
###3.OTHER:.project/ENDVNTR/02_OTHER/SAM/CHR/tes/FINAL_OTHERvntr.csv(BENCH_OTHERvntr.csv)
###这里的REF没有找VNTR的区域（VNTR+STR），其他两个直接统计的VNTR的区域



bedtools subtract -a ins.all -b .project/panVNTR/data/cen_may.bed -A 











superpop= pd.read_csv("./DataPath-PublicAsm.csv", sep=",")
superpop.index=superpop['Assembly ID']
superpopdic=superpop['SP'].to_dict()
allsam=pd.DataFrame(allvntr).T.columns
popcor={}
for sam in allsam:
	ii=0
	for saminit in superpopdic.keys():
		if saminit in sam:
			ii=1
			print(saminit)
			popcor[sam]=superpopdic[saminit]
			break
	if ii==0:
		popcor[sam]="EAS"

superpop_df = pd.DataFrame(list(popcor.items()), columns=['sample', 'Superpopulation'])
df_transposed = vntrdf.T.reset_index()
df_transposed.columns=["sample"]+list(df_transposed.columns)[1:]
merged_df = pd.merge(df_transposed, superpop_df, on='sample')

alldf=pd.DataFrame()
alldf.index=merged_df.iloc[:, 1:-1].columns
for sup in set(merged_df['Superpopulation']):
	spdf=merged_df[merged_df['Superpopulation']==sup]
	numeric_df = spdf.iloc[:, 1:-1].apply(pd.to_numeric, errors='coerce')
	mean_values = numeric_df.mean()
	alldf[sup]=mean_values

alldf_filled = alldf.fillna(0)
alldf_normalized = alldf_filled.div(alldf_filled.sum(axis=1), axis=0)
merged_df = merged_df.sort_values(by='Superpopulation', ascending=True)
ll=merged_df.T
xxxx=pd.concat([ll,alldf_normalized], axis=1)
xxxx.to_csv('xxxx.csv', index=True,sep=",")

row_variances = alldf.var(axis=1)
cho=list(alldf[row_variances > 100].index)
cho.append('Superpopulation')
e=merged_df[list(cho)]
new.to_csv('chr18end.csv', index=True,sep="\t")





########################整合不同样本的T2Tsd数据
import os
import glob
import pandas as pd

# 获取当前工作目录
current_dir = os.getcwd()

# 查找当前目录下所有的 wchange.tsv 文件
tsv_files = glob.glob(os.path.join("./allW/", '*Wchange.tsv'))
edge=pd.read_csv("edge.bed", sep="\t",header=None)
alldf=pd.read_csv("./allW/CHM13.Wchange.tsv", sep="\t")
alldf=alldf[['SDregion']]
for tsv_file in tsv_files:
	file_name = os.path.basename(tsv_file)
	base_name = file_name.replace('.Wchange.tsv', '')
	print(base_name)
	samsd=pd.read_csv(tsv_file, sep="\t")
	samsd = samsd.drop_duplicates(subset=['SDregion'])
	samsd=samsd[['SDregion', 'Sim']]
	samsd.columns=['SDregion',base_name]
	alldf = pd.merge(alldf, samsd, on='SDregion', how='left')

mapped_values = list(map(lambda x: superpopdic.get(x, 'Not Found'), alldf.columns))



numeric_columns = alldf.select_dtypes(include=['float64', 'int64']).columns
alldf[numeric_columns] = alldf[numeric_columns].applymap(lambda x: 1 if x >= 0.95 else x)
alldf = alldf.fillna(-1)
alldf[numeric_columns] = alldf[numeric_columns].applymap(lambda x: 0 if x < 0.95 else x)

from collections import Counter
nodelistsize={}
nodelistpro={}
edgelistpro={}
edgelistother={}
for index, pair in edge.iterrows():
	print(index)
	pair1=pair[0]
	pair2=pair[1]
	pair1df=alldf[alldf['SDregion']==pair1]
	if len(pair1df)!=0:
		node1size=pair1df.iloc[0,1:].sum() ####
		sam1=pair1df.columns[pair1df.iloc[0] == 1].tolist() 
		pro1=sample_pro(sam1) #####
	else:
		node1size=0
		pro1=""
	pair2df=alldf[alldf['SDregion']==pair2]
	if len(pair2df)!=0:
		node2size=pair2df.iloc[0,1:].sum() ####
		sam2=pair2df.columns[pair2df.iloc[0] == 1].tolist()
		pro2=sample_pro(sam2) #####
	else:
		node2size=0
		pro2=""
	df=pd.concat([pair1df, pair2df], ignore_index=True)
	if len(df)==2:
		row_0 = df.iloc[0]
		row_1 = df.iloc[1]
		columns_both_one = df.columns[(row_0 == 1) & (row_1 == 1)] ####
		no=len(df.columns)-len(columns_both_one)
		columns_minis = df.columns[(row_0 == -1) | (row_1 == -1)]  #####
		noma=len(columns_minis)
		have=len(columns_both_one)
		edgepro=edge_pro(columns_both_one) #####
	else:
		no=len(df.columns)
		noma=no
		have=0
		edgepro=""
	nodelistsize[pair1]=node1size
	nodelistsize[pair2]=node2size
	nodelistpro[pair1]=pro1
	nodelistpro[pair2]=pro2
	edgelistpro[pair1+"|"+pair2]=edgepro
	edgelistother[pair1+"|"+pair2]=[no,noma,have]

nodelistsizedf=pd.DataFrame(list(nodelistsize.items()), columns=['Node', 'Size'])
nodelistprodf=pd.DataFrame(nodelistpro).T
edgelistprodf=pd.DataFrame(edgelistpro).T
edgelistotherdf=pd.DataFrame(edgelistother).T
edgelistotherdf.columns=['no', 'noma','have']
edgall=pd.merge(edgelistotherdf, edgelistprodf, left_index=True, right_index=True, how='left')
list_of_lists= edgall.index.str.split('|')
first_elements = [lst[0] for lst in list_of_lists]
sec_elements = [lst[1] for lst in list_of_lists]
edgall['A']=first_elements
edgall['B']=sec_elements
edgall['sorted_pair'] = edgall.apply(lambda row: tuple(sorted([row['A'], row['B']])), axis=1)
edgall_unique = edgall.drop_duplicates(subset='sorted_pair')
edgall_unique = edgall_unique.drop(columns='sorted_pair')


edgall_unique.to_csv('edgall.csv', index=False,sep=",")
nodelistsizedf.to_csv('nodesize.csv', index=False,sep=",")
nodelistprodf.to_csv('nodelistpro.csv', index=False,sep=",")


pd.DataFrame(list(edgelistother.items()), columns=['no', 'noma','have'])


superpop= pd.read_csv("./DataPath-PublicAsm.csv", sep=",")
superpop.index=superpop['Assembly ID']
superpopdic=superpop['SP'].to_dict()

xall={'AFR': 22, 'AMR': 17, 'EAS': 4, 'SAS': 1}


def sample_pro(samnow):
	mapped_values = list(map(lambda x: superpopdic.get(x, 'Not Found'), samnow))
	filtered_values = [value for value in mapped_values if value != 'Not Found']
	value_counts = Counter(filtered_values)
	total_count = sum(value_counts.values())
	value_proportions = {value: count / xall[value] for value, count in value_counts.items()}
	total_samples = sum(value_proportions.values())
	continent_proportions = {continent: count / total_samples for continent, count in value_proportions.items()}
	return continent_proportions

def edge_pro(samnow):
	mapped_values = list(map(lambda x: superpopdic.get(x, 'Not Found'), samnow))
	filtered_values = [value for value in mapped_values if value != 'Not Found']
	value_counts = Counter(filtered_values)
	return value_counts