


for i in {1..22} X Y;do
cd ./project/ENDTE/01_INTE/chr${i}/
for j in {0..9}; do
cat split_${j}.fa |grep ">" >split_${j}.temid 
done
done


python ./project/ENDTE/inte_seq.py  &




for i in {1..22} X Y;do
echo $i
cd ./project/ENDTE/01_INTE/chr${i}/
cat *rm.out  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown">chr${i}.end.out ##有repeats的SDR
cat *fa_rm.alseq >chr${i}.end.alSDR ###全部的SDR
cat ./project/ENDTE/WGS.window.bed  |awk -v chr=chr${i} '$1==chr' >chr${i}.window

done



for i in {1..22} X Y;do
echo $i
cd ./project/ENDTE/01_INTE/chr${i}/

less -S chr${i}.end.alSDR   |cut -f 1-3 |sort |uniq >ref.chr${i}.anchor
less -S chr${i}.end.alSDR   |cut -f 1-4 |sort |uniq >ref.chr${i}.anchor.sam
cp ./project/ENDTE/00_INS/chr${i}/DEL.in ./
cat DEL.in   |tr "@" "\t" |awk -v chr=chr${i} 'OFS="\t"{print chr,$8,$9,$2}'  |sed "s/\.0//g"  >DEL.in.sam
cat DEL.in.sam |cut -f 1-3 >DEL.in
bedtools intersect -b DEL.in.sam -a ./data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed -wa -wb  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.chr${i}.anchor.rm.del
bedtools intersect -b ref.chr${i}.anchor -a ./data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed -wa -wb -f 0.2  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.chr${i}.anchor.rm
cat DEL.in ref.chr${i}.anchor |sort |uniq >ref.chr${i}.anchoral
bedtools intersect -b ref.chr${i}.anchoral -a ./data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed -wa -wb -f 0.2  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.chr${i}.anchoral.rm
cat DEL.in.sam ref.chr${i}.anchor.sam >ref.chr${i}.anchor.samal

done  ###ref.chr${i}.anchoral.rm所有的SDR对应的repeats

#####后面根据 DEL.in和ref.chr${i}.anchor 一起与win.txt取交集 python中获取


for i in {1..22} X Y;do
echo $i
cd ./project/ENDTE/01_INTE/chr${i}/
#bedtools intersect -b ref.chr${i}.anchor -a ./data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed -wa -wb -f 0.2  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.chr${i}.anchor.rm
###取0.2是希望那些INS的原始的数据不要被过滤掉
Rscript -e '
library(data.table)
args <- commandArgs(trailingOnly=TRUE)
chr=args[1]
data<-fread(paste(chr,".end.out",sep=""))
data1=fread(paste("ref.",chr,".anchor.samal",sep=""))
samlis=unique(data$V4)
samlisdf=as.data.frame(samlis)
fwrite(samlisdf,"sam.list",sep="\t",col.names=FALSE)
for(sam in samlis){
	print(sam)
samdf=data[data$V4==sam] ###SDR区域query的repeats组成
samdf1=data1[data1$V4==sam] ###SDR区域坐标
fwrite(samdf,paste(sam,".out",sep=""),sep="\t",col.names=FALSE)
fwrite(samdf1,paste(sam,".outpos",sep=""),sep="\t",col.names=FALSE)

common=paste("bedtools intersect -a ",chr,".window  -b ",sam,".out"," -wa -wb >",sam,".win.out",sep="")
system(common) 
common=paste("bedtools intersect -a ",chr,".window  -b ",sam,".outpos"," -wa -wb >",sam,".win.outpos",sep="")
system(common) 
common=paste("bedtools subtract -a ",chr,".window  -b ",sam,".outpos"," -A >",sam,".win.nointer.outpos",sep="")
system(common) 
common=paste("python ./project/ENDTE/mergeSDR.py --sam ",sam,sep="")
system(common) 

common=paste("bedtools intersect -a ",sam,".win.txt  -b ",sam,".out"," -wa -wb >",sam,".win.outn",sep="")
system(common) 
common=paste("bedtools intersect -a ",sam,".win.txt  -b ",sam,".outpos"," -wa -wb >",sam,".win.outposn",sep="")
system(common) 

}' chr${i} &
done

 
for i in {1..22} X Y;do
echo $i
cd ./project/ENDTE/01_INTE/chr${i}/
bedtools intersect -a chr${i}.window -b ./data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed -wa -wb  |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown" >ref.chr${i}.rm &
done

####generate
####each sample  generate 3 tables

'''
        chm13_LINE1   chm13_SINE    INS_LINE1   INS_SINE   SDR_LINE1  SDR_SINE
window       14          12           3				1        3          5        14-3+3



'''


###py

#### merge SDR的窗口 *win.txt
####  1.生成一个所有refwindow对应的repeats种类的表
for i in {1..22} X Y;do
cd ./project/ENDTE/01_INTE/chr${i}
T2TCHM13rm="./data/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed"
bedtools intersect -a <(cat *win.txt chr${i}.window |sort |uniq) -b $T2TCHM13rm -wa -wb |grep -v "Low_complexity" |grep -v "Satellite" |grep -v "Simple_repeat" |grep -v "Unknown">win.al.rm & ## win+rm
done
####  2.SDR对应的ref区域的repeats种类
ref.chr1.anchoral.rm ##rm+最后几列为SDRcoordinate

####  3.SDR对应的query区域的repeats种类
C001-CHA-E01-Mat.out ##SDRcoordinate+后面的rm

####  4.SDR对应的query区域的repeats种类以及对应的窗口（merge后的新窗口）
C046-CHA-N06-Mat.win.outn ##SDRcoordinate+后面的rm

#### 5.所有SDR区域对应于merge之后的新窗口
C046-CHA-N06-Mat.win.outposn ###ref.chr${i}.anchoral.rm  这个为对应的REF的repeats 配套使用


for i in {1..22} X Y;do
cd ./project/ENDTE/01_INTE/chr${i}
python ./project/ENDTE/inte_allsam.py --chr chr${i}
done


for i in {1..22} X Y;do
cd ./project/ENDTE/01_INTE/chr${i}
python ./project/ENDTE/inte_sam2.py --chr chr${i} &
done

###总结每个窗口的具体情况：
cd /home/jmhan/PANSDEND/02_SAMSD/02_03-VALI/00_PAIR
cat *cor.pos  |awk 'OFS="\t"{print $2,$9,$10,$5}' |grep -v "SAM" >allsampos
#./project/ENDTE/01_INTE/allsampos

for i in {1..22} X Y;do
echo $i
cd ./project/ENDTE/01_INTE/chr${i}/
cat TEal.txt |tail -n +2 |cut -f 1-3 >TE.head

Rscript ./project/ENDTE/filt_dot_cal.r &
done


cd ./project/ENDTE/01_INTE
cat */TEEND |grep -v "row_max" >TE.ALLchr ###传到我们服务器上


