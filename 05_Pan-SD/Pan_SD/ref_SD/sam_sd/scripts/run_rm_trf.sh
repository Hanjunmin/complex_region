  rmtrfin=$1
  line=$2
  RMsoft=$3
  script=$4
  less -S $rmtrfin |bedtools sort |bedtools merge >rm_in.mer
  awk 'OFS="\t"{if ($2 == 0) $2 = 1}1' rm_in.mer > temp && mv temp rm_in.mer
  [ -f rm_in.fa ] && rm rm_in.fa
  fa=$line.fa
  while IFS=$'\t' read -r col1 col2 col3
  do
  samtools faidx  $fa "$col1":"$col2"-"$col3">>rm_in.fa
  done < rm_in.mer
  ## 换名字跑rm和trf
  samtools faidx rm_in.fa
  nl rm_in.fa.fai |awk 'OFS="\t"{print $2,"a"$1}' >newid
  sed -i  's/\./_/g;s/\-/_/g; s/:/_/g; s/|/_/g'  newid
  sed -i  's/\./_/g;s/\-/_/g; s/:/_/g; s/|/_/g'  rm_in.fa
  seqkit replace --ignore-case --kv-file newid --pattern "^(\w+)" --replacement "{kv}" rm_in.fa -o rm_in.new.fa
  RepeatMasker -s -xsmall -e ncbi -species human -dir otherRM rm_in.new.fa -pa 20
  RM2Bed.py -d otherRM otherRM/rm_in.new.fa.out
  trf rm_in.new.fa  2 7 7 80 10 50 500 -h -l 20 -ngs >rm_in.new.dat
  #python /dssg/home/acct-clsmyf/clsmyf-user1/pangenome/Pan-SD/pansam/datafrommao/trf2bed.py --i rm_in.new.dat --o rm_in.new.trf.bed
  python ${script}trf2bed.py --i rm_in.new.dat --o rm_in.new.trf.bed
  
  # less -S otherRM/rm_in.new.fa_rm.bed
  ##下面是把rm和trf对应的区域把名字进行还原
  paste <(cat rm_in.fa.fai |cut -f1) <(cat newid) |less -S |tr ":" "\t" |less -S >allcor.id


  Rscript -e 'library(data.table)
  rmdf<-fread("./otherRM/rm_in.new.fa_rm.bed",header=FALSE)
  trf<-fread("rm_in.new.trf.bed",header=TRUE)
id=fread("allcor.id",header=FALSE)
  colnames(id)=c("sam","s_e","name","chaname")
 
if (nrow(rmdf) != 0) {
    rmdf<-rmdf[,c(1,2,3,4,7,8)]
    colnames(rmdf)=c("chaname","s","e","r1","r2","r3")
    end<-merge(rmdf,id,by="chaname",all.x = TRUE)
  end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), `[`, 1)
  end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), function(x) as.integer(x[1]))
  end$ns=end$s+end$sinit
  end$ne=end$e+end$sinit
  end<-end[,c("sam","ns","ne","r1","r2","r3")]
  fwrite(end, "rm.out", sep = "\t", col.names = FALSE)
}else{file.create("rm.out")}

if (nrow(trf) != 0) {
    print("trf is empty")
    trf<-trf[,c(1,2,3)]
    colnames(trf)=c("chaname","s","e")
    end<-merge(trf,id,by="chaname",all.x = TRUE)
  end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), `[`, 1)
  end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), function(x) as.integer(x[1]))
  end$ns=end$s+end$sinit
  end$ne=end$e+end$sinit
  end<-end[,c("sam","ns","ne")]
  fwrite(end, "trf.out", sep = "\t", col.names = FALSE)
}else{file.create("trf.out")}
'
  cat trf.out <(cat rm.out|cut -f 1-3) >rm_trf.s.bed


  



  # less -S rm_in |bedtools sort |bedtools merge >rm_in.mer
  # [ -f rm_in.fa ] && rm rm_in.fa
  # fa=$line.fa
  # while IFS=$'\t' read -r col1 col2 col3
  # do
  # samtools faidx  $fa "$col1":"$col2"-"$col3">>rm_in.fa
  # done < rm_in.mer
  # ## 换名字跑rm和trf
  # samtools faidx rm_in.fa
  # nl rm_in.fa.fai |awk 'OFS="\t"{print $2,"a"$1}' >newid
  # sed -i  's/\./_/g;s/\-/_/g; s/:/_/g; s/|/_/g'  newid
  # sed -i  's/\./_/g;s/\-/_/g; s/:/_/g; s/|/_/g'  rm_in.fa
  # seqkit replace --ignore-case --kv-file newid --pattern "^(\w+)" --replacement "{kv}" rm_in.fa -o rm_in.new.fa
  # RepeatMasker -s -xsmall -e ncbi -species human -dir otherRM rm_in.new.fa -pa 20
  # /home/Public/software/RepeatMasker-4.1.4/util/RM2Bed.py -d otherRM otherRM/rm_in.new.fa.out
  # trf rm_in.new.fa  2 7 7 80 10 50 500 -h -l 20 -ngs >rm_in.new.dat
  # python /home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/trf2bed.py --i rm_in.new.dat --o rm_in.new.trf.bed
  # less -S otherRM/rm_in.new.fa_rm.bed
  # ##下面是把rm和trf对应的区域把名字进行还原
  # paste <(cat rm_in.fa.fai |cut -f1) <(cat newid) |less -S |tr ":" "\t" |less -S >allcor.id
  # # library(data.table)
  # # rmdf<-fread("./otherRM/rm_in.new.fa_rm.bed",header=FALSE)
  # # trf<-fread("rm_in.new.trf.bed",header=TRUE)
  # # rmdf<-rmdf[,c(1,2,3,4,7,8)]
  # # trf<-trf[,c(1,2,3)]
  # # colnames(rmdf)=c("chaname","s","e","r1","r2","r3")
  # # colnames(trf)=c("chaname","s","e")
  # # id=fread("allcor.id",header=FALSE)
  # # colnames(id)=c("sam","s_e","name","chaname")
  # # end<-merge(rmdf,id,by="chaname",all.x = TRUE)
  # # end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), `[`, 1)
  # # end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), function(x) as.integer(x[1]))
  # # end$ns=end$s+end$sinit
  # # end$ne=end$e+end$sinit
  # # end<-end[,c("sam","ns","ne","r1","r2","r3")]
  # # fwrite(end, "rm.out", sep = "\t", col.names = FALSE)
  # # end<-merge(trf,id,by="chaname",all.x = TRUE)
  # # end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), `[`, 1)
  # # end$sinit <- sapply(strsplit(as.character(end$s_e), "-"), function(x) as.integer(x[1]))
  # # end$ns=end$s+end$sinit
  # # end$ne=end$e+end$sinit
  # # end<-end[,c("sam","ns","ne")]
  # # fwrite(end, "trf.out", sep = "\t", col.names = FALSE)
  # cat trf.out <(cat rm.out|cut -f 1-3) >rm_trf.s.bed
  # cat T2Tpolymor1other.bed |cut -f 1-12 >filin.bed
  # bash /home/jmhan/pangenome/WGSfea/allW/APG3allW/scripts/endfilt.sh rm.out rm_trf.s.bed filin.bed filin.bedend