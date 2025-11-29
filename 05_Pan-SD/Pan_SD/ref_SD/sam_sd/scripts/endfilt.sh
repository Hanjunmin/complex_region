in1=$3
cat $in1 |cut -f 1-12 |grep -v "refsou" >filt1
cat $1 | grep Satellite | cut -f 1-3 | bedtools sort -i - |  bedtools merge -i - | bedtools coverage -header -a <(cat filt1 |tail -n +2) -b -  >end1.statistics.ascov.txt
less -S  end1.statistics.ascov.txt |awk '$10>=0.9 && $11>=1000  && $12>=0.5 && $16<=0.7{print $0}'  >end1.statistics.ascov.end
bedtools coverage -a <(less -S end1.statistics.ascov.end) -b  $2|less -S |awk '{print $0,$19-$18}' |awk '$21>=500{print $0}'   >SDend1.txt
less -S SDend1.txt | awk '!($1==$4 && (($5<=$2 && $6>=$2)||($5<=$3 && $6>=$3)||($2<=$5 && $3>=$5)||($2<=$6 && $3>=$6)))' |sort |uniq >filt1.out

cat filt1.out |cut -f 1-12 |grep -v "refsou" |awk 'OFS="\t"{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12}' >filt1
cat $1 | grep Satellite | cut -f 1-3 | bedtools sort -i - |  bedtools merge -i - | bedtools coverage -header -a <(cat filt1 |tail -n +2) -b -  >end1.statistics.ascov.txt
less -S  end1.statistics.ascov.txt |awk '$10>=0.9 && $11>=1000  && $12>=0.5 && $16<=0.7{print $0}'  >end1.statistics.ascov.end
bedtools coverage -a <(less -S end1.statistics.ascov.end) -b  $2|less -S |awk '{print $0,$19-$18}' |awk '$21>=500{print $0}'   >SDend1.txt
less -S SDend1.txt | awk '!($1==$4 && (($5<=$2 && $6>=$2)||($5<=$3 && $6>=$3)||($2<=$5 && $3>=$5)||($2<=$6 && $3>=$6)))' |sort |uniq >filt1.out

#cat filt1.out |grep -v "refsou"|awk '{tmp=$1; $1=$4; $4=tmp; tmp=$2; $2=$5; $5=tmp; tmp=$3; $3=$6; $6=tmp; print $0}' inputfile
cat filt1.out |grep -v "refsou" |awk 'OFS="\t"{print $4,$5,$6,$1,$2,$3,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' >$4

####最后生成了20列。cluste_end.smk 的整合必须要输入20列