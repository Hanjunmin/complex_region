infile=$1  ##all1.end.statistics.filt.in
rmfile=$2  ##rm.out
rmtrffile=$3 ##rm_trf.s.bed
out=$4
cat $rmfile | grep Satellite | cut -f 1-3 | bedtools sort -i - |  bedtools merge -i - | bedtools coverage -header -a <(cat $infile |tail -n +2) -b -  >end1.statistics.ascov.txt
less -S  end1.statistics.ascov.txt |awk '$10>=0.9 && $11>=1000  && $12>=0.5 && $16<=0.7{print $0}'  >end1.statistics.ascov.end
bedtools coverage -a <(less -S end1.statistics.ascov.end) -b  $rmtrffile|less -S |awk '{print $0,$19-$18}' |awk '$21>=500{print $0}'   >SDend1.txt
less -S SDend1.txt | awk '!($1==$4 && (($5<=$2 && $6>=$2)||($5<=$3 && $6>=$3)||($2<=$5 && $3>=$5)||($2<=$6 && $3>=$6)))' |sort |uniq >$out
rm end1.statistics.ascov.end
rm SDend1.txt
rm end1.statistics.ascov.txt