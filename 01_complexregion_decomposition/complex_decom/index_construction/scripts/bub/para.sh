filename=$1
gfab=$2
gfabase=$3
ref=$4
rm "WGS_bubble_APG.txt"
befst1=0
befen1=0
while read -r line; do
        echo $line
        chr=$(echo "$line" | cut  -f 1)
        start=$(echo "$line" | cut  -f 2)
        end=$(echo "$line" | cut  -f 3)
        if [ "$start" -gt "$befst1" ] && [ "$end" -lt "$befen1" ]; then
            echo -e "$chr\t$start\t$end\t$befst1\t$befen1\t"extend"" >> "WGS_bubble_APG.txt"
            continue
        fi
        $gfabase sub $gfab $ref"#0#"$chr":"$start"-"$end --range --view --cutpoints 1 --guess-ranges -o ${filename}.gfa
        st1=$(less -S ${filename}.gfa  |grep "CHM13v2#" |less -S  |cut -f8 |sed 's/,//g' |tr ":" "\t" |tr "-" "\t" |cut -f4 |sort -k1,1n |head -n 1)
        en1=$(less -S ${filename}.gfa  |grep "CHM13v2#" |less -S  |cut -f8 |sed 's/,//g' |tr ":" "\t" |tr "-" "\t" |cut -f5 |sort -k1,1n |tail -n 1)
        echo -e "$chr\t$start\t$end\t$st1\t$en1\t" >> "WGS_bubble_APG.txt"
        befst1=$st1
        befen1=$en1
done < $filename