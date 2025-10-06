
for i in {1..22} X Y;do
cd ./jmhan/PANSDEND/000_SD_sim/win/
python ./000_SD_sim/chm13_sim.py --c chr$i   --w chr$i".window"  --o chm13.chr${i}.sim
done

for i in {1..22} X Y;do
cd ./000_SD_sim/draft/chr${i}
cat ./000_SD_sim/winal/chr${i}/new.win.al |awk '$5-$4<=2000' >new.win.short &
done



for i in {1..22} X Y;do
cd chr${i}
split -l  50000 new.win.short short_split
for file in short_split*; do
    suffix="${file#short_split}" 
    echo "$suffix"
    echo "$file"  
    bedtools coverage -a chr${i}.vcfcor -b <(cat $file |cut -f 1,4,5 |awk 'OFS="\t"{print $1,$2-1000,$3+1000}' |awk 'OFS="\t"{$2=($2<0)?1:$2} 1') |nl -  |awk '$8>0'  |cut -f1 >"$suffix".vcfid
    awk 'NR==FNR {lines[$1]; next} FNR in lines'   "$suffix".vcfid ~/allvcf/chr${i}.vcf >"$suffix".vcf
    bedtools coverage -a <(cat chr${i}.posCHM.csv |awk 'OFS="\t"{print $7,$4,$5}' |tail -n +2) -b <(cat $file |cut -f 1,4,5 |awk 'OFS="\t"{print $1,$2-1000,$3+1000}' |awk 'OFS="\t"{$2=($2<0)?1:$2} 1') |nl -  |awk '$8>0'  |cut -f1  >"$suffix".posid
    awk 'NR==FNR {lines[$1]; next} FNR==1 || FNR in lines'  <(cat "$suffix".posid |awk '{print $1+1}' )  chr${i}.posCHM.csv >"$suffix".pos
done

for file in short_split*; do
  suffix="${file#short_split}" 
python ./script/seg1_correc_sim_addbord.py --c chr$i --vcf "$suffix".vcf  --w $file  --o $file".out" --seg "$suffix".pos
done
done

###integrate
python ./script/inte_short_split.py --c chr1

for i in {1..22} X Y;do
python ./script/inte_allsim.py --c "chr"${i}
done



while IFS= read -r sam; do
    samfile=$sam
    echo $samfile
    if [ "$sam" = "T2TCHM13" ]; then
    samfile="CHM13v2"
    cat "$sam".bed |cut -f2 >$samfile".out"
  else
    cat $samfile"_REFcor.pos" |cut -f2,9,10 |tail -n +2 >$samfile".temp"
    bedtools coverage -a <(cat "$sam".bed |tr ":" "\t" |tr "-" "\t" |tail -n +2 ) -b $samfile".temp"  |awk '{OFS="\t"}$4=="."{$4="0"}1' |awk '{OFS="\t"}$8==0{$4="."}1' |less -S  |awk 'OFS="\t"{print $4}' |awk -v sam=$sam 'BEGIN{print sam }1'  >$samfile".out"
	fi
    
done < allsamname 

paste <(cat ../chr$i"sim.out" |cut -f 543) *.out  >../chr$i"sim.end" 
done

###calculate the inter-seg homology and sd
library(dplyr)
library(data.table)

dfal <- data.frame()
main_folder_path <- '/PANSDEND/000_SD_sim/2k/'
file_paths=c()
for(i in c(1:22,"X","Y")) {
file_paths=c(file_paths,paste(main_folder_path,"chr",i,"sim.end",sep=""))
}


trimmed_var <- function(x, prop = 0.05) {
  x=x[!is.na(x)]
  lower <- quantile(x, probs = prop, na.rm = TRUE)
  upper <- quantile(x, probs = 1 - prop, na.rm = TRUE)
  x_trimmed <- x[x >= lower & x <= upper]
  if (length(x_trimmed) < 2) {
    NA
  } else {
    sd(x_trimmed)
  }
}

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)


registerDoParallel(cores = 20)
foreach(file_path = file_paths, .packages = c('data.table', 'dplyr')) %dopar% {
  print(file_path)
  df <- fread(file_path)
  win1 <- df$win
  df <- df %>% mutate(across(everything(), as.numeric))
  trimsd <- apply(df, 1, trimmed_var)
  df <- df %>% mutate(sdrow_max = pmax(!!!select(., -c("win")), na.rm = TRUE))
  df %>%
    mutate(win = win1, sdtrimsd = trimsd) %>%
    select(win, sdrow_max, sdtrimsd) %>%
    fwrite(paste0(file_path, ".feature"), sep = "\t")
}



