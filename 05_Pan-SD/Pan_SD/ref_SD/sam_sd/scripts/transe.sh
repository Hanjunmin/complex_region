awk 'BEGIN { OFS="\t" }
{
    # 转换每一列为浮点数并格式化输出
    if ($2 ~ /e/) {
        $2 = sprintf("%.0f", $2)
    }
    if ($3 ~ /e/) {
        $3 = sprintf("%.0f", $3)
    }
   if ($5 ~ /e/) {        
       $5 = sprintf("%.0f", $5)    
    }    
    if ($6 ~ /e/) {        
        $6 = sprintf("%.0f", $6)    
    }
    print $0
}' < $1 |less -S |sort -k2,2n  >$2