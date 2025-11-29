library(data.table)
library(ineq)
ref<-fread("REF.CN")
other<-fread("FINAL_OTHERvntr.csv")

ref[, 2:542] <- lapply(ref[, 2:542, drop = FALSE], as.numeric)
subdata=ref[,c(2:542)]
library(dplyr)
sub_df <- subdata
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
library(ineq)
trimmed_gini <- function(x, prop = 0.05) {
  x=x[!is.na(x)]
  lower <- quantile(x, probs = prop, na.rm = TRUE)
  upper <- quantile(x, probs = 1 - prop, na.rm = TRUE)
  x_trimmed <- x[x >= lower & x <= upper]
  if (length(x_trimmed) < 2) {
    NA
  } else {
    Gini(x_trimmed)
  }
}

trimsd<- apply(sub_df, 1, trimmed_var)
trimgini<- apply(sub_df, 1, trimmed_gini)
end=ref[,1]
end$trimsd=trimsd
end$trimgini=trimgini
end[is.na(trimsd), trimsd := 1]
ref_fea=end
##################Other

library(dplyr)
otherdf <- other %>% 
  select(
    any_of(c("chrom", "start_pos", "end_pos")),
    any_of(colnames(ref))
  )




otherdf[, 4:544] <- lapply(otherdf[, 4:544, drop = FALSE], as.numeric)
subdata=otherdf[,c(4:544)]
sub_df <- subdata
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
library(ineq)
trimmed_gini <- function(x, prop = 0.05) {
  x=x[!is.na(x)]
  lower <- quantile(x, probs = prop, na.rm = TRUE)
  upper <- quantile(x, probs = 1 - prop, na.rm = TRUE)
  x_trimmed <- x[x >= lower & x <= upper]
  if (length(x_trimmed) < 2) {
    NA
  } else {
    Gini(x_trimmed)
  }
}

trimsd<- apply(sub_df, 1, trimmed_var)
trimgini<- apply(sub_df, 1, trimmed_gini)
end=other[,c("chrom", "start_pos", "end_pos")]
end$trimsd=trimsd
end$trimgini=trimgini
end[is.na(trimsd), trimsd := 1]
other_fea=end

library(tidyr)
library(dplyr)
ref_fea <- ref_fea %>%
  separate(V1, into =c("chrom", "start_pos", "end_pos"), sep = "_")

end=rbind(ref_fea,other_fea)



fwrite(end,"VNTRall.fill.VNTR.sd",sep="\t")
