library(data.table)
library(tidyverse)

setwd("/array/psivakumar/Emer_CH/gwas/uk_analysis/")

ch.sig.res <- fread("ldsc/imputed_sig_results_e5_uk_rsIDs.tsv")
ukbb.gwas.sig.paths <- list.files(path = "/array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc", pattern = "*sig_sumstats.txt", full.names = T)
ukbb.gwas.sig.res <- lapply(ukbb.gwas.sig.res, fread)
ukbb.gwas.names <- c("anxiety", "bipolar", "depression", "epilepsy", "head_injury", "headache", "migraine", "nervous_breakdown", "schizophrenia")

ukbbToBed <- function(df) {
  df <- df %>% mutate(start = POS - 1, chrom = gsub("^", "chr", CHR)) %>% select(c(chrom, start, end = POS, ID))
}

ch.sig.bed <- ukbbToBed(ch.sig.res)
ukbb.sig.beds <- lapply(ukbb.gwas.sig.res, ukbbToBed)

fwrite(ch.sig.bed, "coloc-stats/ch_sig_res.bed", sep = '\t', col.names = F)
for (i in 1:length(ukbb.sig.beds)) {
  fwrite(ukbb.sig.beds[[i]], paste0("coloc-stats/ukbb_", ukbb.gwas.names[i], "_sig_res.bed", sep = ""), sep = "\t", col.names = F)
}
