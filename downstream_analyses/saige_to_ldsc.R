library(data.table)
library(tidyverse)
library(stringr)

setwd("/array/psivakumar/Emer_CH/gwas/uk_analysis/")

all.res <- fread("/array/psivakumar/Emer_CH/gwas/uk_analysis/SAIGE_results/imputed_all_results_uk_only.csv")
hrc.ids <- fread("/data/kronos/NGS_Reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
sites.to.use <- fread("/data/kronos/NGS_Software/ldsc/w_hm3.noMHC.snplist/w_hm3.noMHC.snplist")

ukbb.gwas.paths <- list.files(path = "/array/psivakumar/Emer_CH/gwas/ukbb_sumstats", pattern = "imputed", recursive = T, full.names = T)
ukbb.gwas.res <- lapply(ukbb.gwas.paths, fread)
#ukbb.gwas.names <- c("anxiety", "bipolar", "depression", "epilepsy", "head_injury", "headache", "migraine", "nervous_breakdown", "schizophrenia")
ukbb.gwas.names <- str_split(ukbb.gwas.paths, pattern = "/") %>% sapply(`[[`, 7)

all.res$CHR <- as.character(all.res$CHR)
all.res.ids <- left_join(all.res, hrc.ids, by = c("CHR" = "#CHROM", "POS"))
sig.res.ids <- filter(all.res.ids, p.value < 1E-5)
all.res.ids.ldhub <- filter(all.res.ids, ID %in% sites.to.use$SNP, imputationInfo > 0.7) %>% rename(SNP = ID, A1 = Allele2, A2 = Allele1, OTHER_ID = SNPID, PVALUE = p.value, OTHER_PV = p.value.NA, MAF = AF_Allele2)

ukbbGetRSIDs <- function(df) {
  df <- separate(df, "variant", c("CHR", "POS", "REF", "ALT"), sep = ':', remove = F)
  df$POS <- as.integer(df$POS)
  df <- left_join(df, hrc.ids, by = c("CHR" = "#CHROM", "POS"))
}

ukbbToLDSC <- function(df) {
  df <- na.omit(df) %>% filter(ID %in% sites.to.use$SNP) %>% rename(SNP = ID, A1 = ALT.x, A2 = REF.x, PVALUE = pval, BETA = beta, MAF = minor_AF, N = n_complete_samples)
}

ukbb.gwas.ids <- lapply(ukbb.gwas.res, ukbbGetRSIDs)
ukbb.gwas.sig.ids <- lapply(ukbb.gwas.ids, filter, pval < 1E-5)
ukbb.gwas.ids.ldsc <- lapply(ukbb.gwas.ids, ukbbToLDSC)

for (i in 1:length(ukbb.gwas.ids)) {
  fwrite(ukbb.gwas.sig.ids[[i]], paste0("/array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/", ukbb.gwas.names[i], "_sig.txt", sep = ""), sep = "\t")
  fwrite(ukbb.gwas.ids[[i]], paste0("/array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/", ukbb.gwas.names[i], "_full.txt", sep = ""), sep = "\t")
  fwrite(ukbb.gwas.ids.ldsc[[i]], paste0("/array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/", ukbb.gwas.names[i], "_ldsc.txt", sep = ""), sep = " ")
}
