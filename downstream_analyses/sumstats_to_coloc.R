library(snpStats)
library(coloc)
library(data.table)
library(tidyverse)

setwd("/array/psivakumar/Emer_CH/gwas/uk_swedish_alt/")

ch.sumstats <- fread("ldsc/uk_swedish_alt_all_rsids.tsv")
mig.sumstats <- fread("ldsc/ukbb_ukbb_migraine_sumstats_full_sumstats.txt")
sch.sumstats <- fread("ldsc/ukbb_ukbb_schizophrenia_sumstats_full_sumstats.txt")
#ldsnps <- fread("FUMA/ld.txt")
#chr12.ldsnps <- fread("FUMA/chr12_ld_snps.txt") %>% filter(R2 > 0.6, grepl("rs", RS_Number))

ch.sig.regions <- data.frame("CHR" = c(2,2,6,12,1), "start" = c(112656908,200360782,96854444,57497005,222098690), "end" = c(112785237,200513620,97067028,57614217,222227618))
ch.2.112.sumstats <- filter(ch.sumstats, CHR == ch.sig.regions[[1, 1]], POS >= ch.sig.regions[[1, 2]] & POS <= ch.sig.regions[[1, 3]])#, SNPID %in% ldsnps$SNP2)
ch.2.200.sumstats <- filter(ch.sumstats, CHR == ch.sig.regions[[2, 1]], POS >= ch.sig.regions[[2, 2]] & POS <= ch.sig.regions[[2, 3]])#, SNPID %in% ldsnps$SNP2)
ch.6.96.sumstats <- filter(ch.sumstats, CHR == ch.sig.regions[[3, 1]], POS >= ch.sig.regions[[3, 2]] & POS <= ch.sig.regions[[3, 3]])#, SNPID %in% ldsnps$SNP2)
ch.12.57.sumstats <- filter(ch.sumstats, CHR == ch.sig.regions[[4, 1]], POS >= ch.sig.regions[[4, 2]] & POS <= ch.sig.regions[[4, 3]])#, ID %in% chr12.ldsnps$RS_Number)
ch.1.222.sumstats <- filter(ch.sumstats, CHR == ch.sig.regions[[5, 1]], POS >= ch.sig.regions[[5, 2]] & POS <= ch.sig.regions[[5, 3]])#, ID %in% chr12.ldsnps$RS_Number)
#ch.sig.sumstats <- filter(ch.sumstats, CHR == 6)

#common.snps <- intersect(ch.sig.sumstats$ID, mig.sumstats$ID)
ch.2.112.mig.df <- inner_join(ch.2.112.sumstats, mig.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.2.200.mig.df <- inner_join(ch.2.200.sumstats, mig.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.6.96.mig.df <- inner_join(ch.6.96.sumstats, mig.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.12.57.mig.df <- inner_join(ch.12.57.sumstats, mig.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.1.222.mig.df <- inner_join(ch.1.222.sumstats, mig.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)

fwrite(ch.2.112.mig.df, "coloc/ch_mig_chr2_112_sumstats.csv")
fwrite(ch.2.200.mig.df, "coloc/ch_mig_chr2_200_sumstats.csv")
fwrite(ch.6.96.mig.df, "coloc/ch_mig_chr6_96_sumstats.csv")
fwrite(ch.12.57.mig.df, "coloc/ch_mig_chr12_57_sumstats.csv")
fwrite(ch.1.222.mig.df, "coloc/ch_mig_chr1_222_sumstats.csv")

#ch.mig.abf <- coloc.abf(dataset1 = list(beta = ch.mig.df$BETA, varbeta = ch.mig.df$SE, type = "cc", N = ch.mig.df$N, s = 0.142),
#                        dataset2 = list(beta = ch.mig.df$beta, varbeta = ch.mig.df$se, type = "cc", N = ch.mig.df$n_complete_samples, s = 0.029))

ch.2.112.mig.abf <- coloc.abf(dataset1 = list(pvalues = ch.2.112.mig.df$p.value, N = ch.2.112.mig.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.2.112.mig.df$pval, N = ch.2.112.mig.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.2.112.mig.df$minor_AF)
ch.2.200.mig.abf <- coloc.abf(dataset1 = list(pvalues = ch.2.200.mig.df$p.value, N = ch.2.200.mig.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.2.200.mig.df$pval, N = ch.2.200.mig.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.2.200.mig.df$minor_AF)
ch.6.96.mig.abf <- coloc.abf(dataset1 = list(pvalues = ch.6.96.mig.df$p.value, N = ch.6.96.mig.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.6.96.mig.df$pval, N = ch.6.96.mig.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.6.96.mig.df$minor_AF)
ch.12.57.mig.abf <- coloc.abf(dataset1 = list(pvalues = ch.12.57.mig.df$p.value, N = ch.12.57.mig.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.12.57.mig.df$pval, N = ch.12.57.mig.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.12.57.mig.df$minor_AF)
ch.1.222.mig.abf <- coloc.abf(dataset1 = list(pvalues = ch.1.222.mig.df$p.value, N = ch.1.222.mig.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.1.222.mig.df$pval, N = ch.1.222.mig.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.1.222.mig.df$minor_AF)


fwrite(as.data.frame(ch.2.112.mig.abf[1]), "coloc/coloc_abf_mig_chr2_112_res.csv", row.names = T)
fwrite(as.data.frame(ch.2.200.mig.abf[1]), "coloc/coloc_abf_mig_chr2_200_res.csv", row.names = T)
fwrite(as.data.frame(ch.6.96.mig.abf[1]), "coloc/coloc_abf_mig_chr6_96_res.csv", row.names = T)
fwrite(as.data.frame(ch.12.57.mig.abf[1]), "coloc/coloc_abf_mig_chr12_57_res.csv", row.names = T)
fwrite(as.data.frame(ch.1.222.mig.abf[1]), "coloc/coloc_abf_mig_chr1_222_res.csv", row.names = T)


ch.2.112.sch.df <- inner_join(ch.2.112.sumstats, sch.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.2.200.sch.df <- inner_join(ch.2.200.sumstats, sch.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.6.96.sch.df <- inner_join(ch.6.96.sumstats, sch.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.12.57.sch.df <- inner_join(ch.12.57.sumstats, sch.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)
ch.1.222.sch.df <- inner_join(ch.1.222.sumstats, sch.sumstats, by = "ID") %>% distinct(ID, .keep_all = T)

fwrite(ch.2.112.sch.df, "coloc/ch_sch_chr2_112_sumstats.csv")
fwrite(ch.2.200.sch.df, "coloc/ch_sch_chr2_200_sumstats.csv")
fwrite(ch.6.96.sch.df, "coloc/ch_sch_chr6_96_sumstats.csv")
fwrite(ch.12.57.sch.df, "coloc/ch_sch_chr12_57_sumstats.csv")
fwrite(ch.1.222.sch.df, "coloc/ch_sch_chr1_222_sumstats.csv")

#ch.mig.abf <- coloc.abf(dataset1 = list(beta = ch.mig.df$BETA, varbeta = ch.mig.df$SE, type = "cc", N = ch.mig.df$N, s = 0.142),
#                        dataset2 = list(beta = ch.mig.df$beta, varbeta = ch.mig.df$se, type = "cc", N = ch.mig.df$n_complete_samples, s = 0.029))

ch.2.112.sch.abf <- coloc.abf(dataset1 = list(pvalues = ch.2.112.sch.df$p.value, N = ch.2.112.sch.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.2.112.sch.df$pval, N = ch.2.112.sch.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.2.112.sch.df$minor_AF)
ch.2.200.sch.abf <- coloc.abf(dataset1 = list(pvalues = ch.2.200.sch.df$p.value, N = ch.2.200.sch.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.2.200.sch.df$pval, N = ch.2.200.sch.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.2.200.sch.df$minor_AF)
ch.6.96.sch.abf <- coloc.abf(dataset1 = list(pvalues = ch.6.96.sch.df$p.value, N = ch.6.96.sch.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.6.96.sch.df$pval, N = ch.6.96.sch.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.6.96.sch.df$minor_AF)
ch.12.57.sch.abf <- coloc.abf(dataset1 = list(pvalues = ch.12.57.sch.df$p.value, N = ch.12.57.sch.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.12.57.sch.df$pval, N = ch.12.57.sch.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.12.57.sch.df$minor_AF)
ch.1.222.sch.abf <- coloc.abf(dataset1 = list(pvalues = ch.1.222.sch.df$p.value, N = ch.1.222.sch.df$N, type = "cc", s = 0.142),
                        dataset2 = list(pvalues = ch.1.222.sch.df$pval, N = ch.1.222.sch.df$n_complete_samples, type = "cc", s = 0.029), MAF = ch.1.222.sch.df$minor_AF)


fwrite(as.data.frame(ch.2.112.sch.abf[1]), "coloc/coloc_abf_sch_chr2_112_res.csv", row.names = T)
fwrite(as.data.frame(ch.2.200.sch.abf[1]), "coloc/coloc_abf_sch_chr2_200_res.csv", row.names = T)
fwrite(as.data.frame(ch.6.96.sch.abf[1]), "coloc/coloc_abf_sch_chr6_96_res.csv", row.names = T)
fwrite(as.data.frame(ch.12.57.sch.abf[1]), "coloc/coloc_abf_sch_chr12_57_res.csv", row.names = T)
fwrite(as.data.frame(ch.1.222.sch.abf[1]), "coloc/coloc_abf_sch_chr1_222_res.csv", row.names = T)


#(852 / (852 + 5164))
#(4961 / (4961 + 356180))
#(10647 / (10647 + 350494))
#chr6-96854444	97067028
