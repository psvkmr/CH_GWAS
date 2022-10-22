library(data.table)
library(tidyverse)
library(biomaRt)

setwd("I:/psivakumar/Emer_CH/gwas/uk_analysis/")

postgap.sig.sum <- fread("postgap/lead_SNPs_postgap.tsv")
high.mod.so <- fread("I:/psivakumar/high_and_moderate_impact_so_terms.txt", header = F)
sig.res <- fread("SAIGE_results/imputed_sig_results_e5_uk_rsIDs.tsv")
unimp.res <- fread("unimputed_fishers/unimputed_fishers.model")

# for interpreting sum stats
# https://github.com/Ensembl/postgap/wiki/How-do-I-use-POSTGAP-output%3F

gene.hit.counts <- as.data.frame(table(postgap.sig.sum$gene_symbol)) %>% arrange(desc(Freq))
var.hit.counts <- as.data.frame(table(postgap.sig.sum$gwas_snp)) %>% arrange(desc(Freq))
cluster.hit.counts <- as.data.frame(table(postgap.sig.sum$cluster_id)) %>% arrange(desc(Freq))

postgap.sig.low.eu.maf <- filter(postgap.sig.sum, gnomad_nfe < 0.05 | eur < 0.05) %>% distinct_at(.vars = c("ld_snp_rsID", "gene_symbol", "score"), .keep_all = T)
postgap.sig.first.rank.genes <- filter(postgap.sig.sum, rank == 1) %>% distinct_at(.vars = "ld_snp_rsID", .keep_all = T) %>% arrange(desc(score))
postgap.sig.high.mod.so <- postgap.sig.sum[grep(paste(high.mod.so$V1, collapse = "|"), postgap.sig.sum$vep_terms), ] %>% distinct_at(.vars = c("ld_snp_rsID", "gene_symbol", "score"), .keep_all = T)
postgap.sig.brain.gtex.sig <- filter_at(postgap.sig.sum, vars(contains("Brain")), any_vars(. > 0.95)) %>% distinct_at(.vars = c("ld_snp_rsID", "gene_symbol", "score"), .keep_all = T)
mertk.brain.hits <- filter(postgap.sig.brain.gtex.sig, gene_symbol == "MERTK")

fuma.genes <- fread("FUMA/genes.txt")
#fuma.leadsnps <- fread("FUMA/leadSNPs.txt")
fuma.snps <- fread("FUMA/snps.txt")
fuma.eqtl <- fread("FUMA/eqtl.txt")
ukbb.mig.sig <- fread("ldsc/ukbb_migraine_sig_sumstats.txt")
phewas <- fread("I:/psivakumar/Emer_CH/gwas/phewas-catalog.csv")
fuma.ld <- fread("FUMA/ld.txt")
lead.snps <- c("rs4519530", "rs6435024", "rs9386670")
fuma.leadsnps <- filter(fuma.snps, rsID %in% lead.snps)

sig.res.lead <- filter(sig.res, ID %in% lead.snps)
#fwrite(sig.res.lead, "imputed_gwas_sig_res_lead_SNPs.tsv", sep = '\t')
fuma.sig.genes <- c("MERTK", "FHL5", "UFL1", "STAT6", "GABRB1", "TMEM78B")
nearby.genes <- filter(fuma.genes, ((start - 1000000) < fuma.leadsnps[[1, 4]]) & ((end + 1000000) > fuma.leadsnps[[1, 4]]) | ((start - 1000000) < fuma.leadsnps[[2, 4]]) & ((end + 1000000) > fuma.leadsnps[[2, 4]]) | ((start - 1000000) < fuma.leadsnps[[3, 4]]) & ((end + 1000000) > fuma.leadsnps[[3, 4]]))
missense.vars <- postgap.sig.high.mod.so$ld_snp_rsID
missense.vars.fuma.anno <- filter(fuma.snps, rsID %in% missense.vars)
postgap.lead.only <- filter(postgap.sig.sum, ld_snp_rsID %in% fuma.leadsnps$rsID)
postgap.lead.sig.genes <- filter(postgap.lead.only, gene_symbol %in% fuma.sig.genes)
postgap.coronary.artery.sig <- filter(postgap.sig.sum, GTEx_Artery_Coronary > 0.95)
postgap.coronary.artery.sig.lead.snps <- filter(postgap.sig.sum, GTEx_Artery_Coronary > 0.95, ld_snp_rsID %in% fuma.leadsnps$rsID)
postgap.lead.brain.gtex.sig <- filter_at(postgap.sig.sum, vars(contains("Brain")), any_vars(. > 0.995)) %>% filter(ld_snp_rsID %in% fuma.leadsnps$rsID) 
postgap.lead.brain.gtex.sig <- postgap.lead.brain.gtex.sig[, c(1:3, 20:47, grep("Brain", names(postgap.sig.sum)), 92:100)]
fuma.lead.coronary.artery.sig <- filter(fuma.eqtl, tissue == "Artery_Coronary", uniqID %in% fuma.leadsnps$uniqID)
phewas.lead <- phewas[grep(paste(fuma.leadsnps$rsID, collapse = "|"), phewas$snp), ]
phewas.sig <- filter(phewas, snp %in% postgap.sig.sum$ld_snp_rsID)
sig.res.or <- mutate(sig.res, or = exp(BETA))
sig.res.or.lead <- filter(sig.res.or, ID %in% lead.snps)

lead snp 200512641
satb2 200335989
