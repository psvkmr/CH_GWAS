library(data.table)
library(tidyverse)
library(biomaRt)

setwd("I:/psivakumar/Emer_CH/gwas/uk_analysis/SAIGE_results")

# USED DIFFERENT METHOD FOR GETTING RSIDS, JOINED FROM HRC REF FILE

#sig.res <- fread("imputed_sig_res.csv")
#unimp.res <- fread("I:/psivakumar/Emer_CH/gwas/uk_analysis/unimputed_fishers/unimputed_sig.csv")

#mart <- useMart("ENSEMBL_MART_SNP", host = "grch37.ensembl.org")
#dataset <- useDataset("hsapiens_snp", mart = mart)

#variant.identifiers <- list()
#for (i in 1:nrow(sig.res)) {
#  bm.res <- getBM(mart = dataset, 
#                  attributes = c("refsnp_id"),
#                  filters = c("chr_name", "start", "end"), 
#                  values = list(sig.res[[i, 1]], sig.res[[i, 2]], sig.res[[i, 2]])
#  )
#  variant.identifiers[[i]] <- bm.res
#}

# check in multiple hits for which correct
#variant.identifiers.uniq <- sapply(variant.identifiers, function(x) x[[1]][[1]])

#sig.res.ids <- cbind(sig.res, variant.identifiers.uniq)

# change row 362, rsID wrong

#corrected.rsid <- getBM(mart = dataset, 
#                               attributes = c("refsnp_id"),
#                               filters = c("chr_name", "start", "end"), 
#                               values = list(sig.res[[362, 1]], sig.res[[362, 2]], sig.res[[362, 2]])
#                        )[[2, 1]]
#
#sig.res.ids[362, 20] <- corrected.rsid

#fwrite(sig.res.ids, "imputed_sig_res_with_var_ids.csv")
#sig.res.ids <- fread("imputed_sig_res_with_var_ids.csv")

# read summary stats
postgap.sig.sum <- fread("postgap_sig_sumstats.tsv")
high.mod.so <- fread("I:/psivakumar/high_and_moderate_impact_so_terms.txt", header = F)

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
