library(data.table)
library(tidyverse)
library(liftOver)
library(coloc)

setwd("I:/psivakumar/Emer_CH/gwas/uk_analysis")

test.egene <- fread("J:/NGS_Reference/GTEX/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt")
test.egene.format <- separate(test.egene, variant_id, c("chr", "pos", "ref", "alt", NA), sep = "_", remove = F) %>%
  unite("chr_pos", c("chr", "pos"), sep = "_")
test.signif <- fread("J:/NGS_Reference/GTEX/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt")

egenes.list <- list.files(path = "J:/NGS_Reference/GTEX/GTEx_Analysis_v8_eQTL/", pattern = "egenes.txt")
egenes.list.names <- sapply(egenes.list, strsplit, "\\.") %>% sapply(`[[`, 1)
egenes.df <- data.frame()
for (i in 1:length(egenes.list)){
      temp.data <- fread(paste0("J:/NGS_Reference/GTEX/GTEx_Analysis_v8_eQTL/", egenes.list[i], sep = ""))
      temp.data$tissue <- egenes.list.names[i]
      egenes.df <- rbind(egenes.df, temp.data) 
}
egenes.df.format <- separate(egenes.df, variant_id, c("chr", "pos", "ref", "alt", NA), sep = "_", remove = F) %>%
  unite("chr_pos", c("chr", "pos"), sep = "_")

gwas.sig.res <- fread("SAIGE_results/imputed_sig_res_with_var_ids.csv")
lead.snps <- fread("FUMA/leadSNPs.txt")
lead.snps.in.gwas <- filter(gwas.sig.res, POS %in% lead.snps$pos)

hg19.bed <- gwas.sig.res %>%
  filter(POS %in% lead.snps$pos) %>%
  dplyr::select(CHR, POS ,SNPID) %>%
  mutate(end = POS + 1)
hg19.bed$CHR <- gsub("^", "chr", hg19.bed$CHR)

#fwrite(hg19.bed, "sig_res_to_liftover.txt", col.names = F, sep = '\t')

hg38.chain <- import.chain("J:/NGS_Reference/hg19ToHg38.over.chain")

hg19.bed.grange <- makeGRangesFromDataFrame(hg19.bed, ignore.strand = T, seqnames.field = "CHR", start.field = "POS", end.field = "end", keep.extra.columns = T)

hg38.lifted <- liftOver(hg19.bed.grange, hg38.chain)
hg38.lifted.df <- data.frame("chr" = hg38.lifted@unlistData@seqnames, "range" = hg38.lifted@unlistData@ranges, "snpid" = hg38.lifted@unlistData@elementMetadata@listData) %>%
  unite("chr_pos", c(chr, range.start), sep = "_", remove = F)

gwas.in.eqtl <- filter(test.egene.format, gene_chr == as.character(hg38.lifted.df[[1, 2]]), (gene_start - 200000) < hg38.lifted.df[[1, 3]] & (gene_end + 200000) > hg38.lifted.df[[1, 3]]) %>% mutate()
coloc::coloc.abf(dataset1 = list(beta = lead.snps.in.gwas[[1, 10]], varbeta = lead.snps.in.gwas[[1, 11]], type = "cc", s = 0.142), 
                 dataset2 = list(beta = gwas.in.eqtl[[1, 25]], varbeta = gwas.in.eqtl[[1, 26]], type = "quant", MAF = gwas.in.eqtl[[1, 22]]))

gwas.in.eqtl <- rename_all(egenes.df.format, function(x) paste("eqtl", x, sep = "_")) %>%
  inner_join(rename_all(hg38.lifted.df, function(x) paste("gwas", x, sep = "_")), by = c("eqtl_chr_pos" = "gwas_chr_pos")) %>%
  filter(eqtl_qval < 0.05)
#fwrite(gwas.in.eqtl, "gwas_hits_in_eqtl_data.csv")
