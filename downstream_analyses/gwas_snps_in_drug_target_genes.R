library(data.table)
library(tidyverse)
library(biomaRt)

drug.gene.targets <- fread("I:/psivakumar/drug_gene_list.txt")
uk.sw.sig.hits <- fread("I:/psivakumar/Emer_CH/gwas/uk_swedish_alt/SAIGE_results/imputed_sig_res.csv")

gene.mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
gene.dataset <- useDataset(mart = gene.mart, dataset = "hsapiens_gene_ensembl")
genes.pos <- getBM(mart = gene.dataset, 
                   attributes = c("external_gene_name","chromosome_name", "start_position", "end_position", "gene_biotype"), 
                   filters = c("biotype"), 
                   values = list("protein_coding"))

dt.uk.sw.sig.hits <- setDT(uk.sw.sig.hits)
dt.uk.sw.sig.hits$END <- dt.uk.sw.sig.hits$POS
dt.genes.pos <- filter(genes.pos, chromosome_name %in% seq(1, 22, 1)) %>% 
  mutate(chromosome_name = as.integer(chromosome_name)) %>% 
  setDT()
setkey(dt.uk.sw.sig.hits, CHR, POS, END)
uk.sw.sig.hits.genes <- foverlaps(dt.genes.pos, dt.uk.sw.sig.hits, by.x = c("chromosome_name", "start_position", "end_position"), type = "any", nomatch = 0L)

uk.sw.sig.hits.genes.drug.targets <- filter(drug.gene.targets, hgnc_names %in% uk.sw.sig.hits.genes$external_gene_name)


dt.uk.sw.sig.hits.1mb <- mutate(dt.uk.sw.sig.hits, POS - 1000000, END + 1000000) %>% setDT()
setkey(dt.uk.sw.sig.hits.1mb, CHR, POS, END)
uk.sw.sig.hits.genes.1mb <- foverlaps(dt.genes.pos, dt.uk.sw.sig.hits.1mb, by.x = c("chromosome_name", "start_position", "end_position"), type = "any", nomatch = 0L)

uk.sw.sig.hits.genes.1mb.drug.targets <- filter(drug.gene.targets, hgnc_names %in% uk.sw.sig.hits.genes.1mb$external_gene_name)
