library(tidyverse)
library(data.table)
library(gProfileR)
library(biomaRt)

setwd("/array/psivakumar/Emer_CH/gwas/uk_analysis/")

ch.full.res <- fread("SAIGE_results/imputed_all_results_uk_only.csv")
ch.sig.regions <- data.frame("CHR" = c(2,2,6,12), "start" = c(112656908,200360782,96854444,57497005), "end" = c(112785237,200513620,97067028,57614217))
ch.genes <- fread("FUMA/magma.genes.out")
ch.sig.genes <- filter(ch.genes, P < 5E-4)

gene.mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org")
gene.dataset <- useDataset(mart = gene.mart, dataset = "hsapiens_gene_ensembl")
genes.pos <- getBM(mart = gene.dataset,
                   attributes = c("external_gene_name","chromosome_name", "start_position", "end_position", "gene_biotype"),
                   filters = c("biotype"),
                   values = list("protein_coding"))
set.file <- filter(genes.pos, chromosome_name %in% seq(1, 22, 1))

dt.ch.sig.regions <- setDT(ch.sig.regions)
dt.set.file <- setDT(mutate_at(set.file, vars(c(chromosome_name, start_position, end_position)), as.integer))
setkey(dt.ch.sig.regions, CHR, start, end)
gene.overlaps.sig <- foverlaps(dt.set.file, dt.ch.sig.regions, by.x=c("chromosome_name", "start_position", "end_position"), type = "any", nomatch=0L)

dt.ch.sig.reg.1mb <- setDT(mutate(ch.sig.regions, start = start - 1000000, end = end + 1000000))
setkey(dt.ch.sig.reg.1mb, CHR, start, end)
gene.overlaps.sig.1mb <- foverlaps(dt.set.file, dt.ch.sig.reg.1mb, by.x=c("chromosome_name", "start_position", "end_position"), type = "any", nomatch=0L)

dt.ch.sig.genes <- setDT(ch.sig.genes)
setkey(dt.ch.sig.genes, CHR, START, STOP)
gene.overlaps.sig.genes <- foverlaps(dt.set.file, dt.ch.sig.genes, by.x=c("chromosome_name", "start_position", "end_position"), type = "any", nomatch=0L)

goTerms <- function(df) {
  gost(query = df$external_gene_name,
          organism = "hsapiens",
          custom_bg = set.file$external_gene_name,
          correction_method = "fdr")
}

go.res.reg.1mb <- goTerms(gene.overlaps.sig.1mb)
go.res.reg.1mb.df <- filter(go.res.reg.1mb$result, term_size > 5)
go.plot <- gostplot(go.res.reg.1mb, interactive = F)
fwrite(go.res.reg.1mb.df, "gprofiler2/go_res_1mb_sig_loci.csv")
publish_gostplot(go.plot, filename = "gprofiler2/go_manhattan_1mb_sig_loci.png")

go.res.sig.gene <- goTerms(gene.overlaps.sig.genes)
go.res.sig.gene.df <- filter(go.res.sig.gene$result, term_size > 5)
#go.sig.gene.plot <- gostplot(go.res.sig.gene, "gprofiler2/go_res_sig_gene.csv")
fwrite(go.res.sig.gene.df, "gprofiler2/go_res_sig_gene.csv")
