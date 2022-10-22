library(tidyverse)
library(data.table)

setwd("I:/psivakumar/Emer_CH/burden")

all.log <- read.csv("C:/Users/prasa/Google Drive/all_Log_sequencing.csv")
all.log <- mutate(all.log, LogID = paste0("LI", LogID))


### get wgs metrics for samples

#cases.kit <- filter(all.log, LogID %in% filter(new.fam, V6 == 2)$V2)[1, 42]
#control.kit <- filter(all.log, LogID %in% filter(new.fam, V6 == 1)$V2)[1, 42]

#wgs.metrics <- read.csv("C:/Users/prasa/Google Drive/wgsmetrics_8K.csv")
#study.samples <- read.delim("I:/psivakumar/Emer_CH/burden/all_samples.txt", header = F)
#study.samples.metrics <- filter(wgs.metrics, X %in% study.samples$V1) %>%
  mutate(pheno = ifelse(X %in% emer.cases.df$LOGID, 2, 1))

#samples.high.depth <- filter(study.samples.metrics, MEAN_COVERAGE >= 15)$X
#emer.samples.low.depth <- filter(study.samples.metrics, MEAN_COVERAGE <= 15, pheno == 2)$X

#sample.depth.plot <- ggplot(study.samples.metrics, aes(as.factor(pheno), MEAN_COVERAGE)) + 
#  geom_jitter() + 
#  labs(x = "phenotype status") + 
#  theme_classic()

#write.table(samples.high.depth, "I:/psivakumar/Emer_CH/burden/samples_more_15x_keep.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



### create fam file

emer.fam <- read.table("I:/psivakumar/Emer_CH/burden/emer_ch_1.fam")
emer.cases.df <- read.csv("I:/psivakumar/Emer_CH/emer_id_translation.csv")

new.fam <- mutate(emer.fam, 
                  V6 = ifelse(emer.fam$V2 %in% emer.cases.df$LOGID, 2, 1))

#names(new.fam) <- c("fid", "iid", "mother_id", "father_id", "sex", "pheno")

#write.table(new.fam, "I:/psivakumar/Emer_CH/burden/emer_ch_1.fam", quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')



### ISTATS

istats <- read.delim("I:/psivakumar/Emer_CH/burden/pre_QC_istats.txt")

detectOutliers <- function(x) {
  mean.x <- mean(x)
  sd.x <- sd(x)
  outliers <- c()
  for (i in 1:length(x)){
    if (x[i] > mean.x + 3*sd.x) {
      outliers <- c(outliers, x[i])
    } else if (x[i] < mean.x - 3*sd.x) {
      outliers <- c(outliers, x[i])
    }
  }
  outliers.rows <- match(outliers, x)
  return(outliers.rows)
}

filter.columns <- c(2:7)
istats.outliers <- lapply(istats[, filter.columns], detectOutliers)
outlier.samples <- lapply(istats.outliers, function(x) istats[unlist(x), 1])
all.outlier.samples <- unique(unlist(outlier.samples))
outlier.samples.df <- data.frame(fid = 0, iid = all.outlier.samples)
#write.table(outlier.samples.df, "I:/psivakumar/Emer_CH/burden/istats_outliers_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


## het outliers

hets <- read.table("I:/psivakumar/Emer_CH/burden/het_statistics.het", header = T)
hets$HET_RATE = (hets$N.NM. - hets$O.HOM.) / hets$N.NM.
het.fail = subset(hets, (hets$HET_RATE < mean(hets$HET_RATE)-3*sd(hets$HET_RATE)) | (hets$HET_RATE > mean(hets$HET_RATE)+3*sd(hets$HET_RATE)))
het.fail$HET_DST = (het.fail$HET_RATE - mean(hets$HET_RATE)) / sd(hets$HET_RATE)
het.fail <- het.fail[, 1:2]
#dont remove hets for now
het.fail <- data.frame()
#write.table(het.fail, "I:/psivakumar/Emer_CH/burden/hets_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


## relatives

rel.check <- read.table("I:/psivakumar/Emer_CH/burden/relatives_check.genome", header = T)
duplicate.samples <- filter(rel.check, PI_HAT > 0.95)$IID1
rel.to.remove <- filter(rel.check, PI_HAT > 0.2)[, 1:2] %>% unique.data.frame()
#write.table(rel.to.remove, "I:/psivakumar/Emer_CH/burden/related_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


## sex check

sex.check <- read.table("I:/psivakumar/Emer_CH/burden/sex_check.sexcheck", header = T)
# imputed sexes all look fine so taken no further


## pop stratification

pop.strat <- read.csv("I:/psivakumar/Emer_CH/burden/peddy_emer.het_check.csv")
pop.dist <- table(pop.strat$ancestry.prediction)
european <- filter(pop.strat, ancestry.prediction %in% c("EUR", "UNKNOWN"))$sample_id
#write.table(european, "I:/psivakumar/Emer_CH/burden/european_to_keep.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


### PCA

eigenvecs <- read.table("I:/psivakumar/Emer_CH/burden/pca.eigenvec") 
eigenvals <- read.table("I:/psivakumar/Emer_CH/burden/pca.eigenval")

sum.eig <- sum(eigenvals$V1)

sum.eigs <- lapply(eigenvals$V1, function(x){
  rt<-(x/sum.eig)*100
  rt<-round(rt)
  return(rt)
})

pca.plot <- ggplot(eigenvecs, aes(V3, V4)) +
  geom_point(size=1) +
  xlab(paste0("PC1: ",sum.eigs[[1]],"% variance")) +
  ylab(paste0("PC2: ",sum.eigs[[2]],"% variance"))

pca.filter.columns <- c(3:22)    
pc.outliers <- lapply(eigenvecs[, pca.filter.columns], detectOutliers)
pca.outlier.samples <- lapply(pc.outliers, function(x) eigenvecs[unlist(x), 2])
pca.all.outlier.samples <- unique(unlist(pca.outlier.samples))
pca.outlier.samples.df <- data.frame(fid = 0, iid = pca.all.outlier.samples)
#write.table(pca.outlier.samples.df, "I:/psivakumar/Emer_CH/burden/pca_samples_to_remove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


# create separate case and control samples list

all.samples.final <- read.table("I:/psivakumar/Emer_CH/burden/samples_final_set.txt")
case.samples <- filter(all.samples.final, V1 %in% filter(new.fam, V6 == 2)$V2)
control.samples <- filter(all.samples.final, V1 %in% filter(new.fam, V6 == 1)$V2)
#write.table(case.samples, "I:/psivakumar/Emer_CH/burden/case_samples_final_set.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(control.samples, "I:/psivakumar/Emer_CH/burden/control_samples_final_set.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


## pre vs post QC PCA

pre_qc_pca <- left_join(pop.strat, new.fam, by = c("sample_id" = "V2")) %>%
  ggplot(aes(PC1, PC2, colour = as.factor(V6))) + geom_point()
post_qc_pca <- left_join(pop.strat, new.fam, by = c("sample_id" = "V2")) %>%
  filter(sample_id %in% all.samples.final$V1) %>%
  ggplot(aes(PC1, PC2, colour = as.factor(V6))) + geom_point()

### create list of phenotype specific low quality variants to remove

cases.depth <- read.table("I:/psivakumar/Emer_CH/burden/case_new_depth.ldepth", header = T)
cases.depth.pass <- filter(cases.depth, (SUM_DEPTH / nrow(case.samples)) > 10)
controls.depth <- read.table("I:/psivakumar/Emer_CH/burden/control_new_depth.ldepth", header = T)
controls.depth.pass <- filter(controls.depth, (SUM_DEPTH / nrow(control.samples)) > 10)
depth.pass <- bind_rows(cases.depth.pass[, 1:2], controls.depth.pass[, 1:2]) %>% unique.data.frame()
#write.table(depth.pass, "I:/psivakumar/Emer_CH/burden/read_depth_pass_to_include.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

### create refflat and set files for rvtests burden

library(biomaRt)

gene.mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
gene.dataset <- useDataset(mart = gene.mart, dataset = "hsapiens_gene_ensembl")
genes.pos <- getBM(mart = gene.dataset, 
                   attributes = c("external_gene_name","chromosome_name", "start_position", "end_position"), 
                   filters = "biotype", 
                   values = c("protein_coding"))
tmp.set.file <- dplyr::filter(genes.pos, chromosome_name %in% seq(1, 22, 1))
tmp.set.file$chromosome_name <- gsub("^", "chr", tmp.set.file$chromosome_name)
set.file <- unite(tmp.set.file, "tmp", c("chromosome_name", "start_position"), sep = ":") %>%
  unite("pos", c("tmp", "end_position"), sep = "-")

#write.table(set.file, "I:/psivakumar/Emer_CH/burden/rvtests/hg38_genes.set", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

#refflat <- read.table("J:/NGS_Reference/refFLAT/refFlat_hg38.txt")
#refflat22 <- dplyr::filter(refflat, V3 %in% tmp.set.file$chromosome_name)

#write.table(refflat22, "I:/psivakumar/Emer_CH/burden/rvtests/hg38_refFlat_chr1to22.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

### RESULTS
### Bonferroni correct results

#protein.trunc.assoc <- read.delim("I:/psivakumar/Emer_CH/burden/rvtest_rare_protein_truncating_variants.CMCFisherExact.assoc")
#missense.assoc <- read.delim("I:/psivakumar/Emer_CH/burden/rvtest_rare_missense_variants.CMCFisherExact.assoc")
#syn.assoc <- read.delim("I:/psivakumar/Emer_CH/burden/rvtests/rvtest_rare_syn_variants.CMCFisherExact.assoc")
#nonsyn.assoc <- read.delim("I:/psivakumar/Emer_CH/burden/rvtests/rvtest_rare_nonsyn_variants.CMCFisherExact.assoc")

#protein.trunc.res <- protein.trunc.assoc[!is.na(protein.trunc.assoc$PvalueTwoSide), ] %>% 
#  mutate(padj = p.adjust(PvalueTwoSide, "bonferroni"))
#missense.res <- missense.assoc[!is.na(missense.assoc$PvalueTwoSide), ] %>% 
#  mutate(padj = p.adjust(PvalueTwoSide, "bonferroni"))
#syn.res <- syn.assoc[!is.na(syn.assoc$PvalueTwoSide), ] %>% 
#  mutate(padj = p.adjust(PvalueTwoSide, "bonferroni"))
#nonsyn.res <- nonsyn.assoc[!is.na(nonsyn.assoc$PvalueTwoSide), ] %>%
#  mutate(padj = p.adjust(PvalueTwoSide, "bonferroni"))


## read in results vcfs

# Annotated out inefficient method, code still useful

#getGenotypes <- function(df) {
#  sample.cols <- c(10:ncol(df))
#  gen.mat <- apply(df[, sample.cols], 2, function(x) strsplit(x, split = ":") %>% sapply('[[', 1))
#  row.names(gen.mat) <- paste0(df$X.CHROM, "-", df$POS, "-", df$REF, "-", df$ALT)
#  return(gen.mat)
#}
#
#nonsyn.genotypes <- getGenotypes(nonsyn.vcf)
#syn.genotypes <- getGenotypes(syn.vcf)
#case.nonsyn.genotypes <- nonsyn.genotypes[, unlist(case.samples)]
#control.nonsyn.genotypes <- nonsyn.genotypes[, unlist(control.samples)]
#case.syn.genotypes <- syn.genotypes[, unlist(case.samples)]
#control.syn.genotypes <- syn.genotypes[, unlist(control.samples)]
#
#tableSummarise <- function(m) {
#  var.summary <- apply(m, 1, table)
#  sample.summary <- apply(m, 2, table)
#  sum.list <- list(var.summary, sample.summary)
#  return(sum.list)
#}

#addSummaryCol <- function(df) {
#  sum_list <- tableSummarise(df)
#  df2 <- data.frame(df, as.character(sum_list[[1]]), stringsAsFactors = F)
#  df3 <- rbind.data.frame(df2, as.character(sum_list[[2]]))
#  return(df3)
#}

#nonsyn.genotypes.summary <- addSummaryCol(nonsyn.genotypes)
#syn.genotypes.summary <- addSummaryCol(syn.genotypes)
#case.nonsyn.genotypes.summary <- addSummaryCol(case.nonsyn.genotypes)
#control.nonsyn.genotypes.summary <- addSummaryCol(control.nonsyn.genotypes)
#case.syn.genotypes.summary <- addSummaryCol(case.syn.genotypes)
#control.syn.genotypes.summary <- addSummaryCol(control.syn.genotypes)

#altCounts <- function(df) {
#  alt.counts <- list()
#  hom.counts <- list()
#  het.counts <- list()
#  for (i in 1:nrow(df)) {
#    alt.count <- sum(grepl(paste(c("0/1", "1/0", "1/1"), collapse = "|"), df[i, 1:ncol(df)]))
#    hom.count <- sum(grepl("1/1", df[i, 1:ncol(df)]))
#    het.count <- sum(grepl(paste(c("0/1", "1/0"), collapse = "|"), df[i, 1:ncol(df)]))
#    alt.counts[[i]] <- alt.count
#    hom.counts[[i]] <- hom.count
#    het.counts[[i]] <- het.count
#  }
#  data.frame(df, "alt" = I(alt.counts), "het" = I(het.counts), "hom" = I(hom.counts))
#}

#nonsyn.genotypes.table.counts <- altCounts(nonsyn.genotypes)
#syn.genotypes.table.counts <- altCounts(syn.genotypes)
#case.nonsyn.genotypes.table.counts <- altCounts(case.nonsyn.genotypes)
#control.nonsyn.genotypes.table.counts <- altCounts(control.nonsyn.genotypes)
#case.syn.genotypes.table.counts <- altCounts(case.syn.genotypes)
#control.syn.genotypes.table.counts <- altCounts(control.syn.genotypes)

#addVarGenes <- function(df, setfile) {
#  df <- rownames_to_column(df, "tmp") %>%
#    separate(tmp, c("X.CHROM", "POS", "REF", "ALT"))
#  genes.list <- list()
#  for (i in 1:nrow(df)) {
#    gene <- filter(setfile, chromosome_name == df[i, 1], start_position < df[i, 2], end_position > df[i, 2])$external_gene_name
#    genes.list[[i]] <- gene
#  }
#  #data.frame(df, "genes" = I(genes.list))
#  genes.df <- do.call(rbind, genes.list)
#  df2 <- data.frame(df, genes.df)
#  return(df2)
#}

#nonsyn.genotypes.table.counts.genes <- addVarGenes(nonsyn.genotypes.table.counts, tmp.set.file)
#syn.genotypes.table.counts.genes <- addVarGenes(syn.genotypes.table.counts, tmp.set.file)
#case.nonsyn.genotypes.table.counts.genes <- addVarGenes(case.nonsyn.genotypes.table.counts, tmp.set.file)
#control.nonsyn.genotypes.table.counts.genes <- addVarGenes(control.nonsyn.genotypes.table.counts, tmp.set.file)
#case.syn.genotypes.table.counts.genes <- addVarGenes(case.syn.genotypes.table.counts, tmp.set.file)
#control.syn.genotypes.table.counts.genes <- addVarGenes(control.syn.genotypes.table.counts, tmp.set.file)

#genesVarCounts <- function(df, assoc) {
#  genes.data <- data.frame()
#  for (gene in assoc$Range) {
#    df2 <- filter_at(df, vars(starts_with("X")), any_vars(. == gene))
#    df3 <- data.frame(gene, nvar = nrow(df2), alt_count = sum(unlist(df2$alt)), het_count = sum(unlist(df2$het)), hom_count = sum(unlist(df2$hom)))
#    genes.data <- rbind.data.frame(genes.data, df3)
#  }
#  df4 <- left_join(assoc, genes.data, by = c("Range" = "gene"))
#  return(df4)
#}

#nonsyn.genotypes.full.summary.df <- genesVarCounts(nonsyn.genotypes.table.counts.genes, nonsyn.assoc)
#syn.genotypes.full.summary.df <- genesVarCounts(syn.genotypes.table.counts.genes, syn.assoc)
#case.nonsyn.genotypes.full.summary.df <- genesVarCounts(case.nonsyn.genotypes.table.counts.genes, nonsyn.assoc)
#control.nonsyn.genotypes.full.summary.df <- genesVarCounts(control.nonsyn.genotypes.table.counts.genes, nonsyn.assoc)
#case.syn.genotypes.full.summary.df <- genesVarCounts(case.syn.genotypes.table.counts.genes, syn.assoc)
#control.syn.genotypes.full.summary.df <- genesVarCounts(control.syn.genotypes.table.counts.genes, syn.assoc)

library(bigreadr)


#test.vcf <- fread2("I:/psivakumar/Emer_CH/burden/tmp15_rare_nonsyn.vcf", sep = '\t', nrows = 1000, skip = 3431)
nonsyn.vcf <- fread2("I:/psivakumar/Emer_CH/burden/rare_nonsyn.vcf", sep = '\t', skip = 3432)
syn.vcf <- fread2("I:/psivakumar/Emer_CH/burden/rare_syn.vcf", sep = '\t', skip = 3431)
#full.vcf <- fread2("I:/psivakumar/Emer_CH/burden/tmp12_var_filt.vcf", sep = '\t', skip = 3428)

codeGenotypes <- function(x) {
  y <- strsplit(x, split = ":")[[1]][1]
  if (y == "./.") {
    return(0)
  } else if (y == "0/0") {
    return(0)
  } else if (y == "0/1") {
    return(1)
  } else if (y == "1/0") {
    return(1)
  } else if (y == "1/1") {
    return(2)
  } else {
    return(NA)
  }
}

#test.mat <- apply(test.vcf[, 10:ncol(test.vcf)], c(1, 2), codeGenotypes)
nonsyn.mat <- apply(nonsyn.vcf[, 10:ncol(nonsyn.vcf)], c(1, 2), codeGenotypes)
syn.mat <- apply(syn.vcf[, 10:ncol(syn.vcf)], c(1, 2), codeGenotypes)
#full.mat <- apply(full.vcf[, 10:ncol()])

varAddGenes <- function(vcf, setfile) {
  gene.mat.list <- list(length = nrow(setfile))
  for (i in 1:nrow(setfile)) {
    gene.mat <- which(vcf$`#CHROM` == setfile[i, 2] & vcf$POS > setfile[i, 3] & vcf$POS < setfile[i, 4])
    gene.mat.list[[i]] <- gene.mat
  }
  return(gene.mat.list)
}

#test.gene.list <- varAddGenes(test.vcf, tmp.set.file)
nonsyn.gene.list <- varAddGenes(nonsyn.vcf, tmp.set.file)
syn.gene.list <- varAddGenes(syn.vcf, tmp.set.file)

#getAltCounts <- function(mat, genelist) {
#  gene.data <- data.frame()
#  case.mat <- mat[, unlist(case.samples)]
#  control.mat <- mat[, unlist(control.samples)]
#  for (i in 1:length(genelist)) {
#    nvar <- length(genelist[[i]])
#    alt.count <- sum(rowSums(mat[genelist[[i]], , drop = FALSE] != 0))
#    het.count <- sum(rowSums(mat[genelist[[i]], , drop = FALSE] == 1))
#    hom.count <- sum(rowSums(mat[genelist[[i]], , drop = FALSE] == 2))
#    control.alt <- sum(rowSums(control.mat[genelist[[i]], , drop = FALSE] != 0))
#    control.het <- sum(rowSums(control.mat[genelist[[i]], , drop = FALSE] == 1))
#    control.hom <- sum(rowSums(control.mat[genelist[[i]], , drop = FALSE] == 2))
#    case.alt <- sum(rowSums(case.mat[genelist[[i]], , drop = FALSE] != 0))
#    case.het <- sum(rowSums(case.mat[genelist[[i]], , drop = FALSE] == 1))
#    case.hom <- sum(rowSums(case.mat[genelist[[i]], , drop = FALSE] == 2))
#    gene.df <- c(nvar, alt.count, het.count, hom.count, control.alt, control.het, control.hom, case.alt, case.het, case.hom)
#    gene.data <- rbind(gene.data, gene.df)
#  }
#  names(gene.data) <- c("nvar.count", "alt.count", "het.count", "hom.count", "control.alt.count", "control.het.count", "control.hom.count", "case.alt.count", "case.het.count", "case.hom.count")
#  #mutate(gene.data, ref.only.count = nrow(all.samples.final) - alt.count, control.ref.only.count = nrow(control.samples) - control.alt.count, case.ref.only.count = nrow(case.samples) - case.alt.count)
#  #dplyr::select(gene.data, c(1, 11, 2:4, 12, 5:7, 13, 8:10))
#  return(gene.data)
#}

getAltCounts <- function(mat, genelist) {
  gene.data <- data.frame()
  case.mat <- mat[, unlist(case.samples)]
  control.mat <- mat[, unlist(control.samples)]
  for (i in 1:length(genelist)) {
    nvar <- length(genelist[[i]])
    ref.only.count <- sum(colSums(mat[genelist[[i]], , drop = FALSE] != 0) == 0)
    has.alt.count <- sum(colSums(mat[genelist[[i]], , drop = FALSE] != 0) != 0)
    control.ref.only.count <- sum(colSums(control.mat[genelist[[i]], , drop = FALSE] != 0) == 0)
    control.has.alt.count <- sum(colSums(control.mat[genelist[[i]], , drop = FALSE] != 0) != 0)
    case.ref.only.count <- sum(colSums(case.mat[genelist[[i]], , drop = FALSE] != 0) == 0)
    case.has.alt.count <- sum(colSums(case.mat[genelist[[i]], , drop = FALSE] != 0) != 0)
    het.count <- length(which(mat[genelist[[i]], , drop = FALSE] == 1))
    hom.count <- length(which(mat[genelist[[i]], , drop = FALSE] == 2))
    control.het.count <- length(which(control.mat[genelist[[i]], , drop = FALSE] == 1))
    control.hom.count <- length(which(control.mat[genelist[[i]], , drop = FALSE] == 2))
    case.het.count <- length(which(case.mat[genelist[[i]], , drop = FALSE] == 1))
    case.hom.count <- length(which(case.mat[genelist[[i]], , drop = FALSE] == 2))
    gene.df <- c(nvar, ref.only.count, has.alt.count, control.ref.only.count, control.has.alt.count, case.ref.only.count, case.has.alt.count, het.count, hom.count, control.het.count, control.hom.count, case.het.count, case.hom.count)
    gene.data <- rbind(gene.data, gene.df)
  }
  names(gene.data) <- c("nvar_count", "ref_only_count", "alt_count", "control_ref_only_count", "control_alt_count", "case_ref_only_count", "case_alt_count", "het_count", "hom_count", "control_het_count", "control_hom_count", "case_het_count", "case_hom_count")
  return(gene.data)
}

#test.df <- getAltCounts(test.mat, test.gene.list)
nonsyn.df <- getAltCounts(nonsyn.mat, nonsyn.gene.list)
syn.df <- getAltCounts(syn.mat, syn.gene.list)

#test.df.join <- cbind(tmp.set.file, test.df)
nonsyn.df.join <- cbind(tmp.set.file, nonsyn.df)
syn.df.join <- cbind(tmp.set.file, syn.df)

#test.assoc.join <- left_join(test.df.join, nonsyn.assoc, by = c("external_gene_name" = "Range"))
#nonsyn.assoc.join <- left_join(nonsyn.df.join, nonsyn.assoc, by = c("external_gene_name" = "Range"))
#syn.assoc.join <- left_join(syn.df.join, syn.assoc, by = c("external_gene_name" = "Range"))

runFishers <- function(df_join) {
  fishers.res <- data.frame()
  for (i in 1:nrow(df_join)) {
    test.res <- fisher.test(matrix(c(df_join[i, 8], df_join[i, 9], df_join[i, 10], df_join[i, 11]), nrow = 2))
    res.vec <- c(test.res[[1]], test.res[[3]], test.res[[2]][1], test.res[[2]][2])
    fishers.res <- rbind(fishers.res, res.vec)
  }
  names(fishers.res) <- c("pvalue", "or", "lower_ci", "upper_ci")
  fishers.res <- cbind(df_join, fishers.res)
  return(fishers.res)
}

#test.fishers.res <- runFishers(test.df.join)
nonsyn.fishers.res <- runFishers(nonsyn.df.join)
syn.fishers.res <- runFishers(syn.df.join)

# genes with duplicate location regions
#problem.genes <- filter(nonsyn.assoc.join, paste0(chromosome_name, ":", start_position, "-", end_position) != RANGE)$external_gene_name %>% unique()

nonsyn.cleaned <- nonsyn.fishers.res %>% #[!nonsyn.fishers.res$external_gene_name %in% problem.genes, ] %>%
  filter(alt_count != 0) %>%
  #select(c(1:11, 18:19, 23:26)) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>%
  arrange(pvalue)

#fwrite(nonsyn.cleaned, "rare_nonsynonymous_burden_results.csv")

syn.cleaned <- syn.fishers.res %>% #[!syn.fishers.res$external_gene_name %in% problem.genes, ] %>%
  filter(alt_count != 0) %>%
  #select(c(1:11, 18:19, 23:26)) %>%
  mutate(padj = p.adjust(pvalue, method = "bonferroni")) %>%
  arrange(pvalue)

#fwrite(syn.cleaned, "rare_synonymous_burden_results.csv")

## qqplots

library(ggrepel)

qqplotter <- function(df) {
  df <- df[!is.na(df$pvalue), ]
  df <- mutate(df, exp.pvalues = (rank(df$pvalue, ties.method="first")+.5)/(length(df$pvalue)+1))
  ggplot(df, aes(-log10(exp.pvalues), -log10(pvalue), label = external_gene_name)) + 
    geom_point() +
    labs(x = "Expected logP", y = "Observed logP") + 
    geom_abline(slope = 1, intercept = 0) +
    geom_label_repel(data = arrange(df, pvalue)[1:5, ], box.padding = 0.5) +
    theme_classic() 
}


nonsyn.qq <- qqplotter(nonsyn.cleaned)
syn.qq <- qqplotter(syn.cleaned)

nonsyn.chisq <- qchisq(1-nonsyn.cleaned$pvalue, 1)
nonsyn.lambda <- median(nonsyn.chisq) / qchisq(0.5,1)

syn.chisq <- qchisq(1-syn.cleaned$pvalue, 1)
syn.lambda <- median(syn.chisq) / qchisq(0.5,1)

# go terms

library(gProfileR)

goTerms <- function(df_cleaned) {
  go.terms <- gprofiler(query = df_cleaned[1:2000, 1], 
                        organism = "hsapiens", 
                        underrep = F, 
                        min_set_size = 5, 
                        custom_bg = df_cleaned$external_gene_name)
}

nonsyn.go <- goTerms(nonsyn.cleaned)
syn.go <- goTerms(syn.cleaned)

#write.table(nonsyn.cleaned, "I:/psivakumar/Emer_CH/burden/final_results/nonsynonymous_rare_variant_burden_results.tab", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(syn.cleaned, "I:/psivakumar/Emer_CH/burden/final_results/synonymous_rare_variants_burden_results.tab", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(nonsyn.go, "I:/psivakumar/Emer_CH/burden/final_results/nonsynonymous_rare_variant_burden_go_terms.tab", row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(syn.go, "I:/psivakumar/Emer_CH/burden/final_results/synonymous_rare_variants_burden_go_terms.tab", row.names = FALSE, quote = FALSE, sep = '\t')
