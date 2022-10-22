# sig hits in vcf
library(data.table)
library(tidyverse)

setwd("I:/psivakumar/Emer_CH/gwas/uk_analysis/")

sig.hits.vcf <- fread("all.sig.dose.vcf", skip = 42)
all.20.fam <- fread("all_20_noDups.fam")

codeGenotypes <- function(x) {
  y <- strsplit(x, split = ":")[[1]][1]
  if (y == "./.") {
    return(0)
  } else if (y == "0|0") {
    return(0)
  } else if (y == "0/0") {
    return(0)
  } else if (y == "0|1") {
    return(1)
  } else if (y == "0/1") {
    return(1)
  } else if (y == "1|0") {
    return(1)
  } else if (y == "1/0") {
    return(1)
  } else if (y == "1|1") {
    return(2)
  } else if (y == "1/1") {
    return(2)
  } else {
    return(NA)
  }
}

sig.mat <- apply(sig.hits.vcf[, 10:ncol(sig.hits.vcf)], c(1, 2), codeGenotypes)

sig.var.list <- sig.hits.vcf[, "ID", drop = F]
vcf.sample.names <- names(sig.hits.vcf)[-c(1:9)]
vcf.case.samples <- filter(all.20.fam, V6 == 2)$V2 %>% .[. %in% vcf.sample.names]
vcf.control.samples <- filter(all.20.fam, V6 == 1)$V2 %>% .[. %in% vcf.sample.names]


getAltCounts <- function(mat) {
  gene.data <- data.frame(ref_only_count=integer(), ref_only_names=character(), alt_count=integer(), alt_only_names=character(), control_ref_only_count=integer(), control_ref_only_names=character(), control_alt_count=integer(), control_alt_names=character(), case_ref_only_count=integer(), case_ref_only_names=character(), case_alt_count=integer(), case_alt_names=character(), het_count=integer(), het_names=character(), hom_count=integer(), hom_names=character(), control_het_count=integer(), control_het_names=character(), control_hom_count=integer(), control_hom_names=character(), case_het_count=integer(), case_het_names=character(), case_hom_count=integer(), case_hom_names=character())
  df.names <- names(gene.data)
  case.mat <- mat[, vcf.case.samples]
  control.mat <- mat[, vcf.control.samples]
  for (i in 1:nrow(mat)) {
    ref.only.count <- rowSums(mat[i, , drop = FALSE] == 0)
    ref.only.names <- colnames(mat)[which(mat[i, , drop = FALSE] == 0)]
    has.alt.count <- rowSums(mat[i, , drop = FALSE] != 0)
    has.alt.names <- colnames(mat)[which(mat[i, , drop = FALSE] != 0)]
    control.ref.only.count <- rowSums(control.mat[i, , drop = FALSE] == 0)
    control.ref.only.names <- colnames(control.mat)[which(control.mat[i, , drop = FALSE] == 0)]
    control.has.alt.count <- rowSums(control.mat[i, , drop = FALSE] != 0)
    control.has.alt.names <- colnames(control.mat)[which(control.mat[i, , drop = FALSE] != 0)]
    case.ref.only.count <- rowSums(case.mat[i, , drop = FALSE] == 0)
    case.ref.only.names <- colnames(case.mat)[which(case.mat[i, , drop = FALSE] == 0)]
    case.has.alt.count <- rowSums(case.mat[i, , drop = FALSE] != 0)
    case.has.alt.names <- colnames(case.mat)[which(case.mat[i, , drop = FALSE] != 0)]
    het.count <- length(which(mat[i, , drop = FALSE] == 1))
    het.names <- colnames(mat)[which(mat[i, , drop = FALSE] == 1)]
    hom.count <- length(which(mat[i, , drop = FALSE] == 2))
    hom.names <- colnames(mat)[which(mat[i, , drop = FALSE] == 2)]
    control.het.count <- length(which(control.mat[i, , drop = FALSE] == 1))
    control.het.names <- colnames(control.mat)[which(control.mat[i, , drop = FALSE] == 1)]
    control.hom.count <- length(which(control.mat[i, , drop = FALSE] == 2))
    control.hom.names <- colnames(control.mat)[which(control.mat[i, , drop = FALSE] == 2)]
    case.het.count <- length(which(case.mat[i, , drop = FALSE] == 1))
    case.het.names <- colnames(case.mat)[which(case.mat[i, , drop = FALSE] == 1)]
    case.hom.count <- length(which(case.mat[i, , drop = FALSE] == 2))
    case.hom.names <- colnames(case.mat)[which(case.mat[i, , drop = FALSE] == 2)]
    gene.df <- data.frame(ref.only.count, paste(ref.only.names, collapse = ' '), has.alt.count, paste(has.alt.names, collapse = ' '), control.ref.only.count, paste(control.ref.only.names, collapse = ' '), control.has.alt.count, paste(control.has.alt.names, collapse = ' '), case.ref.only.count, paste(case.ref.only.names, collapse = ' '), case.has.alt.count, paste(case.has.alt.names, collapse = ' '), het.count, paste(het.names, collapse = ' '), hom.count, paste(hom.names, collapse = ' '), control.het.count, paste(control.het.names, collapse = ' '), control.hom.count, paste(control.hom.names, collapse = ' '), case.het.count, paste(case.het.names, collapse = ' '), case.hom.count, paste(case.hom.names, collapse = ' '))
    gene.data <- rbind.data.frame(gene.data, gene.df)
  }
  names(gene.data) <- df.names
  return(gene.data)
}

sig.hits.vcf.df <- getAltCounts(sig.mat)
sig.hits.vcf.df <- data.frame(SNP = sig.var.list$ID, sig.hits.vcf.df)
#fwrite(sig.hits.vcf.df, "gwas_sig_hits_samples_summary.csv", quote = T)
sig.hits.vcf.df.counts.only <- dplyr::select(sig.hits.vcf.df, c(1, seq(2, 24, 2)))
