library(data.table)
library(tidyverse)

setwd("I:/psivakumar/Emer_CH/gwas/euro_uk_analysis/SAIGE_results")

saige.results <- list.files(pattern = "*_output.txt") %>%
  lapply(fread)

saige.cleaned <- lapply(saige.results, function(df) filter(df, AC_Allele2 > 3, imputationInfo > 0.3))

saige.df <- bind_rows(saige.cleaned)
#fwrite(saige.df, "imputed_all_results.csv")
saige.sig <- filter(saige.df, p.value < 5E-6)
#fwrite(saige.sig, "imputed_sig_results.csv")
saige.sig.vcf.format <- saige.sig[, 1:2]
#write.table(saige.sig.vcf.format, "imputed_sig_results_vcf_format.txt", col.names = F, quote = F, row.names = F, sep = '\t')
                               
res_log10 <- saige.df %>%
  mutate(log10Praw = -1*log10(p.value)) %>%
  mutate(log10P = ifelse(log10Praw > 40, 40, log10Praw)) %>%
  mutate(log10Plevel = ifelse(p.value < 5E-08, "possible", NA)) %>%
  mutate(log10Plevel = ifelse(p.value < 5E-09, "likely", log10Plevel))

res_filt <- filter(res_log10, log10P > 3.114074) %>%
  arrange(CHR, POS)

plotting <- res_filt %>%
  group_by(CHR) %>%
  summarise(chr_len=max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(res_filt, ., by = "CHR") %>%
  arrange(ordered(CHR), POS) %>%
  mutate(POScum=POS+tot)

axis_df <- plotting %>%
  group_by(CHR) %>%
  summarise(center = (max(POScum) + min(POScum)) / 2)

manhattan <-
  ggplot(plotting, aes(POScum, log10P)) +
  geom_point(aes(colour=as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("light grey", "dark grey"), 22)) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(x = "CHR", y = "-log10P") +
  theme(legend.position = "none") +
  geom_abline(intercept = -log10(5E-08), slope = 0, linetype = 2)

qqplotter <- function(df) {
  df <- df[!is.na(df$p.value), ]
  df <- mutate(df, exp.pvalues = (rank(df$p.value, ties.method="first")+.5)/(length(df$p.value)+1))
  ggplot(df, aes(-log10(exp.pvalues), -log10(p.value))) + 
    geom_point() +
    labs(x = "Expected logP", y = "Observed logP") + 
    geom_abline(slope = 1, intercept = 0) +
    #geom_label_repel(data = arrange(df, p.value)[1:5, ], box.padding = 0.5) +
    theme_classic() 
}

qqplt <- qqplotter(saige.df)

chisq <- qchisq(1-saige.df$p.value, 1)
lambda <- median(chisq) / qchisq(0.5,1)
alt.chisq <- qnorm(saige.df$p.value/2)
alt.lambda <- median(alt.chisq^2, na.rm=T)/qchisq(0.5,df=1)

#get snps around regions of interest
#library(biomaRt)
#mart <- useMart("ENSEMBL_MART_SNP", host = "grch37.ensembl.org")
#dataset <- useDataset("hsapiens_snp_som", mart = mart)
#snps.df <- getBM(mart = dataset, 
#                 filters = c("chr_name", "start", "end"), 
#                 values = list(8, 139050000, 139200000), 
#                 attributes = c("refsnp_id", "chr_name", "chrom_start"))

#create files for locuszoom
locuszoom_data <- lapply(saige.cleaned, function(df) unite(df, "ID", c("CHR", "POS"), sep = ":") %>% dplyr::select("ID", "p.value") %>% arrange(p.value))

chr.order <- c(1, seq(10, 19, 1), 2, seq(20, 22, 1), seq(3, 9, 1))

#for(i in 1:length(locuszoom_data)) {
#  fwrite(locuszoom_data[[i]], paste0("I:/psivakumar/Emer_CH/gwas/euro_uk_analysis/locuszoom/lz_chr", chr.order[i], ".txt", sep = ""), sep = '\t')
#}

# sig hits in vcf

sig.hits.vcf <- fread("all.sig.dose.vcf", skip = 42)

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


#unimputed manhattan
unimp.model <- fread("I:/psivakumar/Emer_CH/gwas/euro_uk_analysis/unimputed_fishers/unimputed_fishers.model")
unimp.model.pos <- unimp.model %>%
  filter(TEST == "ALLELIC", CHR %in% seq(1, 22, 1)) %>% 
  separate(SNP, c(NA, "POS", NA), sep = ":", remove = FALSE) %>%
  mutate(POS = as.numeric(gsub("([0-9]+).*$", "\\1", .$POS))) %>%
  mutate(log10Praw = -1*log10(P)) %>%
  mutate(log10P = ifelse(log10Praw > 40, 40, log10Praw)) %>%
  mutate(cases_AF = apply(.["AFF"], 1, function(x) eval(parse(text = x)))) %>%
  mutate(controls_AF = apply(.["UNAFF"], 1, function(x) eval(parse(text = x))))

unimp.model.plotting <- unimp.model.pos %>%
  group_by(CHR) %>%
  summarise(chr_len=max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(unimp.model.pos, by = "CHR") %>%
  arrange(ordered(CHR), POS) %>%
  mutate(POScum=POS+tot)

unimp.manhattan <-  
  ggplot(unimp.model.plotting, aes(POScum, log10P)) +
  geom_point(aes(colour=as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("light grey", "dark grey"), 22)) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(x = "CHR", y = "-log10P") +
  theme(legend.position = "none") +
  geom_abline(intercept = -log10(5E-08), slope = 0, linetype = 2)

unimp.qq <- rename(unimp.model.plotting, "p.value" = P) %>% qqplotter()
unimp.alt.chisq <- qnorm(unimp.model.plotting$P/2)
unimp.alt.lambda <- median(unimp.alt.chisq^2, na.rm=T)/qchisq(0.5,df=1)

unimp.sig <- filter(unimp.model.plotting, P < 5E-06)
#fwrite(unimp.sig, "I:/psivakumar/Emer_CH/gwas/euro_uk_analysis/unimputed_fishers/unimputed_sig.csv")

unimp.lz.data <- unimp.model.plotting %>% dplyr::select("CHR", "POS", "P") %>% unite("ID", c("CHR", "POS"), sep = ":") %>% arrange(P)
#fwrite(unimp.lz.data, "I:/psivakumar/Emer_CH/gwas/euro_uk_analysis/locuszoom/unimputed_lz_data.txt", sep = '\t')
