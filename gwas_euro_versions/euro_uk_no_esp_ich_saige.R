library(data.table)
library(tidyverse)

setwd("I:/psivakumar/Emer_CH/gwas/euro_uk_no_esp_ich/SAIGE_results")

saige.results <- list.files(pattern = "*_output.txt") %>%
  lapply(fread)

saige.cleaned <- lapply(saige.results, function(df) filter(df, AC_Allele2 > 3, imputationInfo > 0.3))

saige.df <- bind_rows(saige.cleaned)
#fwrite(saige.df, "imputed_all_results.csv")
saige.sig <- filter(saige.df, p.value < 5E-6)
#fwrite(saige.sig, "imputed_sig_results.csv")

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

#unimputed manhattan
unimp.model <- fread("I:/psivakumar/Emer_CH/gwas/euro_uk_no_esp_ich/unimputed_fishers/unimputed_fishers.model")
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

unimp.axis.df <- unimp.model.plotting %>%
  group_by(CHR) %>%
  summarise(center = (max(POScum) + min(POScum)) / 2)

unimp.manhattan <-  
  ggplot(unimp.model.plotting, aes(POScum, log10P)) +
  geom_point(aes(colour=as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_color_manual(values = rep(c("light grey", "dark grey"), 22)) +
  scale_x_continuous(label = unimp.axis.df$CHR, breaks = unimp.axis.df$center) +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(x = "CHR", y = "-log10P") +
  theme(legend.position = "none") +
  geom_abline(intercept = -log10(5E-08), slope = 0, linetype = 2)

unimp.qq <- rename(unimp.model.plotting, "p.value" = P) %>% qqplotter()
unimp.alt.chisq <- qnorm(unimp.model.plotting$P/2)
unimp.alt.lambda <- median(unimp.alt.chisq^2, na.rm=T)/qchisq(0.5,df=1)
#write.table(unimp.alt.lambda, "I:/psivakumar/Emer_CH/gwas/euro_uk_no_esp_ich/unimputed_fishers/lambda_value.txt")

unimp.sig <- filter(unimp.model.plotting, P < 5E-06)
#fwrite(unimp.sig, "I:/psivakumar/Emer_CH/gwas/euro_uk_no_esp_ich/unimputed_fishers/unimputed_sig.csv")

unimp.lz.data <- unimp.model.plotting %>% dplyr::select("CHR", "POS", "P") %>% unite("ID", c("CHR", "POS"), sep = ":") %>% arrange(P)
#fwrite(unimp.lz.data, "I:/psivakumar/Emer_CH/gwas/euro_uk_no_esp_ich/locuszoom/unimputed_lz_data.txt", sep = '\t')
