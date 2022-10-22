library(data.table)
library(tidyverse)

setwd("I:/psivakumar/Emer_CH/gwas/uk_analysis/")

ch.sumstats <- fread("ldsc/imputed_all_results_uk_rsIDS.tsv")
mig.sumstats <- fread("ldsc/ukbb_ukbb_migraine_sumstats_full_sumstats.txt")
ldsnps <- fread("FUMA/ld.txt")

ch.chr6 <- filter(ch.sumstats, CHR == 6, POS > 95000000, POS < 98000000)
mig.chr6 <- filter(mig.sumstats, CHR == 6, POS > 95000000, POS < 98000000)

ch.mig.chr6 <- inner_join(ch.chr6, mig.chr6, by = "ID")

ch.mig.chr6.overlap.plt <- ch.mig.chr6 %>% 
  select(POS.x, p.value, pval) %>% 
  gather(whichP, P, -POS.x) %>% 
  mutate(whichP = ifelse(whichP == "p.value", "CH", "MIG")) %>%
  ggplot(aes(x = POS.x, y = -log10(P), colour = whichP)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth()

ch.chr12 <- filter(ch.sumstats, CHR == 12)
mig.chr12 <- filter(mig.sumstats, CHR == 12)
ch.mig.chr12 <- inner_join(ch.chr12, mig.chr12, by = "ID")
ch.mig.chr12.overlap.plt <- ch.mig.chr12 %>% 
  select(POS.x, p.value, pval) %>% 
  gather(whichP, P, -POS.x) %>% 
  mutate(whichP = ifelse(whichP == "p.value", "CH", "MIG")) %>%
  ggplot(aes(x = POS.x, y = -log10(P), colour = whichP)) + 
  geom_point(alpha = 0.1) + 
  geom_smooth()

ch.mig.2.112.sum <- fread("coloc/ch_mig_chr2_112_sumstats.csv")
ch.mig.2.200.sum <- fread("coloc/ch_mig_chr2_200_sumstats.csv")
ch.mig.6.96.sum <- fread("coloc/ch_mig_chr6_96_sumstats.csv")
ch.mig.12.57.sum <- fread("coloc/ch_mig_chr12_57_sumstats.csv")

pltCommonSums <- function(df) {
  df %>%
  select(POS.x, p.value, pval) %>% 
  gather(whichP, P, -POS.x) %>% 
  mutate(whichP = ifelse(whichP == "p.value", "CH", "MIG")) %>%
  ggplot(aes(x = POS.x, y = -log10(P), colour = whichP)) + 
  geom_point() + 
  geom_smooth() +
  theme_classic()
}

ch.mig.2.112.plt <- pltCommonSums(ch.mig.2.112.sum)
ch.mig.2.200.plt <- pltCommonSums(ch.mig.2.200.sum)
ch.mig.6.96.plt <- pltCommonSums(ch.mig.6.96.sum)
ch.mig.12.57.plt <- pltCommonSums(ch.mig.12.57.sum)
