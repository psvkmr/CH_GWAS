library(tidyverse)
library(data.table)

setwd("I:/psivakumar/Emer_CH/gwas/greek_only_gwas/")

fam_cols <- c(rep("character", 2), rep("integer", 4))
gsa2019.fam <- fread("GSA2019_432_025_CC.fam", colClasses = fam_cols)

greek.cases <- fread("I:/psivakumar/Emer_CH/gwas/greek_cases.csv", header = F)
greek.controls <- fread("I:/psivakumar/Emer_CH/gwas/greek_controls.csv", header = F)

all.ids <- c(greek.cases$V1, greek.controls$V1)

greek.in.fam <- all.ids %in% gsa2019.fam$V2

# remove monomorphic and ambiguous strand

gsa2019.frq <- fread("alleleFreqs_gsa2019.frq")
frq.list <- mget(ls(pattern = ".frq"))

getMonomorphics <- function(df) {
  filter(df, MAF < 0.0000001)$SNP
}

getAmbiguous <- function(df) {
  filter(df, 
         (A1 == "A" & A2 == "T") | 
           (A1 == "T" & A2 == "A") |
           (A1 == "C" & A2 == "G") | 
           (A1 == "G" & A2 == "C"))$SNP
}

monomorphic.list <- lapply(frq.list, getMonomorphics)
ambiguous.list <- lapply(frq.list, getAmbiguous)

mono.or.ambig <- c(unlist(monomorphic.list), unlist(ambiguous.list)) %>% unique()

#write.table(mono.or.ambig, "monomorphic_or_ambig_strand_variants.txt", col.names = F, row.names = F, quote = F)

# check map files
gsa2019.premerge.map <- fread("gsa2019_3_biallelic.map", header = F)

getOtherChr <- function(df) {
  snps <- df[df$V1 %in% c(0, 24, 25, 26), 2]
}

addUniquePos <- function(df){
  mutate(df, V5 = paste0(V1, ":", V4, "-", V4))
}

getCommonPos <- function(lst) {
  compos <- Reduce(intersect, lapply(lst, `[[`, 5))
  snplst <- lapply(lst, function(x) x[which(x$V5 %in% compos), 2])
  comsnp <- Reduce(intersect, snplst)
}

premerge.map.list <- mget(ls(pattern = ".map"))

premerge.common.snps <- getCommonPos(lapply(premerge.map.list, addUniquePos)) %>% .[!. %in% mono.or.ambig]

#write.table(premerge.common.snps, "commonSnpsToExtract.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


# read in all samples fam file

all.fam <- fread("gsa2019_4_commonSnps.fam")
all.fam$V6 <- ifelse(all.fam$V1 %in% greek.cases$V1, 2, 
                     ifelse(all.fam$V1 %in% greek.controls$V1, 1, -9))
#fwrite(all.fam, "gsa2019_4_commonSnps.fam", col.names = F, sep = '\t')


# look at snps to be flipped

flipscan <- fread("postMerge.flipscan", colClasses = c("integer", "character", "integer", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "character"), header = T, fill = T)
non.empty.scan <- filter(flipscan, R_NEG > 0)

#write.table(non.empty.scan$SNP, "strandIssueSnpsToRemove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# check imiss stats

cases.imiss <- fread("casesMissingStats.imiss", header = T)
controls.imiss <- fread("controlsMissingStats.imiss", header = T)

#check sex check

cases.sex <- fread("caseSex.sexcheck", header = T)
controls.sex <- fread("controlSex.sexcheck", header = T)
case.problem.sex <- filter(cases.sex, F < 0.7 & F > 0.3)
control.problem.sex <- filter(controls.sex, F < 0.7 & F > 0.3)
#fwrite(case.problem.sex[, 1:2], "caseFailSexCheck.txt", col.names = F, sep = '\t')
#fwrite(control.problem.sex[, 1:2], "controlFailSexCheck.txt", col.names = F, sep = '\t')

# check re-merged allele freqs
control.recalc.freq <- fread("controlRecalcFreq.frq", header = T)
case.recalc.freq <- fread("casesRecalcFreq.frq", header = T)
all.recalc.freq <- left_join(control.recalc.freq, case.recalc.freq, by = "SNP") %>%
  mutate(DIFF = abs(MAF.y - MAF.x), logDIFF = log10(MAF.y) - log10(MAF.x)) %>%
  arrange(desc(DIFF))
recalc.freq.plot <- ggplot(all.recalc.freq, aes(MAF.x, MAF.y)) + geom_point()


# check het
all.hets <- fread("allHets.het", header = T)
all.hets$HET_RATE = (all.hets$`N(NM)` - all.hets$`O(HOM)`) / all.hets$`N(NM)`
#het.fail = subset(all.hets, (all.hets$HET_RATE < mean(all.hets$HET_RATE)-3*sd(all.hets$HET_RATE)) | (all.hets$HET_RATE > mean(all.hets$HET_RATE)+3*sd(all.hets$HET_RATE)))
#het.fail$HET_DST = (het.fail$HET_RATE - mean(all.hets$HET_RATE)) / sd(all.hets$HET_RATE)
#het.fail <- het.fail[, 1:2]
het.fail <- filter(all.hets, HET_RATE > 0.5 | HET_RATE < 0.1)
#write.table(het.fail, "hetOutliersToRemove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


# check fam file for duplicate ids
all.20.fam <- fread("all_20_noDups.fam")


# peddy pop strat

pop.strat <- fread("peddy_gwas.het_check.csv")
pop.dist <- table(pop.strat$`ancestry-prediction`)
european <- filter(pop.strat, `ancestry-prediction` == "EUR" | `ancestry-prediction` == "UNKNOWN")$sample_id
european.to.keep <- filter(all.20.fam, V2 %in% european)[, 1:2]
#fwrite(european.to.keep, "european_to_keep.txt", col.names = F, sep = '\t')


# pca

eigenvecs <- fread("pca.eigenvec") 
eigenvals <- fread("pca.eigenval")

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

pca.filter.columns <- c(3:22)    
pc.outliers <- lapply(eigenvecs[, ..pca.filter.columns], detectOutliers)
pca.outlier.samples <- lapply(pc.outliers, function(x) eigenvecs[unlist(x), 2])
pca.all.outlier.samples <- unique(unlist(pca.outlier.samples))
#pca.outlier.samples.df <- data.frame(fid = pca.all.outlier.samples, iid = pca.all.outlier.samples)
# not excluding outliers for now
#pca.outlier.samples.df <- filter(all.20.fam, V2 %in% pca.all.outlier.samples)[, 1:2]
#not removing outliers
pca.outlier.samples.df <- data.frame()
#write.table(pca.outlier.samples.df, "pcaSamplesToRemove.txt", col.names = F, sep = '\t')

# check batch effects

info.fam <- all.20.fam %>%
  separate(V2, c("old_WT_FID", "old_FID", "old_IID"), sep = "_", remove = F, fill = "left") %>%
  #mutate(pheno = ifelse(old_IID %in% euro.cases, 2, 1)) %>%
  mutate(batch = ifelse(old_IID %in% euro.samples.list[[1]]$V1, yes = "spanish", no = 
                          ifelse(old_IID %in% euro.samples.list[[2]]$V1, yes = "belgian", no = 
                                   ifelse(old_IID %in% euro.samples.list[[3]]$V1, yes = "danish", no = 
                                            ifelse(old_IID %in% euro.samples.list[[4]]$V1, yes = "greek", no = 
                                                     ifelse(old_IID %in% euro.samples.list[[5]]$V1, yes = "greek_control", no = 
                                                              ifelse(old_IID %in% ich.controls.fam$V2, yes = "ich_control", no = 
                                                                       ifelse(old_IID %in% gsa2016.fam$V2, yes = "gsa2016", no = 
                                                                                ifelse(old_IID %in% gsa2018.fam$V2, yes = "gsa2018", no = 
                                                                                         ifelse(old_IID %in% reseq.fam$V2, yes = "reseq", no = 
                                                                                                  ifelse(old_IID %in% c1958.fam$V2, yes = "1958_control", no = 
                                                                                                           ifelse(old_IID %in% nbs.fam$V2, yes = "nbs_control", no = "other"))))))))))))


# pre and post pop  strat/pca qc

pre.qc.pca <- left_join(eigenvecs, info.fam, by = "V2") %>%
  ggplot(aes(V3.x, V4.x, colour = as.factor(batch), shape = as.factor(V6.y))) +
  geom_point(size=2, alpha = 0.7) +
  scale_color_brewer(palette = "Paired") +
  theme_classic()

post.qc.pca <- left_join(eigenvecs, info.fam, by = "V2") %>%
  .[!.$V2 %in% pca.all.outlier.samples, ] %>%
  ggplot(aes(V3.x, V4.x, colour = as.factor(batch))) + 
  geom_point(size = 2, alpha = 0.4) +
  theme_classic()
#ggsave("pre_qc_PC1_PC2.png", pre.qc.pca, device = "png")
#ggsave("post_qc_PC1_PC2.png", post.qc.pca, device = "png")

# SAIGE
# pheno 

final.samples <- fread("ID_order_chrX.txt", header = F)
saige.pheno <- left_join(all.20.fam, eigenvecs, by = "V2") %>%
  filter(V2 %in% final.samples$V1) %>%
  mutate(new_pheno = ifelse(V6.x == 2, 1, ifelse(V6.x == 1, 0, NA))) %>%
  select(FID = V1.x, IID = V2, PATID = V3.x, MATID = V4.x, SEX = V5.x, Phenotype = new_pheno, PC1 = V3.y, PC2 = V4.y, PC3 = V5.y, PC4 = V6.y, PC5 = V7, PC6 = V8, PC7 = V9, PC8 = V10, PC9 = V11, PC10 = V12, PC11 = V13, PC12 = V14, PC13 = V15, PC14 = V16, PC15 = V17, PC16 = V18, PC17 = V19, PC18 = V20, PC19 = V21, PC20 = V22)
#fwrite(saige.pheno, "saige.pheno", sep = '\t')
