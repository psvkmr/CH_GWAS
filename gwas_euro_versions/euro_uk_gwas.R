library(tidyverse)
library(data.table)

# uk samples

reseq.samplesheet <- fread("I:/psivakumar/Emer_CH/Plate 3_samplesheet.csv", skip = 16)
reseq.emer.samples <- fread("I:/psivakumar/Emer_CH/Emer_resampled.csv")
uk.reseq.samples <- filter(reseq.samplesheet, Sample_ID %in% reseq.emer.samples$ID_2) %>%
  unite("id", c(SentrixBarcode_A, SentrixPosition_A), sep = "_") %>%
  .[['id']]

uk.samples <- read.csv("I:/psivakumar/Emer_CH/gwas/samplesheet_UK_only_final.csv", stringsAsFactors = FALSE)
fam_cols <- c(rep("character", 2), rep("integer", 4))
gsa2016.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2016_142_025_CC/PLINK_210318_0250/GSA2016_142_025_CC.fam", colClasses = fam_cols)
gsa2018.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2018_310_025-CC/PLINK_150319_1054/GSA2018_310_025-CC.fam", colClasses = fam_cols)
gsa2019.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2019_432_025_CC/PLINK_250719_0647/GSA2019_432_025_CC.fam", colClasses = fam_cols)
c1958.fam <- read.table("I:/psivakumar/Emer_CH/gwas/1958_controls/1958controls_14.09.2018.fam", colClasses = fam_cols)
nbs.fam <- read.table("I:/psivakumar/Emer_CH/gwas/NBS_controls/NBScontrols_14.09.2018.fam", colClasses = fam_cols)
reseq.fam <- read.table("I:/psivakumar/Emer_CH/gwas/Plate_3_Emer/Plate_3_Emer.fam", colClasses = fam_cols)

fam.list <- list(gsa2016.fam, gsa2018.fam, c1958.fam, nbs.fam, reseq.fam)

full.fam <- bind_rows(fam.list)

duplicate.ids <- full.fam[duplicated(full.fam$V2), 2]
missing.uk <- uk.samples[!uk.samples$dna %in% full.fam$V2, 2]

bad.qual.samples <- list(read.delim("I:/psivakumar/Emer_CH/gwas/GSA2016_142_025_CC/bad_qual_samples.txt", colClasses = rep("character", 2)), 
                         read.delim("I:/psivakumar/Emer_CH/gwas/GSA2018_310_025-CC/badqual.txt", colClasses = rep("character", 2))) %>% 
  bind_rows() %>%
  .[[1]]

missing.uk.unexp <- missing.uk[!toupper(missing.uk) %in% bad.qual.samples]
bad.qual.uk.cases <- missing.uk[!missing.uk %in% missing.uk.unexp]

#write.table(duplicate.ids, "I:/psivakumar/Emer_CH/gwas/uk_analysis/duplicated_samples.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
#write.table(bad.qual.uk.cases, "I:/psivakumar/Emer_CH/gwas/uk_analysis/bad_quality_seq_not_included.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
#write.table(missing.uk.unexp, "I:/psivakumar/Emer_CH/gwas/uk_analysis/unexplained_missing_samples.txt.", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

uk.cases <- c(uk.samples[uk.samples$dna %in% full.fam$V2, 1], uk.reseq.samples) # %>% .[!. %in% duplicate.ids]
uk.controls <- bind_rows(c1958.fam, nbs.fam) %>% .[[2]]

#write.table(data.frame(uk.cases, uk.cases), "I:/psivakumar/Emer_CH/gwas/uk_analysis/ukCases.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
#write.table(data.frame(uk.controls, uk.controls), "I:/psivakumar/Emer_CH/gwas/uk_analysis/ukControls.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

# euro samples

setwd("I:/psivakumar/Emer_CH/gwas/euro_analysis")

spanish.cases <- fread("spanish_cases.csv", header = F, colClasses = "character")
belgian.cases <- fread("belgian_cases.csv", header = F, colClasses = "character")
danish.cases <- fread("danish_cases.csv", header = F, colClasses = "character")
greek.cases <- fread("greek_cases.csv", header = F, colClasses = "character")
greek.controls <- fread("greek_controls.csv", header = F, colClasses = "character")

gsa2019.fam <- fread("GSA2019_432_025_CC.fam", header = F, colClasses = c("character", "character", "integer", "integer", "integer", "integer"))
ich.controls.fam <- fread("GSA2016_142_025_CC_04-10.fam", header = F, colClasses = c("character", "character", "integer", "integer", "integer", "integer"))

euro.samples.list <- list(spanish.cases, belgian.cases, danish.cases, greek.cases, greek.controls)

missing.euro.samples <- lapply(euro.samples.list, function(x) x[!unlist(x[[1]]) %in% gsa2019.fam$V2])

duplicated.gsa2019.samples <- gsa2019.fam[duplicated(gsa2019.fam$V2), 1:2]
duplicated.ich.samples <- ich.controls.fam[duplicated(ich.controls.fam$V2), 1:2]
#fwrite(duplicated.case.samples, "duplicated_gsa2019_ids_remove.txt", col.names = F, sep = '\t')

#ich.renamed.fam <- fread("ICH_renamed.fam", header = F, colClasses = c("character", "character", "integer", "integer", "integer", "integer"))
#gsa2019.renamed.fam <- fread("GSA2019_renamed.fam", header = F, colClasses = c("character", "character", "integer", "integer", "integer", "integer"))

euro.cases <- unlist(bind_rows(spanish.cases, belgian.cases, danish.cases, greek.cases), use.names = F) %>% .[!. %in% unlist(missing.euro.samples, use.names = F)] %>% .[!duplicated(.)]
euro.controls <- unlist(bind_rows(greek.controls, ich.controls.fam[, "V2"]), use.names = F) %>% .[!. %in% unlist(missing.euro.samples, use.names = F)] %>% .[!duplicated(.)]

euro.cases.fam <- filter(gsa2019.fam, V2 %in% euro.cases)[, 1:2]
euro.controls.fam <- bind_rows(filter(ich.controls.fam, V2 %in% euro.controls)[, 1:2], filter(gsa2019.fam, V2 %in% euro.controls)[, 1:2])
#fwrite(euro.cases.fam, "euroCases.txt", col.names = F, sep = '\t')
#fwrite(euro.controls.fam, "euroControls.txt", col.names = F, sep = '\t')

all.ids <- c(euro.cases, euro.controls) #%>% .[!duplicated(.)]

gsa2019.samples.other <- gsa2019.fam[!gsa2019.fam$V2 %in% all.ids, 1:2]
#fwrite(gsa2019.samples.other, "unknown_gsa2019_ids_remove.txt", col.names = F, sep = '\t')

# merge uk and euro

fam.list <- list(gsa2016.fam, gsa2018.fam, gsa2019.fam, reseq.fam, nbs.fam, c1958.fam, ich.controls.fam)
full.fam <- bind_rows(fam.list)
all.cases <- c(uk.cases, euro.cases)
all.controls <- c(uk.controls, euro.controls)
all.ids <- c(all.cases, all.controls)

# check for monomorphic or ambiguous variants

setwd("I:/psivakumar/Emer_CH/gwas/euro_uk_analysis")

gsa2016.frq <- fread("alleleFreqs_GSA2016.frq")
gsa2018.frq <- fread("alleleFreqs_GSA2018.frq")
reseq.frq <- fread("alleleFreqs_reseq.frq") 
c1958.frq <- fread("alleleFreqs_1958c.frq")
nbs.frq <- fread("alleleFreqs_nbs.frq")
gsa2019.frq <- fread("alleleFreqs_gsa2019.frq")
ich.frq <- fread("alleleFreqs_ich.frq")
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

gsa2016.premerge.map <- fread("gsa2016_3_biallelic.map", header = F)
gsa2018.permerge.map <- fread("gsa2018_3_biallelic.map", header = F)
reseq.premerge.map <- fread("reseq_3_biallelic.map", header = F)
nbs.premerge.map <- fread("nbs_3_biallelic.map", header = F)
c1958.premerge.map <- fread("1958c_3_biallelic.map", header = F)
gsa2019.premerge.map <- fread("gsa2019_3_biallelic.map", header = F)
ich.premerge.map <- fread("ich_3_biallelic.map", header = F)

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

all.fam <- fread("all_5_merged.fam")
all.fam$V6 <- ifelse(all.fam$V2 %in% all.cases, 2, 1)
#all.fam$V1 <- all.fam$V2

#fwrite(all.fam, "all_5_merged.fam", col.names = F, sep = '\t')


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
