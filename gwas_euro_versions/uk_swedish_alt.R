library(tidyverse)
library(data.table)

setwd("I:/psivakumar/Emer_CH/gwas/uk_swedish_alt")

reseq.samplesheet <- fread("I:/psivakumar/Emer_CH/Plate 3_samplesheet.csv", skip = 16)
reseq.emer.samples <- fread("I:/psivakumar/Emer_CH/Emer_resampled.csv")
reseq.emer.samplesheet <- filter(reseq.samplesheet, Sample_ID %in% reseq.emer.samples$ID_2) %>%
  unite("id", c(SentrixBarcode_A, SentrixPosition_A), sep = "_", remove = F)

uk.samples <- fread("I:/psivakumar/Emer_CH/gwas/samplesheet_UK_only_final.csv", stringsAsFactors = FALSE, header = T)
fam_cols <- c(rep("character", 2), rep("integer", 4))
gsa2016.fam <- fread("I:/psivakumar/Emer_CH/gwas/GSA2016_142_025_CC/PLINK_210318_0250/GSA2016_142_025_CC.fam", colClasses = fam_cols)
gsa2018.fam <- fread("I:/psivakumar/Emer_CH/gwas/GSA2018_310_025-CC/PLINK_150319_1054/GSA2018_310_025-CC.fam", colClasses = fam_cols)
c1958.fam <- fread("I:/psivakumar/Emer_CH/gwas/1958_controls/1958controls_14.09.2018.fam", colClasses = fam_cols)
nbs.fam <- fread("I:/psivakumar/Emer_CH/gwas/NBS_controls/NBScontrols_14.09.2018.fam", colClasses = fam_cols)
reseq.fam <- fread("I:/psivakumar/Emer_CH/gwas/Plate_3_Emer/Plate_3_Emer.fam", colClasses = fam_cols)
swedish.fam <- fread("I:/psivakumar/Emer_CH/gwas/swedish_data/Dataset_UK_newIDord.fam", colClasses = fam_cols)

fam.list <- list(c1958.fam, gsa2016.fam, gsa2018.fam, nbs.fam, reseq.fam, swedish.fam)

full.fam <- bind_rows(fam.list)

duplicate.ids <- full.fam[duplicated(full.fam$V2), ]
missing.uk <- uk.samples[!uk.samples$dna %in% full.fam$V2, 2]

bad.qual.samples <- list(read.delim("I:/psivakumar/Emer_CH/gwas/GSA2016_142_025_CC/bad_qual_samples.txt", colClasses = rep("character", 2)), 
                         read.delim("I:/psivakumar/Emer_CH/gwas/GSA2018_310_025-CC/badqual.txt", colClasses = rep("character", 2))) %>% 
  bind_rows()

#missing.uk.unexp <- missing.uk[!toupper(missing.uk) %in% bad.qual.samples]
#bad.qual.uk.cases <- missing.uk[!missing.uk %in% missing.uk.unexp]

#fwrite(duplicate.ids[, "V2"], "duplicated_samples.txt", col.names = F)
#fwrite(bad.qual.uk.cases, "bad_quality_seq_not_included.txt", col.names = F)
#fwrite(missing.uk.unexp, "unexplained_missing_samples.txt.", col.names = F)

all.samples <- full.fam$V2
swedish.samples <- swedish.fam$V2
swedish.cases <- filter(swedish.fam, V6 == 2)$V2
swedish.controls <- filter(swedish.fam, V6 == 1)$V2
uk.samples.vec <- all.samples[!all.samples %in% swedish.samples]
uk.cases <- c(filter(uk.samples, dna %in% full.fam$V2)$dna, filter(reseq.emer.samplesheet, Sample_ID %in% full.fam$V2)$Sample_ID)
uk.controls <- uk.samples.vec[!uk.samples.vec %in% uk.cases]
all.cases <- c(uk.cases, swedish.cases)
all.controls <- c(uk.controls, swedish.controls)

#fwrite(as.data.frame(uk.cases), "uk_cases_pre_QC.txt", col.names = F)
#fwrite(as.data.frame(uk.controls), "uk_controls_pre_QC.txt", col.names = F)


# Use frq files to get monomorphic variants, and variants with ambiguous strand A/T C/G

gsa2016.frq <- fread("alleleFreqs_GSA2016.frq", header = T)
gsa2018.frq <- fread("alleleFreqs_GSA2018.frq", header = T)
c1958.frq <- fread("alleleFreqs_1958c.frq", header = T)
nbs.frq <- fread("alleleFreqs_nbs.frq", header = T)
reseq.frq <- fread("alleleFreqs_reseq.frq", header = T)
swedish.frq <- fread("alleleFreqs_swedish.frq", header = T)

frq.list <- list(gsa2016.frq, gsa2018.frq, c1958.frq, nbs.frq, reseq.frq, swedish.frq)

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

mono.snps <- unlist(monomorphic.list) %>% unique()
ambig.snps <- unlist(ambiguous.list) %>% unique()
mono.or.ambig <- c(unlist(monomorphic.list), unlist(ambiguous.list)) %>% unique()

#fwrite(as.data.frame(mono.snps), "monomorphic_variants.txt", col.names = F)
#fwrite(as.data.frame(ambig.snps), "ambigious_strand_variants", col.names = F)


# get common snps

getOtherChr <- function(df) {
  snps <- df[df$V1 %in% c(0, 24, 25, 26), 2]
}

#other.chr.snps <- lapply(preQC.map.list, getOtherChr)

addUniquePos <- function(df){
  mutate(df, V5 = paste0(V1, ":", V4, "-", V4))
}

getCommonPos <- function(lst) {
  compos <- Reduce(intersect, lapply(lst, `[[`, 5))
  snplst <- lapply(lst, function(x) x[which(x$V5 %in% compos), 2])
  comsnp <- Reduce(intersect, snplst)
}

gsa2016.premerge.map <- fread("GSA2016_3_biallelic.map")
gsa2018.premerge.map <- fread("GSA2018_3_biallelic.map")
c1958.premerge.map <- fread("1958c_3_biallelic.map")
nbs.premerge.map <- fread("nbs_3_biallelic.map")
reseq.premerge.map <- fread("reseq_3_biallelic.map")
swedish.premerge.map <- fread("swedish_3_biallelic.map")

premerge.map.list <- list(gsa2016.premerge.map, gsa2018.premerge.map, c1958.premerge.map, nbs.premerge.map, reseq.premerge.map, swedish.premerge.map)

premerge.common.snps <- getCommonPos(lapply(premerge.map.list, addUniquePos))

#fwrite(as.data.frame(premerge.common.snps), "commonSnpsToExtract.txt", col.names = F)


# read in all samples fam file

all.5.fam <- fread("all_5_merged.fam")
all.5.df <- separate(all.5.fam, V2, c("SET", "FID", "IID"), sep = ":", remove = F)
# 9 duplicated samples removed during merge (when in same dataset, not from reseq)
all.5.duplicated <- all.5.df[duplicated(all.5.df$IID), ]
full.fam.duplicated <- full.fam[duplicated(full.fam$V2), ]

all.5.fam$V6 <- ifelse(all.5.fam$V1 %in% all.cases, 2, 
                       ifelse(all.5.fam$V1 %in% all.controls, 1, -9))

#fwrite(all.5.fam, "all_5_merged.fam", col.names = F, sep = '\t')


# look at snps to be flipped

flipscan <- fread("postMerge.flipscan", colClasses = c("integer", "character", "integer", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "character"), header = T, fill = T)
non.empty.scan <- filter(flipscan, R_NEG > 0)

#fwrite(non.empty.scan["SNP"], "strandIssueSnpsToRemove.txt", col.names = F)


# check imiss stats

cases.imiss <- fread("casesMissingStats.imiss", header = T)
controls.imiss <- fread("controlsMissingStats.imiss", header = T)

#check sex check

cases.sex <- fread("caseSex.sexcheck", header = T)
controls.sex <- fread("controlSex.sexcheck", header = T)
case.problem.sex <- filter(cases.sex, F < 0.6 & F > 0.4)
control.problem.sex <- filter(controls.sex, F < 0.6 & F > 0.4)
#fwrite(case.problem.sex[, 1:2], "caseFailSexCheck.txt", col.names = F, sep = '\t')
#fwrite(control.problem.sex[, 1:2], "controlFailSexCheck.txt", col.names = F, sep = '\t')


# check re-merged allele freqs
control.recalc.freq <- fread("controlRecalcFreq.frq", header = T)
case.recalc.freq <- fread("casesRecalcFreq.frq", header = T)
all.recalc.freq <- left_join(control.recalc.freq, case.recalc.freq, by = "SNP") %>%
  mutate(DIFF = log10(MAF.y) - log10(MAF.x)) %>%
  arrange(desc(DIFF))


# check het
all.hets <- fread("allHets.het", header = T)
all.hets$HET_RATE = (all.hets$`N(NM)` - all.hets$`O(HOM)`) / all.hets$`N(NM)`
#het.fail = subset(all.hets, (all.hets$HET_RATE < mean(all.hets$HET_RATE)-3*sd(all.hets$HET_RATE)) | (all.hets$HET_RATE > mean(all.hets$HET_RATE)+3*sd(all.hets$HET_RATE)))
#het.fail$HET_DST = (het.fail$HET_RATE - mean(all.hets$HET_RATE)) / sd(all.hets$HET_RATE)
#het.fail <- het.fail[, 1:2]
het.fail <- filter(all.hets, HET_RATE > 0.5 | HET_RATE < 0.1)
#fwrite(het.fail[, 1:2], "hetOutliersToRemove.txt", col.names = F, sep = '\t')


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

pca.filter.columns <- c(3:6)    
pc.outliers <- lapply(eigenvecs[, ..pca.filter.columns], detectOutliers)
pca.outlier.samples <- lapply(pc.outliers, function(x) eigenvecs[unlist(x), 2])
pca.all.outlier.samples <- unique(unlist(pca.outlier.samples))
#pca.outlier.samples.df <- data.frame(fid = pca.all.outlier.samples, iid = pca.all.outlier.samples)
# not excluding outliers for now
pca.outlier.samples.df <- data.frame(NA)
#fwrite(pca.outlier.samples.df, "pcaSamplesToRemove.txt", col.names = F, sep = '\t')

# check batch effects

info.fam <- all.20.fam %>%
  separate(V2, c("BATCH", "FID", "IID"), sep = ":", remove = F)


# pre and post pop  strat/pca qc

pre.qc.pca <- left_join(eigenvecs, info.fam, by = "V2") %>%
  ggplot(aes(V3.x, V4.x, colour = as.factor(BATCH))) +
  geom_point(size=2, alpha = 0.5) +
  theme_classic()

post.qc.pca <- left_join(eigenvecs, info.fam, by = "V2") %>%
  .[!.$V2 %in% pca.all.outlier.samples, ] %>%
  ggplot(aes(V3.x, V4.x, colour = as.factor(BATCH), alpha = 0.1)) + 
  geom_point() +
  theme_classic()

final.fam <- fread("I:/psivakumar/Emer_CH/gwas/uk_swedish_alt/all_24_hhMissing-updated.fam")

# SAIGE
# pheno 

final.samples <- fread("ID_order_chrX.txt", header = F)
final.cases <- dplyr::filter(all.20.fam, V2 %in% final.samples$V1, V6 == 2)[, 1:2]
#fwrite(final.cases, "final_cases.csv")
saige.pheno <- left_join(all.20.fam, eigenvecs, by = "V2") %>%
  filter(V2 %in% final.samples$V1) %>%
  mutate(new_pheno = ifelse(V6.x == 2, 1, ifelse(V6.x == 1, 0, NA))) %>%
  select(FID = V1.x, IID = V2, PATID = V3.x, MATID = V4.x, SEX = V5.x, Phenotype = new_pheno, PC1 = V3.y, PC2 = V4.y, PC3 = V5.y, PC4 = V6.y, PC5 = V7, PC6 = V8, PC7 = V9, PC8 = V10, PC9 = V11, PC10 = V12, PC11 = V13, PC12 = V14, PC13 = V15, PC14 = V16, PC15 = V17, PC16 = V18, PC17 = V19, PC18 = V20, PC19 = V21, PC20 = V22)
#fwrite(saige.pheno, "saige.pheno", sep = '\t')
