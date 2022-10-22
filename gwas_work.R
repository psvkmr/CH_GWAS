library(tidyverse)
library(bigreadr)
library(data.table)

# check all sample IDs are in plink files

setwd("I:/psivakumar/Emer_CH/gwas/uk_analysis")

reseq.samplesheet <- fread2("I:/psivakumar/Emer_CH/Plate 3_samplesheet.csv", skip = 16)
reseq.emer.samples <- fread2("I:/psivakumar/Emer_CH/Emer_resampled.csv")
uk.reseq.samples <- filter(reseq.samplesheet, Sample_ID %in% reseq.emer.samples$ID_2) %>%
  unite("id", c(SentrixBarcode_A, SentrixPosition_A), sep = "_") %>%
  .[['id']]

uk.samples <- read.csv("I:/psivakumar/Emer_CH/gwas/samplesheet_UK_only_final.csv", stringsAsFactors = FALSE)
fam_cols <- c(rep("character", 2), rep("integer", 4))
gsa2016.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2016_142_025_CC/PLINK_210318_0250/GSA2016_142_025_CC.fam", colClasses = fam_cols)
gsa2018.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2018_310_025-CC/PLINK_150319_1054/GSA2018_310_025-CC.fam", colClasses = fam_cols)
c1958.fam <- read.table("I:/psivakumar/Emer_CH/gwas/1958_controls/1958controls_14.09.2018.fam", colClasses = fam_cols)
nbs.fam <- read.table("I:/psivakumar/Emer_CH/gwas/NBS_controls/NBScontrols_14.09.2018.fam", colClasses = fam_cols)
reseq.fam <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/Plate_3_Emer.fam", colClasses = fam_cols)

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


# check map files

#gsa2016.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/GSA2016_142_025_CC.map")
#gsa2018.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/GSA2018_310_025-CC.map")
#c1958.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/1958controls_14.09.2018.map")
#nbs.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/NBScontrols_14.09.2018.map")

#preQC.map.list <- list(gsa2016.map, gsa2018.map, c1958.map, nbs.map)

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

#preQC.common.snps <- getCommonPos(lapply(preQC.map.list, addUniquePos))

#library(biomaRt)
#
#mart <- useMart("ENSEMBL_MART_SNP", host = "grch37.ensembl.org")
#dataset <- useDataset(mart = mart, dataset = "hsapiens_snp_som")
#snp.names <- getBM(mart = dataset, 
#                   attributes = c("refsnp_id", "refsnp_source", "chr_name", "chrom_start", "chrom_end", "allele", "minor_allele", "variation_names", "associated_gene"), 
#                   filters = c("chromosomal_region"), 
#                   values = common.snps[1])

# Use frq files to get monomorphic variants, and variants with ambiguous strand A/T C/G

#gsa2016.frq <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/alleleFreqs_GSA2016.frq", header = T)
#gsa2018.frq <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/alleleFreqs_GSA2018.frq", header = T)
#c1958.frq <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/alleleFreqs_1958c.frq", header = T)
#nbs.frq <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/alleleFreqs_nbs.frq", header = T)

#frq.list <- list(gsa2016.frq, gsa2018.frq, c1958.frq, nbs.frq)

#getMonomorphics <- function(df) {
#  filter(df, MAF < 0.0000001)$SNP
#}

#getAmbiguous <- function(df) {
#  filter(df, 
#         (A1 == "A" & A2 == "T") | 
#           (A1 == "T" & A2 == "A") |
#           (A1 == "C" & A2 == "G") | 
#           (A1 == "G" & A2 == "C"))$SNP
#}

#monomorphic.list <- lapply(frq.list, getMonomorphics)
#ambiguous.list <- lapply(frq.list, getAmbiguous)

#mono.or.ambig <- c(unlist(monomorphic.list), unlist(ambiguous.list)) %>% unique()

#write.table(unlist(monomorphic.list) %>% unique(), "I:/psivakumar/Emer_CH/gwas/uk_analysis/monomorphic_variants.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
#write.table(unlist(ambiguous.list) %>% unique(), "I:/psivakumar/Emer_CH/gwas/uk_analysis/ambigious_strand_variants", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

# get new map list after run_premerge_QC.sh has been run on each dataset

gsa2016.premerge.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/GSA2016_3_biallelic.map")
gsa2018.premerge.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/GSA2018_3_biallelic.map")
c1958.premerge.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/1958c_3_biallelic.map")
nbs.premerge.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/nbs_3_biallelic.map")
reseq.premerge.map <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/reseq_3_biallelic.map")

premerge.map.list <- list(gsa2016.premerge.map, gsa2018.premerge.map, c1958.premerge.map, nbs.premerge.map, reseq.premerge.map)

premerge.common.snps <- getCommonPos(lapply(premerge.map.list, addUniquePos))

#write.table(premerge.common.snps, "I:/psivakumar/Emer_CH/gwas/uk_analysis/commonSnpsToExtract.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)


# read in all samples fam file

all.fam <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/all_5_merged.fam")
all.fam$V6 <- ifelse(all.fam$V1 %in% uk.samples$dna | all.fam$V1 %in% uk.reseq.samples, 2, 1)

#write.table(all.fam, "I:/psivakumar/Emer_CH/gwas/uk_analysis/all_5_merged.fam", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = '\t')


# look at snps to be flipped

flipscan <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/postMerge.flipscan", colClasses = c("integer", "character", "integer", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "character"), header = T, fill = T)
non.empty.scan <- filter(flipscan, R_NEG > 0)

#write.table(non.empty.scan$SNP, "I:/psivakumar/Emer_CH/gwas/uk_analysis/strandIssueSnpsToRemove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# check imiss stats

cases.imiss <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/casesMissingStats.imiss", header = T)
controls.imiss <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/controlsMissingStats.imiss", header = T)

#check sex check

cases.sex <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/caseSex.sexcheck", header = T)
controls.sex <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/controlSex.sexcheck", header = T)
case.problem.sex <- filter(cases.sex, F < 0.6 & F > 0.4)
control.problem.sex <- filter(controls.sex, F < 0.6 & F > 0.4)
#fwrite(case.problem.sex[, 1:2], "caseFailSexCheck.txt", col.names = F, sep = '\t')
#fwrite(control.problem.sex[, 1:2], "controlFailSexCheck.txt", col.names = F, sep = '\t')


# check re-merged allele freqs
control.recalc.freq <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/controlRecalcFreq.frq", header = T)
case.recalc.freq <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/casesRecalcFreq.frq", header = T)
all.recalc.freq <- left_join(control.recalc.freq, case.recalc.freq, by = "SNP") %>%
  mutate(DIFF = log10(MAF.y) - log10(MAF.x)) %>%
  arrange(desc(DIFF))


# check het
all.hets <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/allHets.het", header = T)
all.hets$HET_RATE = (all.hets$N.NM. - all.hets$O.HOM.) / all.hets$N.NM.
#het.fail = subset(all.hets, (all.hets$HET_RATE < mean(all.hets$HET_RATE)-3*sd(all.hets$HET_RATE)) | (all.hets$HET_RATE > mean(all.hets$HET_RATE)+3*sd(all.hets$HET_RATE)))
#het.fail$HET_DST = (het.fail$HET_RATE - mean(all.hets$HET_RATE)) / sd(all.hets$HET_RATE)
#het.fail <- het.fail[, 1:2]
het.fail <- filter(all.hets, HET_RATE > 0.5 | HET_RATE < 0.1)
#write.table(het.fail, "I:/psivakumar/Emer_CH/gwas/uk_analysis/hetOutliersToRemove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


# check fam file for duplicate ids
all_20_fam <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/all_20_noDups.fam")


# peddy pop strat

pop.strat <- read.csv("I:/psivakumar/Emer_CH/gwas/uk_analysis/peddy_gwas.het_check.csv")
pop.dist <- table(pop.strat$ancestry.prediction)
european <- filter(pop.strat, ancestry.prediction == "EUR" | ancestry.prediction == "UNKNOWN")$sample_id
european.to.keep <- filter(all_20_fam, V2 %in% european)[, 1:2]
#write.table(european.to.keep, "I:/psivakumar/Emer_CH/gwas/uk_analysis/european_to_keep.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# pca

eigenvecs <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/pca.eigenvec") 
eigenvals <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/pca.eigenval")

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
pc.outliers <- lapply(eigenvecs[, pca.filter.columns], detectOutliers)
pca.outlier.samples <- lapply(pc.outliers, function(x) eigenvecs[unlist(x), 2])
pca.all.outlier.samples <- unique(unlist(pca.outlier.samples))
#pca.outlier.samples.df <- data.frame(fid = pca.all.outlier.samples, iid = pca.all.outlier.samples)
# not excluding outliers for now
pca.outlier.samples.df <- data.frame()
#write.table(pca.outlier.samples.df, "I:/psivakumar/Emer_CH/gwas/uk_analysis/pcaSamplesToRemove.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')

# check batch effects

info.fam <- all.fam %>%
  mutate(pheno = ifelse(V2 %in% uk.samples$dna, 2, 1)) %>%
  mutate(batch = ifelse(V2 %in% gsa2016.fam$V2, yes = "gsa2016", no = 
                          ifelse(V2 %in% gsa2018.fam$V2, yes = "gsa2018", no = 
                                   ifelse(V2 %in% c1958.fam$V2, yes = "c1958", no = 
                                            ifelse(V2 %in% nbs.fam$V2, yes = "nbs", no = "other")))))


# pre and post pop  strat/pca qc

pre.qc.pca <- left_join(eigenvecs, info.fam, by = "V2") %>%
  ggplot(aes(V3.x, V4.x, colour = as.factor(batch))) +
  geom_point(size=2, alpha = 0.5) +
  theme_classic()

post.qc.pca <- left_join(eigenvecs, info.fam, by = "V2") %>%
  .[!.$V2 %in% pca.all.outlier.samples, ] %>%
  ggplot(aes(V3.x, V4.x, colour = as.factor(batch), alpha = 0.1)) + 
  geom_point() +
  theme_classic()


# SAIGE
# pheno 

final.samples <- fread("ID_order_chrX.txt", header = F)
final.cases <- dplyr::filter(all_20_fam, V2 %in% final.samples$V1, V6 == 2)[, 1:2]
#fwrite(final.cases, "final_cases.csv")
saige.pheno <- left_join(all_20_fam, eigenvecs, by = "V2") %>%
  filter(V2 %in% final.samples$V1) %>%
  mutate(new_pheno = ifelse(V6.x == 2, 1, ifelse(V6.x == 1, 0, NA))) %>%
  select(FID = V1.x, IID = V2, PATID = V3.x, MATID = V4.x, SEX = V5.x, Phenotype = new_pheno, PC1 = V3.y, PC2 = V4.y, PC3 = V5.y, PC4 = V6.y, PC5 = V7, PC6 = V8, PC7 = V9, PC8 = V10, PC9 = V11, PC10 = V12, PC11 = V13, PC12 = V14, PC13 = V15, PC14 = V16, PC15 = V17, PC16 = V18, PC17 = V19, PC18 = V20, PC19 = V21, PC20 = V22)
#fwrite(saige.pheno, "saige.pheno", sep = '\t')
