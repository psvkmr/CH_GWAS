library(snpStats)
library(coloc)

setwd("I:/psivakumar/Emer_CH/gwas/")

euro.bed <- "euro_analysis/all_24_hhMissing.bed"
euro.bim <- "euro_analysis/all_24_hhMissing.bim"
euro.fam <- "euro_analysis/all_24_hhMissing.fam"

euro.all.24 <- read.plink(euro.bed, euro.bim, euro.fam)

uk.bed <- "uk_analysis/strict_pca/all_24_hhMissing.bed"
uk.bim <- "uk_analysis/strict_pca/all_24_hhMissing.bim"
uk.fam <- "uk_analysis/strict_pca/all_24_hhMissing.fam"

uk.all.24 <- read.plink(uk.bed, uk.bim, uk.fam)

euro.pheno <- euro.all.24[[2]]$affected
uk.pheno <- uk.all.24[[2]]$affected

euro.df <- as.data.frame(cbind(Y = euro.pheno, as(euro.all.24[[1]], "numeric")))
uk.df <- as.data.frame(cbind(Y = uk.pheno, as(uk.all.24[[1]], "numeric")))

ct.abf <- coloc.abf.datasets(uk.df, euro.df, response1="Y", response2="Y", type1="quant", type2="quant")

uk.df.pc <- uk.df[, -1]
uk.df.pc <- uk.df.pc[, which(apply(uk.df.pc, 2, var) !=0)]

euro.df.pc <- euro.df[, -1]
euro.df.pc <- euro.df.pc[, which(apply(euro.df.pc, 2, var) !=0)]

pcs <- pcs.prepare(uk.df.pc, euro.df.pc)
pcs.1 <- pcs.model(pcs, group=1, Y=uk.df[,1], threshold = 0.8)
pcs.2 <- pcs.model(pcs, group = 2, Y = euro.df[,1], threshold = 0.8)
ct.pcs <- coloc.test(pcs.1, pcs.2)
ct.pcs.sum <- str(summary(ct.pcs))
ct.pcs.bayes <- coloc.test(pcs.1, pcs.2, bayes=TRUE)
ct.pcs.bayes.ci <- ci(ct.pcs.bayes)
