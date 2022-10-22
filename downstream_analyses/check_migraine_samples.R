library(data.table)
library(tidyverse)

fam_cols <- c(rep("character", 2), rep("integer", 4))
gsa2016.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2016_142_025_CC/PLINK_210318_0250/GSA2016_142_025_CC.fam", colClasses = fam_cols)
gsa2018.fam <- read.table("I:/psivakumar/Emer_CH/gwas/GSA2018_310_025-CC/PLINK_150319_1054/GSA2018_310_025-CC.fam", colClasses = fam_cols)
c1958.fam <- read.table("I:/psivakumar/Emer_CH/gwas/1958_controls/1958controls_14.09.2018.fam", colClasses = fam_cols)
nbs.fam <- read.table("I:/psivakumar/Emer_CH/gwas/NBS_controls/NBScontrols_14.09.2018.fam", colClasses = fam_cols)
reseq.fam <- read.table("I:/psivakumar/Emer_CH/gwas/uk_analysis/Plate_3_Emer.fam", colClasses = fam_cols)

sample.var.sums <- fread("I:/psivakumar/Emer_CH/gwas/uk_analysis/gwas_sig_hits_samples_summary.csv")
migraine.samples <- fread("I:/psivakumar/Emer_CH/gwas/all_case_samples_migraine_status.csv")

fam.list <- list(gsa2016.fam, gsa2018.fam, c1958.fam, nbs.fam, reseq.fam)
fam.source.names <- c("GSA2016", "GSA2018", "1958c", "nbs", "reseq")

for (i in 1:length(fam.source.names)){
 fam.list[[i]] <- mutate(fam.list[[i]], source = fam.source.names[i])
}

full.fam <- bind_rows(fam.list)
full.fam.ids <- unite(full.fam, "new_ids", c("source", V1, V2), sep = "_", remove = F)

migraine.samples <- mutate(migraine.samples, source = ifelse(ARRAY == "GSA2018_310_025-CC", "GSA2018", 
                                                             ifelse(ARRAY == "GSA2016_142_025_CC", "GSA2016", 
                                                                    ifelse(ARRAY == "Plate_3_Emer", "reseq", "other"))))
migraine.samples$Migraine <- gsub("#N/A", NA, migraine.samples$Migraine)
migraine.samples$Migraine <- gsub("na", NA, migraine.samples$Migraine)
migraine.samples$Migraine <- gsub("n", 0, migraine.samples$Migraine)
migraine.samples$Migraine <- gsub("N", "0", migraine.samples$Migraine)
migraine.samples$Migraine <- gsub("y", 1, migraine.samples$Migraine)
migraine.samples$Migraine <- gsub("Y", 1, migraine.samples$Migraine)

full.fam.migraine.status <- left_join(migraine.samples, full.fam.ids, by = c("source", "SAMPLE" = "V2"))

migraine.pos <- filter(full.fam.migraine.status, Migraine == 1, source != "other")$new_ids
migraine.neg <- filter(full.fam.migraine.status, Migraine == 0, source != "other")$new_ids
migraine.na <- filter(full.fam.migraine.status, is.na(Migraine))$new_ids

# var.int <- "6:97058553:A:G"
var.int.mig.checker <- function(var_int) {
  var.int <- var_int
  var.int.sums <- filter(sample.var.sums, SNP == var.int)
  var.int.case.ref.names <- str_split(var.int.sums[[1, 11]], pattern = " ") %>% unlist(use.names = F)
  var.int.case.alt.names <- str_split(var.int.sums[[1, 13]], pattern = " ") %>% unlist(use.names = F)
  var.int.case.ref.mig.neg.count <- length(intersect(migraine.neg, var.int.case.ref.names))
  var.int.case.alt.mig.neg.count <- length(intersect(migraine.neg, var.int.case.alt.names))
  var.int.case.ref.mig.pos.count <- length(intersect(migraine.pos, var.int.case.ref.names))
  var.int.case.alt.mig.pos.count <- length(intersect(migraine.pos, var.int.case.alt.names))
  var.int.mat <- matrix(c(var.int.case.ref.mig.neg.count, var.int.case.alt.mig.neg.count, var.int.case.ref.mig.pos.count, var.int.case.alt.mig.pos.count), nrow = 2, byrow = T)
  var.int.fishers <- fisher.test(var.int.mat)
  var.int.stats <- list(var.int.mat, var.int.fishers)
  return(var.int.stats)
}

case.names <- c(sample.var.sums[[1, 11]], sample.var.sums[[1, 13]]) %>% str_split(pattern = " ") %>% unlist()
gwas.cases.not.in.migraine.file <- case.names[!case.names %in% full.fam.migraine.status$new_ids]
mig.cases.not.in.gwas <- migraine.samples[!migraine.samples$SAMPLE %in% full.fam$V2, ]$SAMPLE

unimp.res <- fread("I:/psivakumar/Emer_CH/gwas/uk_analysis/unimputed_fishers/unimputed_sig.csv")
unimp.res.sig.snps <- sample.var.sums[grep(paste(unimp.res$POS, collapse = "|"), sample.var.sums$SNP), ]$SNP
unimp.res.mig.res <- lapply(unimp.res.sig.snps, var.int.mig.checker)
names(unimp.res.mig.res) <- unimp.res.sig.snps

sig.res.mig.res <- lapply(sig.res$SNPID, var.int.mig.checker)

fishers.res.to.vec <- function(fisher_res){
  ctrl_ref <- fisher_res[[1]][[1]]
  case_ref <- fisher_res[[1]][[2]]
  ctrl_alt <- fisher_res[[1]][[3]]
  case_alt <- fisher_res[[1]][[4]]
  p <- fisher_res[[2]][[1]]
  or <- fisher_res[[2]][[3]][[1]]
  lci <- fisher_res[[2]][[2]][[1]]
  uci <- fisher_res[[2]][[2]][[2]]
  res <- data.frame(ctrl_ref, ctrl_alt, case_ref, case_alt, p, or, lci, uci)
  return(res)
}

sig.res.mig.res.df <- bind_rows(lapply(sig.res.mig.res, fishers.res.to.vec)) %>% cbind(sig.res$SNPID)
sig.res.mig.sig.res <- filter(sig.res.mig.res.df, p < 0.05)
#fwrite(sig.res.mig.res.df, "I:/psivakumar/Emer_CH/gwas/uk_analysis/migraine_vs_nonmigraine_in_gwas_sig_res_snps.csv")

# check ref alt order diff between vcfs and saige output, maybe unimp plink is opposite
