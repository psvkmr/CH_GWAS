# cases
# CH_GWAS contains 2016 data, lots of idats, maybe plink files but no permissions
# 2018 contains 2018 data, lots of idats again but also an xml file and .sdf file, also CC directory with QC reports, Data etc, also plink files?
# 2018 also contains 2019 data, no UK
# UK only in CH_GWAS - 2016, 2018
# all folders (2016, 2018)? also contain 'control' ie Isobel cases
# 2019 folder contains all non-UK, cases and controls for most

# controls
# 1958 contains vcf files split by chromosome, and also plink files, and also bgen files with script to convert to vcf
# ICH_controls contains plink files and some resulting QC outputs
# NBS contains vcf files - labelled imputed - split by chromosome, and some sex chr plink files and bgen files
# 1958 and NBS Wellcome controls are UK, look out for batch effect
# controls_raw WTCCC_nonlmp.zip
# ICH batch effect controls in ICH_controls

# gwas guide suggests starting with plink files, find where originals are
# alternatively, find raw data types and convert to plink
# do they have the illumina final report files?

# created binary files from plink text files for each of 2016, 2018, 2019
# in uk_analysis, re-created text files for each plus controls

# prior to merge:
# set ids uniquely with plink2 set-all-var-ids may have to use bfile '@:#$r:$a'
# only keep chrosomosomes 1-23
# remove duplicates with --rm-dup
# keep only biallelic with biallelic only
# find monomorphic alleles with maf == 0 in frq
# remove A/T and C/G variants with those in frq?
# snps-only 'just-acgt' ?


# run run_premerge_QC.sh on each of the 2016, 2018, NBS and 1958 datasets
bash run_premerge_QC.sh GSA2016_142_025_CC GSA2016
bash run_premerge_QC.sh GSA2018_310_025-CC GSA2018
bash run_premerge_QC.sh 1958controls_14.09.2018 1958c
bash run_premerge_QC.sh NBScontrols_14.09.2018 nbs
bash run_premerge_QC.sh Plate_3_Emer reseq
# In R, intersect maps to create common snp list, called commonSnpsToExtract.txt, and extract
plink --bfile GSA2016_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out GSA2016_4_commonSnps
plink --bfile GSA2018_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out GSA2018_4_commonSnps
plink --bfile 1958c_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out 1958c_4_commonSnps
plink --bfile nbs_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out nbs_4_commonSnps
plink --bfile reseq_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out reseq_4_commonSnps
# rename samples
awk '{print $2,"GSA2016_"$1"_"$2,$3,$4,$5,$6}' GSA2016_4_commonSnps.fam > rename_GSA2016.txt
mv rename_GSA2016.txt GSA2016_4_commonSnps.fam
#plink --bfile GSA2016_4_commonSnps --update-ids rename_GSA2016.txt --make-bed --out GSA2016_5_newIds
awk '{print $2,"GSA2018_"$1"_"$2,$3,$4,$5,$6}' GSA2018_4_commonSnps.fam > rename_GSA2018.txt
mv rename_GSA2018.txt GSA2018_4_commonSnps.fam
awk '{print $2,"1958c_"$1"_"$2,$3,$4,$5,$6}' 1958c_4_commonSnps.fam > rename_1958c.txt
mv rename_1958c.txt 1958c_4_commonSnps.fam
awk '{print $2,"nbs_"$1"_"$2,$3,$4,$5,$6}' nbs_4_commonSnps.fam > rename_nbs.txt
mv rename_nbs.txt nbs_4_commonSnps.fam
awk '{print $2,"reseq_"$1"_"$2,$3,$4,$5,$6}' reseq_4_commonSnps.fam > rename_reseq.txt
mv rename_reseq.txt reseq_4_commonSnps.fam
# merge
plink --bfile GSA2016_4_commonSnps --bmerge GSA2018_4_commonSnps --merge-mode 4 --out tmp1
plink --bfile tmp1 --bmerge 1958c_4_commonSnps --merge-mode 4 --out tmp2
plink --bfile tmp2 --bmerge nbs_4_commonSnps --merge-mode 4 --out tmp3
plink --bfile tmp3 --bmerge reseq_4_commonSnps --merge-mode 4 --out all_5_merged
plink --bfile all_5_merged --recode --out all_5_merged
# edit fam file to include phenotypes in R
# check for any strand errors
plink --bfile all_5_merged --allow-no-sex --flip-scan --out postMerge
# create list of strand issue variants to remove in R
plink --bfile all_5_merged --exclude strandIssueSnpsToRemove.txt --allow-no-sex --make-bed --recode --out all_6_noStrandIssue
# split by case and control
#plink --bfile all_6_noStrandIssue --keep ukCases.txt --make-bed --recode --out cases_7_split
#plink --bfile all_6_noStrandIssue --keep ukControls.txt --make-bed --recode --out controls_7_split
# split by case and control
awk '$6==2{print $1,$2}' all_6_noStrandIssue.fam > ukCases.txt
awk '$6==1{print $1,$2}' all_6_noStrandIssue.fam > ukControls.txt
plink --bfile all_6_noStrandIssue --keep ukCases.txt --make-bed --recode --out cases_7_split
plink --bfile all_6_noStrandIssue --keep ukControls.txt --make-bed --recode --out controls_7_split
# generate missingness stats
plink --bfile controls_7_split --missing --allow-no-sex --out controlsMissingStats
plink --bfile cases_7_split --missing --allow-no-sex --out casesMissingStats
# remove high missing SNPs
plink --bfile controls_7_split --geno 0.05 --make-bed --out controls_8_snps95
plink --bfile cases_7_split --geno 0.05 --make-bed --out cases_8_snps95
# remove any high missing ind
plink --bfile controls_8_snps95 --mind 0.02 --make-bed --out controls_9_ind98
plink --bfile cases_8_snps95 --mind 0.02 --make-bed --out cases_9_ind98
# sex check
plink --bfile controls_9_ind98 --check-sex --out controlSex
plink --bfile cases_9_ind98 --check-sex --out caseSex
# get list of problem sex
# F statistic dist looked fine, didn't remove any samples
# grep PROBLEM controlSex.sexcheck | awk '{print $1,$2}' > controlFailSexCheck.txt
# grep PROBLEM caseSex.sexcheck | awk '{print $1,$2}' > caseFailSexCheck.txt
# remove sex problem samples
plink --bfile controls_9_ind98 --remove controlFailSexCheck.txt --make-bed --out controls_10_sexChecked
plink --bfile cases_9_ind98 --remove caseFailSexCheck.txt --make-bed --out cases_10_sexChecked
# impute sex instead
plink --bfile controls_9_ind98 --impute-sex --make-bed --out controls_10_sexChecked
plink --bfile cases_9_ind98 --impute-sex --make-bed --out cases_10_sexChecked
# set higher missing snp threshold
plink --bfile controls_10_sexChecked --geno 0.02 --make-bed --out controls_11_snps98
plink --bfile cases_10_sexChecked --geno 0.02 --make-bed --out cases_11_snps98
# generate hwe stats
plink --bfile controls_11_snps98 --hardy --out controlsHWEcheck
plink --bfile cases_11_snps98 --hardy --out casesHWEcheck
# remove hwe deviants
plink --bfile controls_11_snps98 --hwe 1e-6 --make-bed --out controls_12_hwe
plink --bfile cases_11_snps98 --hwe 1e-6 --make-bed --out cases_12_hwe
#generate allele freq stats
plink --bfile controls_11_snps98 --freq --out controlRecalcFreq
plink --bfile cases_11_snps98 --freq --out casesRecalcFreq
# re-merge
plink --bfile controls_12_hwe --bmerge cases_12_hwe --allow-no-sex --make-bed --out all_13_remerge
# snps 98 over all
plink --bfile all_13_remerge --geno 0.02 --make-bed --out all_14_snps98
# missing ind 98 all
plink --bfile all_14_snps98 --mind 0.02 --make-bed --out all_15_ind98
# check no monomorphic
plink --bfile all_15_ind98 --maf 0.0000001 --make-bed --out all_16_noMonomoprhs
# remove more hwe
plink --bfile all_16_noMonomoprhs --hwe 1e-6 --make-bed --out all_17_hwe
# remove het outliers
plink --bfile all_17_hwe --het --out allHets
# In R, find samples with outlier number of hets
# Remove
plink --bfile all_17_hwe --remove hetOutliersToRemove.txt --make-bed --out all_18_noHetOutliers
# prune for ld
# got set of high ld regions from https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
awk '$4=(FNR FS $4)' highLDRegions.txt > highLDRegionsIds.txt
plink --bfile all_18_noHetOutliers --make-set highLDRegionsIds.txt --write-set --out highLDR
plink --bfile all_18_noHetOutliers --exclude highLDR.set  --make-bed --out all_19_noHighLD
plink --bfile all_19_noHighLD --indep-pairwise 50 5 0.2 --out LDprune
# prune and calculate relatedness
plink --bfile all_19_noHighLD --extract LDprune.prune.in --genome --out relCheck
# remove duplicates/monozygotic twins
awk '$10>0.2 {print $2,$4,$10}' relCheck.genome > duplicateInfo.txt
awk '$10>0.2 {print $1,$2}' relCheck.genome | sed '1d' > duplicateSamples.txt
# maybe should do from all_18_noHetOutliers??
plink --bfile all_19_noHighLD --remove duplicateSamples.txt --make-bed --out all_20_noDups
# use peddy for pop strat
plink --bfile all_20_noDups --recode --out all_20_noDups
plink --bfile all_20_noDups --keep-allele-order --recode vcf-iid --out for_peddy
bgzip for_peddy.vcf
tabix -p vcf for_peddy.vcf.gz
python3.6 -m peddy -p 4 --plot --prefix peddy_gwas for_peddy.vcf.gz all_20_noDups.ped
# create list of non-european to exclude in R
plink --bfile all_20_noDups --keep european_to_keep.txt --make-bed --out all_21_european
# pca
plink --bfile all_21_european --pca 20 --out pca
# In R, get list of samples which are outliers in any of first 20 PCs
# nonclean analysis: no PC outliers removed, first 4 PCs covariates
# alt pca: first 4 PC outliers removed, PCs 5-8 as covariates
# strict pca : 20 PC outliers removed, so no covariates
# remove
plink --bfile all_21_european --remove pcaSamplesToRemove.txt --make-bed --out all_22_pcaCleaned
# genotyping differences
plink --bfile all_22_pcaCleaned --test-missing --out genotypingDiffs
# remove significant different snps
awk '$5<0.00001 {print $2}' genotypingDiffs.missing > sigDiffGenoSnps.txt
plink --bfile all_22_pcaCleaned --exclude sigDiffGenoSnps.txt --make-bed --recode --out all_23_noSigDiffGeno
# manually exclude questionable snps with silly pvalue, no reason why yet
plink --bfile all_23_noSigDiffGeno --exclude excessPvalSnpsRemove.txt --make-bed --out all_23_noSigDiffGeno
# set hh to missing
plink --bfile all_23_noSigDiffGeno --set-hh-missing --make-bed --out all_24_hhMissing
# get final freq file
plink --bfile all_24_hhMissing --freq --out finalFreq
# run pre-imputation tool
/data/kronos/NGS_Software/HRC-1000G-check-bim.pl -b all_24_hhMissing.bim -f finalFreq.frq -r /data/kronos/NGS_Reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
# above perl script creates sh script for cleaning prior to imputation, run the sh script
chmod u+x Run-plink.sh
# change plink to plink2 because doesn't recognise --real-ref-alleles arg
sed -i 's/plink/\/data\/kronos\/NGS_Software\/plink2\/plink2/g' Run-plink.sh
bash Run-plink.sh
# unimputed fishers exact assoc test
mkdir unimputed_fishers
plink --bfile all_24_hhMissing-updated --fisher --model --out unimputed_fishers/unimputed_fishers
# convert plink to vcf
# bash convertToVcf
for i in `seq 1 22`
do
        plink --bfile all_24_hhMissing-updated-chr${i} --keep-allele-order --recode vcf-iid --out for_imputation_chr${i}
        bgzip for_imputation_chr${i}.vcf
        tabix -p vcf for_imputation_chr${i}.vcf.gz
done
# do X chr alone to change 23 to X first
plink --bfile all_24_hhMissing-updated-chr23 --keep-allele-order --recode vcf-iid --out for_imputation_chr23
awk -F$'\t' 'BEGIN {OFS = FS} {sub("23","X",$1)}1' for_imputation_chr23.vcf > for_imputation_chrX.vcf
bgzip for_imputation_chrX.vcf
tabix -p vcf for_imputation_chrX.vcf.gz
# imputed on Michigan imputation server as per GWAS guideline
# download data using wget links provided
# hades
bash getImputed.sh
# unzip per-chromosome files
for i in `seq 1 22`
do
        unzip -P'<qTdVPSgy7$Yq9' chr_${i}.zip &
done
unzip chr_X.zip
# check sample orders
bcftools query -l chr22.dose.vcf.gz > ID_order_autosomes.txt
bcftools query -l chrX.dose.vcf.gz > ID_order_chrX.txt
diff ID_order_autosomes.txt ID_order_chrX.txt
# get lists of monomorphic
for i in `seq 1 22`;
do
bcftools view -q1.0:major -S ID_order_chrX.txt chr${i}.dose.vcf.gz -O v | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' >> imputed_monomorphic.txt
done
# create vcf with only polymophics
# made err dir, had file truncation issues
mkdir err
for i in `seq 1 22`;
do
bcftools view -c1.0:minor -S ID_order_chrX.txt chr${i}.dose.vcf.gz -O z > chr${i}.imputed.poly.vcf.gz 2> err/chr${i}.imputed.poly.err &
done
# remove rare and low qual
for i in `seq 1 22`;
do
bcftools filter -i 'AF > 0.01 & R2 > 0.3' -O z chr${i}.imputed.poly.vcf.gz > chr${i}.imputed.filt.vcf.gz 2> err/chr${i}.imputed.filt.err &
done
# and index
for i in `seq 1 22`;
do
tabix -p vcf chr${i}.imputed.filt.vcf.gz 2> err/chr${i}.imputed.filt.tabix.err &
done
# saige prep
# create phenotype file in R using fam file + eigenvecs file
# sample file
cp ID_order_chrX.txt saige.sample
# pruned plink
awk '$4=(FNR FS $4)' highLDRegions.txt > highLDRegionsIds.txt
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_24_hhMissing-updated --make-set highLDRegionsIds.txt --write-set --out highLDR_v2
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_24_hhMissing-updated --exclude highLDR_v2.set  --make-bed --out all_25_nohighLDR
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_25_nohighLDR --indep-pairwise 50 5 0.2 --out LDprune_v2
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_25_nohighLDR --extract LDprune_v2.prune.in --make-bed --out saige_all_26
# run SAIGE
# scp psivakumar@144.82.49.64:/array/psivakumar/Emer_CH/gwas/uk_analysis/chr*.imputed.poly.vcf.gz* hades@144.82.49.41:/data01/hades/psivakumar/Emer_CH/gwas
# step 1
Rscript /hades/Software/NGS_Software/SAIGE/extdata/step1_fitNULLGLMM.R \
        --plinkFile=/hades/psivakumar/Emer_CH/uk_gwas/saige_all_26 \
        --phenoFile=/hades/psivakumar/Emer_CH/uk_gwas/saige.pheno \
        --phenoCol=Phenotype \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --skipModelFitting=FALSE \
        --outputPrefix=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step1 \
        --nThreads=2 \
        &> /hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step1.log

# step 2
Rscript /hades/Software/NGS_Software/SAIGE/extdata/step2_SPAtests.R \
        --vcfFile=/hades/psivakumar/Emer_CH/uk_gwas/chr22.imputed.filt.vcf.gz \
        --vcfFileIndex=/hades/psivakumar/Emer_CH/uk_gwas/chr22.imputed.filt.vcf.gz.tbi \
        --vcfField=DS \
        --chrom=22 \
        --minMAC=1 \
        --sampleFile=/hades/psivakumar/Emer_CH/uk_gwas/saige.sample \
        --GMMATmodelFile=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step1.rda \
        --varianceRatioFile=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step1.varianceRatio.txt \
        --SAIGEOutputFile=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step2_output.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE \
        &> /hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step2.log
# as loop
#!/bin/bash
saige_step2(){
  Rscript /hades/Software/NGS_Software/SAIGE/extdata/step2_SPAtests.R \
          --vcfFile=/hades/psivakumar/Emer_CH/uk_gwas/chr$1.imputed.filt.vcf.gz \
          --vcfFileIndex=/hades/psivakumar/Emer_CH/uk_gwas/chr$1.imputed.filt.vcf.gz.tbi \
          --vcfField=DS \
          --chrom=$1 \
          --minMAC=1 \
          --sampleFile=/hades/psivakumar/Emer_CH/uk_gwas/saige.sample \
          --GMMATmodelFile=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step1.rda \
          --varianceRatioFile=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step1.varianceRatio.txt \
          --SAIGEOutputFile=/hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step2_chr$1_output.txt \
          --numLinesOutput=2 \
          --IsOutputAFinCaseCtrl=TRUE
}

for i in `seq 1 22`
do
    saige_step2 ${i} &> /hades/psivakumar/Emer_CH/uk_gwas/SAIGE_step2_chr${i}.log &
done

# make plink files for SMR
for i in `seq 1 22`
do
  /hades/Software/NGS_Software/plink2_linux_x86_64_20190221/plink2 --vcf chr${i}.imputed.filt.vcf.gz --set-missing-var-ids '@:#\$r:\$a'  --new-id-max-allele-len 1000 --make-bed --allow-extra-chr --out chr${i}_imputed_filt &
done

ls -1 chr*.bed > bed_list_imputed_filt.txt
ls -1 chr*.bim > bim_list_imputed_filt.txt
ls -1 chr*.fam > fam_list_imputed_filt.txt
paste bed_list_imputed_filt.txt bim_list_imputed_filt.txt fam_list_imputed_filt.txt > plink_list_imputed_filt.txt
sed -i '/^chr1_/d' plink_list_imputed_filt.txt

plink --bfile chr1_imputed_filt --allow-no-sex --merge-list plink_list_imputed_filt.txt --make-bed --out all_imputed_filt
