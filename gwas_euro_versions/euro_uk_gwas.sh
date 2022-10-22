#!/bin/bash

# check samples in sample sheets and fam files

# remove samples with suspicious greek ids
#plink --bfile ICH_controls --remove suspected_greek_ids.txt --make-bed --out ICH_trimmed

# rename ICH controls to avoid duplicated ids
#awk '{print $1,$2,"control_"$2,"control_"$2}' GSA2016_142_025_CC_04-10.fam > update_ICH_ids.txt
#plink --bfile GSA2016_142_025_CC_04-10 --update-ids update_ICH_ids.txt --make-bed --out ICH_renamed

# remove gsa2019 duplicate ids, checked in R
#awk '!seen[$2]++' GSA2019_432_025_CC.ped > GSA2019_trimmed.ped
#awk '!seen[$0]++' GSA2019_432_025_CC.nosex > GSA2019_trimmed.nosex
#scp GSA2019_432_025_CC.map GSA2019_trimmed.map
#plink --file GSA2019_trimmed --make-bed --out GSA2019_trimmed

# make gsa2019 fids same as iids
#awk '{print $1,$2,$2,$2}' GSA2019_trimmed.fam > update_GSA2019_ids.txt
#plink --bfile GSA2019_trimmed --update-ids update_GSA2019_ids.txt --make-bed --out GSA2019_renamed

# remove unaccounted for gsa2019 samples
#plink --bfile GSA2019_renamed --remove unknown_gsa2019_ids_remove.txt --make-bed --out GSA2019_known

# change var ids
bash run_premerge_QC.sh GSA2019_432_025_CC gsa2019
bash run_premerge_QC.sh GSA2016_142_025_CC_04-10 ich
bash run_premerge_QC.sh GSA2016_142_025_CC GSA2016
bash run_premerge_QC.sh GSA2018_310_025-CC GSA2018
bash run_premerge_QC.sh 1958controls_14.09.2018 1958c
bash run_premerge_QC.sh NBScontrols_14.09.2018 nbs
bash run_premerge_QC.sh Plate_3_Emer reseq
# exlcude monomorphic and ambiguous strand variants, checked in R
plink --bfile gsa2019_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out gsa2019_4_commonSnps
plink --bfile ich_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out ich_4_commonSnps
plink --bfile GSA2016_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out GSA2016_4_commonSnps
plink --bfile GSA2018_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out GSA2018_4_commonSnps
plink --bfile 1958c_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out 1958c_4_commonSnps
plink --bfile nbs_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out nbs_4_commonSnps
plink --bfile reseq_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out reseq_4_commonSnps
# merge
plink --bfile GSA2016_4_commonSnps --bmerge GSA2018_4_commonSnps --merge-mode 4 --out tmp1
plink --bfile tmp1 --bmerge 1958c_4_commonSnps --merge-mode 4 --out tmp2
plink --bfile tmp2 --bmerge nbs_4_commonSnps --merge-mode 4 --out tmp3
plink --bfile tmp3 --bmerge reseq_4_commonSnps --merge-mode 4 --out tmp4
plink --bfile tmp4 --bmerge gsa2019_4_commonSnps --merge-mode 4 --out tmp5
plink --bfile tmp5 --bmerge ich_4_commonSnps --merge-mode 4 --out all_5_merged
plink --bfile all_5_merged --recode --out all_5_merged
# edit fam file to include phenotypes in R
# check for any strand errors
plink --bfile all_5_merged --allow-no-sex --flip-scan --out postMerge
# create list of strand issue variants to remove in R
plink --bfile all_5_merged --exclude strandIssueSnpsToRemove.txt --allow-no-sex --make-bed --recode --out all_6_noStrandIssue
# split by case and control
awk '$6==2{print $1,$2}' all_6_noStrandIssue.fam > allCases.txt
awk '$6==1{print $1,$2}' all_6_noStrandIssue.fam > allControls.txt
plink --bfile all_6_noStrandIssue --keep allCases.txt --make-bed --recode --out cases_7_split
plink --bfile all_6_noStrandIssue --keep allControls.txt --make-bed --recode --out controls_7_split
# add case and control prefix to sample names
awk '{print $1,$2,$1,$1"_"$2}' cases_7_split.fam > rename_cases.txt
plink --bfile cases_7_split --update-ids rename_cases.txt --make-bed --out cases_7_split
awk '{print $1,$2,$1,$1"_"$2}' controls_7_split.fam > rename_controls.txt
plink --bfile controls_7_split --update-ids rename_controls.txt --make-bed --out controls_7_split
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
plink --bfile controls_9_ind98 --remove controlFailSexCheck.txt --make-bed --out controls_10_sexCleaned
plink --bfile cases_9_ind98 --remove caseFailSexCheck.txt --make-bed --out cases_10_sexCleaned
# impute sex instead
plink --bfile controls_10_sexCleaned --impute-sex --make-bed --out controls_10_sexChecked
plink --bfile cases_10_sexCleaned --impute-sex --make-bed --out cases_10_sexChecked
# set higher missing snp threshold
plink --bfile controls_10_sexChecked --geno 0.02 --make-bed --out controls_11_snps98
plink --bfile cases_10_sexChecked --geno 0.02 --make-bed --out cases_11_snps98
# generate hwe stats
plink --bfile controls_11_snps98 --hardy --out controlsHWEcheck
plink --bfile cases_11_snps98 --hardy --out casesHWEcheck
# remove hwe deviants
plink --bfile controls_11_snps98 --hwe 1e-6/1e-10 --make-bed --out controls_12_hwe
plink --bfile cases_11_snps98 --hwe 1e-6/1e-10 --make-bed --out cases_12_hwe
#generate allele freq stats
plink --bfile controls_11_snps98 --freq --out controlRecalcFreq
plink --bfile cases_11_snps98 --freq --out casesRecalcFreq
# re-merge
plink --bfile controls_12_hwe --bmerge cases_12_hwe --allow-no-sex --make-bed --out all_13_remerge
# change to unique fids
awk '$1=(FNR FS $1)' all_13_remerge.fam > tmp.fam
awk '{print $2,$3,$1,$3}' tmp.fam > new_fids.txt
plink --bfile all_13_remerge --update-ids new_fids.txt --make-bed --out all_13_remerge
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
chmod u+x getImputed.sh
sed -i 's/zip$/zip \&/g' getImputed.sh
bash getImputed.sh
# unzip per-chromosome files
for i in `seq 1 22`
do
        unzip -P-omF4JI77yTiLm chr_${i}.zip &
done
unzip chr_X.zip
# check sample orders
bcftools query -l chr22.dose.vcf.gz > ID_order_autosomes.txt
bcftools query -l chrX.dose.vcf.gz > ID_order_chrX.txt
diff ID_order_autosomes.txt ID_order_chrX.txt
# create full vcf file
ls chr[0-9]*.dose.vcf.gz | sort -n -t r -k 2 > imputed_files.txt
bcftools concat -f imputed_files.txt -O z -o all.dose.vcf.gz
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
awk '$4=(FNR FS $4)' /hades/psivakumar/highLDRegions.txt > highLDRegionsIds.txt
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_24_hhMissing-updated --make-set highLDRegionsIds.txt --write-set --out highLDR_v2
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_24_hhMissing-updated --exclude highLDR_v2.set  --make-bed --out all_25_nohighLDR
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_25_nohighLDR --indep-pairwise 50 5 0.2 --out LDprune_v2
/hades/Software/NGS_Software/plink_linux_x86_64_dev/plink --bfile all_25_nohighLDR --extract LDprune_v2.prune.in --make-bed --out saige_all_26
# run SAIGE
# hades
#scp psivakumar@144.82.49.64:/array/psivakumar/Emer_CH/euro_uk_gwas/uk_analysis/chr*.imputed.filt.vcf.gz* hades@144.82.49.41:/data01/hades/psivakumar/Emer_CH/euro_uk_gwas
# step 1
Rscript /hades/Software/NGS_Software/SAIGE/extdata/step1_fitNULLGLMM.R \
        --plinkFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/saige_all_26 \
        --phenoFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/saige.pheno \
        --phenoCol=Phenotype \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --skipModelFitting=FALSE \
        --outputPrefix=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step1 \
        --nThreads=2 \
        &> /hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step1.log

# step 2
Rscript /hades/Software/NGS_Software/SAIGE/extdata/step2_SPAtests.R \
        --vcfFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/chr22.imputed.filt.vcf.gz \
        --vcfFileIndex=/hades/psivakumar/Emer_CH/euro_uk_gwas/chr22.imputed.filt.vcf.gz.tbi \
        --vcfField=DS \
        --chrom=22 \
        --minMAC=1 \
        --sampleFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/saige.sample \
        --GMMATmodelFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step1.rda \
        --varianceRatioFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step1.varianceRatio.txt \
        --SAIGEOutputFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step2_output.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE \
        &> /hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step2.log
# as loop
#!/bin/bash
saige_step2(){
  Rscript /hades/Software/NGS_Software/SAIGE/extdata/step2_SPAtests.R \
          --vcfFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/chr$1.imputed.filt.vcf.gz \
          --vcfFileIndex=/hades/psivakumar/Emer_CH/euro_uk_gwas/chr$1.imputed.filt.vcf.gz.tbi \
          --vcfField=DS \
          --chrom=$1 \
          --minMAC=1 \
          --sampleFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/saige.sample \
          --GMMATmodelFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step1.rda \
          --varianceRatioFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step1.varianceRatio.txt \
          --SAIGEOutputFile=/hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step2_chr$1_output.txt \
          --numLinesOutput=2 \
          --IsOutputAFinCaseCtrl=TRUE
}

for i in `seq 1 22`
do
    saige_step2 ${i} &> /hades/psivakumar/Emer_CH/euro_uk_gwas/SAIGE_step2_chr${i}.log &
done

#get sig hits from all vcf
# sig hits file made from R
for i in `seq 1 22`;
do
tabix -p vcf chr${i}.dose.vcf.gz 2> err/chr${i}.dose.tabix.err &
done
bcftools concat -a -f imputed_files.txt -R imputed_sig_results_vcf_format.txt -O z -o all.sig.dose.vcf.gz
