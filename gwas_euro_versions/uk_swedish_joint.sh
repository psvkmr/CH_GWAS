#!/bin/bash

SAMPLES=( "1958c" "GSA2016" "GSA2018" "nbs" "reseq" "swedish")
# run run_premerge_QC.sh on each of the 2016, 2018, NBS and 1958 datasets
bash run_premerge_QC.sh 1958controls_14.09.2018 1958c
bash run_premerge_QC.sh GSA2016_142_025_CC GSA2016
bash run_premerge_QC.sh GSA2018_310_025-CC GSA2018
bash run_premerge_QC.sh NBScontrols_14.09.2018 nbs
bash run_premerge_QC.sh Plate_3_Emer reseq
bash run_premerge_QC.sh Dataset_UK_newIDord swedish
# In R, intersect maps to create common snp list, called commonSnpsToExtract.txt, and extract
#plink --bfile GSA2016_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out GSA2016_4_commonSnps
#plink --bfile GSA2018_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out GSA2018_4_commonSnps
#plink --bfile 1958c_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out 1958c_4_commonSnps
#plink --bfile nbs_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out nbs_4_commonSnps
#plink --bfile reseq_3_biallelic --extract commonSnpsToExtract.txt --make-bed --recode --out reseq_4_commonSnps
# edit fam file phenotypes in R
# rename samples
awk '{print $2,"1958c:"$1":"$2,$3,$4,$5,$6}' 1958c_3_biallelic.fam > rename_1958c.txt
mv rename_1958c.txt 1958c_3_biallelic.fam
awk '{print $2,"GSA2016:"$1":"$2,$3,$4,$5,$6}' GSA2016_3_biallelic.fam > rename_GSA2016.txt
mv rename_GSA2016.txt GSA2016_3_biallelic.fam
awk '{print $2,"GSA2018:"$1":"$2,$3,$4,$5,$6}' GSA2018_3_biallelic.fam > rename_GSA2018.txt
mv rename_GSA2018.txt GSA2018_3_biallelic.fam
awk '{print $2,"nbs:"$1":"$2,$3,$4,$5,$6}' nbs_3_biallelic.fam > rename_nbs.txt
mv rename_nbs.txt nbs_3_biallelic.fam
awk '{print $2,"reseq:"$1":"$2,$3,$4,$5,$6}' reseq_3_biallelic.fam > rename_reseq.txt
mv rename_reseq.txt reseq_3_biallelic.fam
awk '{print $2,"swedish:"$1":"$2,$3,$4,$5,$6}' swedish_3_biallelic.fam > rename_swedish.txt
mv rename_swedish.txt swedish_3_biallelic.fam
# merge
#plink --bfile GSA2016_4_commonSnps --bmerge GSA2018_4_commonSnps --merge-mode 4 --out tmp1
#plink --bfile tmp1 --bmerge 1958c_4_commonSnps --merge-mode 4 --out tmp2
#plink --bfile tmp2 --bmerge nbs_4_commonSnps --merge-mode 4 --out tmp3
#plink --bfile tmp3 --bmerge reseq_4_commonSnps --merge-mode 4 --out all_5_merged
#plink --bfile all_5_merged --recode --out all_5_merged
# edit fam file to include phenotypes in R
# check for any strand errors
# remove monomorphic snps
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_3_biallelic --exclude ${SAMPLES[i]}_monomorphic_snps.txt --make-bed --out ${SAMPLES[i]}_4_noMonomorph
done
# remove ambiguous strand snps
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_4_noMonomorph --exclude ${SAMPLES[i]}_ambiguous_strand_snps.txt --make-bed --out ${SAMPLES[i]}_5_noAmbigStrand
done
# postmerge run
#for i in ${!SAMPLES[@]}
#do
#  plink --bfile ${SAMPLES[i]}_5_noAmbigStrand --allow-no-sex --flip-scan --out ${SAMPLES[i]}_postMerge
#done
# create list of strand issue variants to remove in R
#plink --bfile all_5_merged --exclude strandIssueSnpsToRemove.txt --allow-no-sex --make-bed --recode --out all_6_noStrandIssue
# split by case and control
#awk '$6==2{print $1,$2}' all_6_noStrandIssue.fam > ukCases.txt
#awk '$6==1{print $1,$2}' all_6_noStrandIssue.fam > ukControls.txt
#plink --bfile all_6_noStrandIssue --keep ukCases.txt --make-bed --recode --out cases_7_split
#plink --bfile all_6_noStrandIssue --keep ukControls.txt --make-bed --recode --out controls_7_split
# generate missingness stats
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_5_noAmbigStrand --missing --allow-no-sex --out ${SAMPLES[i]}_missingnessStats
done
# remove high missing SNPs
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_5_noAmbigStrand --geno 0.02 --allow-no-sex --make-bed --out ${SAMPLES[i]}_6_snps98
done
# remove any high missing ind
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_6_snps98 --mind 0.02 --allow-no-sex --make-bed --out ${SAMPLES[i]}_7_ind98
done
# sex check
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_7_ind98 --check-sex --out ${SAMPLES[i]}_sexCheck
done
# get list of problem sex
# F statistic dist looked fine, didn't remove any samples
# grep PROBLEM controlSex.sexcheck | awk '{print $1,$2}' > controlFailSexCheck.txt
# grep PROBLEM caseSex.sexcheck | awk '{print $1,$2}' > caseFailSexCheck.txt
# edit duplicate ids (from duplicated_ids_to_edit in R)
sed -e '0,/GSA2016:6:96843/ s/GSA2016:6:96843/GSA2016:6:96843_edit/' GSA2016_7_ind98.fam > tmp.fam
sed -e '0,/GSA2016:97:65920/ s/GSA2016:97:65920/GSA2016:97:65920_edit/' tmp.fam > tmp2.fam
sed -e '0,/GSA2016:37:68048/ s/GSA2016:37:68048/GSA2016:37:68048_edit/' tmp2.fam > tmp.fam
sed -e '0,/GSA2016:237:67695/ s/GSA2016:237:67695/GSA2016:237:67695_edit/' tmp.fam > tmp2.fam
sed -e '0,/GSA2016:217:64247/ s/GSA2016:217:64247/GSA2016:217:64247_edit/' tmp2.fam > tmp.fam
sed -e '0,/GSA2016:201:96844/ s/GSA2016:201:96844/GSA2016:201:96844_edit/' tmp.fam > tmp2.fam
sed -e '0,/GSA2016:341:96845/ s/GSA2016:341:96845/GSA2016:341:96845_edit/' tmp2.fam > GSA2016_7_ind98.fam
sed -e '0,/GSA2018:434:101172/ s/GSA2018:434:101172/GSA2018:434:101172_edit/' GSA2018_7_ind98.fam > tmp.fam
sed -e '0,/GSA2018:606:103485/ s/GSA2018:606:103485/GSA2018:606:103485_edit/' tmp.fam > GSA2018_7_ind98.fam
# remove sex problem samples
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_7_ind98 --remove ${SAMPLES[i]}_fail_sexcheck.txt --make-bed --out ${SAMPLES[i]}_8_sexChecked
done
# impute sex instead
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_7_ind98 --impute-sex --make-bed --out ${SAMPLES[i]}_9_sexImputed
done
# generate hwe stats
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_9_sexImputed --hardy --out ${SAMPLES[i]}_HWEcheck
done
# remove hwe deviants
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_9_sexImputed --hwe 1e-6 --make-bed --out ${SAMPLES[i]}_10_hwe
done
# remove het outliers
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_10_hwe --het --out ${SAMPLES[i]}_het
done
# In R, find samples with outlier number of hets
# Remove
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_10_hwe --remove ${SAMPLES[i]}_fail_hetcheck.txt --make-bed --out ${SAMPLES[i]}_11_noHetOutliers
done
# prune for ld
# got set of high ld regions from https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)
awk '$4=(FNR FS $4)' highLDRegions.txt > highLDRegionsIds.txt
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_11_noHetOutliers --make-set highLDRegionsIds.txt --write-set --out ${SAMPLES[i]}_highLDR
done
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_11_noHetOutliers --exclude ${SAMPLES[i]}_highLDR.set --make-bed --out ${SAMPLES[i]}_12_noHighLd
done
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_12_noHighLd --indep-pairwise 50 5 0.2 --out ${SAMPLES[i]}_LDprune
done
# prune and calculate relatedness
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_12_noHighLd --extract ${SAMPLES[i]}_LDprune.prune.in --genome --out ${SAMPLES[i]}_relCheck
done
# remove duplicates/monozygotic twins
for i in ${!SAMPLES[@]}
do
  awk '$10>0.2 {print $2,$4,$10}' ${SAMPLES[i]}_relCheck.genome > ${SAMPLES[i]}_relatednessInfo.txt
  awk '$10>0.2 {print $1,$2}' ${SAMPLES[i]}_relCheck.genome | sed '1d' > ${SAMPLES[i]}_relatedSamples.txt
  plink --bfile ${SAMPLES[i]}_11_noHetOutliers --remove ${SAMPLES[i]}_relatedSamples.txt --make-bed --out ${SAMPLES[i]}_13_noRelated
done
# use peddy for pop strat
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_13_noRelated --keep-allele-order --recode vcf-iid --out ${SAMPLES[i]}_forPeddy
  plink --bfile ${SAMPLES[i]}_13_noRelated --recode --out ${SAMPLES[i]}_13_noRelated
  bgzip ${SAMPLES[i]}_forPeddy.vcf
  tabix -p vcf ${SAMPLES[i]}_forPeddy.vcf.gz
  python3.6 -m peddy -p 4 --plot --prefix ${SAMPLES[i]}_peddy_gwas ${SAMPLES[i]}_forPeddy.vcf.gz ${SAMPLES[i]}_13_noRelated.ped
done
# create list of non-european to exclude in R
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_13_noRelated --keep ${SAMPLES[i]}_european_to_keep.txt --make-bed --out ${SAMPLES[i]}_14_european
done
# pca
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_14_european --pca --out ${SAMPLES[i]}_pca
done
# set hh to missing
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_14_european --set-hh-missing --make-bed --out ${SAMPLES[i]}_15_hhMissing
done
# get final freq file
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_15_hhMissing --freq --out ${SAMPLES[i]}_finalFreq
done
# run pre-imputation tool
for i in ${!SAMPLES[@]}
do
  /data/kronos/NGS_Software/HRC-1000G-check-bim.pl -b ${SAMPLES[i]}_15_hhMissing.bim -f ${SAMPLES[i]}_finalFreq.frq -r /data/kronos/NGS_Reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
  chmod u+x Run-plink.sh
  sed -i 's/plink/\/data\/kronos\/NGS_Software\/plink2\/plink2/g' Run-plink.sh
  mv Run-plink.sh ${SAMPLES[i]}_Run-plink.sh
done
# run plink scripts
bash
for i in ${!SAMPLES[@]}
do
  bash "${SAMPLES[i]}_Run-plink.sh"
done
# unimputed fishers exact assoc test
#mkdir unimputed_fishers
#for i in ${!SAMPLES[@]}
#do
#  plink --bfile ${SAMPLES[i]}_15_hhMissing-updated --fisher --model --out unimputed_fishers/${SAMPLES[i]}_unimputed_fishers
#done
# convert plink to vcf
# bash convertToVcf
for i in ${!SAMPLES[@]}
do
  for j in `seq 1 22`
  do
          plink --bfile ${SAMPLES[i]}_15_hhMissing-updated-chr${j} --keep-allele-order --recode vcf-iid --out ${SAMPLES[i]}_for_imputation_chr${j}
          bgzip ${SAMPLES[i]}_for_imputation_chr${j}.vcf
          tabix -p vcf ${SAMPLES[i]}_for_imputation_chr${j}.vcf.gz
  done
done
# do X chr alone to change 23 to X first
for i in ${!SAMPLES[@]}
do
  plink --bfile ${SAMPLES[i]}_15_hhMissing-updated-chr23 --keep-allele-order --recode vcf-iid --out ${SAMPLES[i]}_for_imputation_chr23
  awk -F$'\t' 'BEGIN {OFS = FS} {sub("23","X",$1)}1' ${SAMPLES[i]}_for_imputation_chr23.vcf > ${SAMPLES[i]}_for_imputation_chrX.vcf
  bgzip ${SAMPLES[i]}_for_imputation_chrX.vcf
  tabix -p vcf ${SAMPLES[i]}_for_imputation_chrX.vcf.gz
done
# imputed on Michigan imputation server as per GWAS guideline
# download data using wget links provided
# hades
for i in ${!SAMPLES[@]}
do
  bash "${SAMPLES[i]}_Run-plink.sh"
done
chmod u+x getImputed.sh
sed -i 's/zip$/zip \&/g' getImputed.sh
bash getImputed.sh
# unzip per-chromosome files
for i in `seq 1 22`
do
        unzip -P'o;KnJ8cGX9as$K' chr_${i}.zip &
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
        --plinkFile=/hades/psivakumar/Emer_CH/uk_only_gwas/saige_all_26 \
        --phenoFile=/hades/psivakumar/Emer_CH/uk_only_gwas/saige.pheno \
        --phenoCol=Phenotype \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --skipModelFitting=FALSE \
        --outputPrefix=/hades/psivakumar/Emer_CH/uk_only_gwas/SAIGE_step1 \
        --nThreads=2 \
        &> /hades/psivakumar/Emer_CH/uk_only_gwas/SAIGE_step1.log

# step 2
# as loop
#!/bin/bash
saige_step2(){
  Rscript /hades/Software/NGS_Software/SAIGE/extdata/step2_SPAtests.R \
          --vcfFile=/hades/psivakumar/Emer_CH/uk_only_gwas/chr$1.imputed.filt.vcf.gz \
          --vcfFileIndex=/hades/psivakumar/Emer_CH/uk_only_gwas/chr$1.imputed.filt.vcf.gz.tbi \
          --vcfField=DS \
          --chrom=$1 \
          --minMAC=1 \
          --sampleFile=/hades/psivakumar/Emer_CH/uk_only_gwas/saige.sample \
          --GMMATmodelFile=/hades/psivakumar/Emer_CH/uk_only_gwas/SAIGE_step1.rda \
          --varianceRatioFile=/hades/psivakumar/Emer_CH/uk_only_gwas/SAIGE_step1.varianceRatio.txt \
          --SAIGEOutputFile=/hades/psivakumar/Emer_CH/uk_only_gwas/SAIGE_step2_chr$1_output.txt \
          --numLinesOutput=2 \
          --IsOutputAFinCaseCtrl=TRUE
}

for i in `seq 1 22`
do
    saige_step2 ${i} &> /hades/psivakumar/Emer_CH/uk_only_gwas/SAIGE_step2_chr${i}.log &
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
