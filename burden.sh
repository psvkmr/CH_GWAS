#!/bin/bash

# /array/psivakumar/Emer_CH/burden
mkdir err
# split vcf records which contain multiple alt alleles into separate records
# remove uninformative alt alleles listed as *, as plink cannot handle these
bcftools norm --remove-duplicates --multiallelics - --output-type u ../samples_for_emer.vep.out.vcf.gz  | \
bcftools view -i'ALT!="*"' -O u | \
bcftools view -c1.0:minor -O z -o tmp1_no_multiallelics.vcf.gz 2> err/tmp1_no_multiallelics.err
#QC
# for pseq, need to edit vcf header because pseq cannot handle vcfv4.2, only 4.1:
bgzip -cd tmp1_no_multiallelics.vcf.gz | sed 's/VCFv4\.2/VCFv4.1/' | bgzip -c > tmp2_for_pseq.vcf.gz
tabix -p vcf tmp2_for_pseq.vcf.gz
# get individual stats
/data/kronos/NGS_Software/plinkseq-0.10/pseq tmp2_for_pseq.vcf.gz i-stats > pre_QC_istats.txt
#get variant stats
/data/kronos/NGS_Software/plinkseq-0.10/pseq tmp2_for_pseq.vcf.gz v-stats > pre_QC_vstats.txt
# create file of istats samples to remove in R
bcftools stats -s - tmp1_no_multiallelics.vcf.gz > pre_QC_bcftools_stats.txt

# filter by sample read depth
# get sample read depth data from wgs_metrics file
#bcftools view -S samples_more_15x_keep.txt -O z -o tmp3_sample_depth_filt.vcf.gz tmp1_no_multiallelics.vcf.gz

# convert vcf to plink, allow extra chromosomes for unplaced/unlocalised scaffolds, use plink2 to set unique ids for variants
/data/kronos/NGS_Software/plink2/plink2 --vcf tmp1_no_multiallelics.vcf.gz --set-missing-var-ids '@:#\$r:\$a'  --new-id-max-allele-len 1000 --make-bed --keep-allele-order --allow-extra-chr --out emer_ch_1
# from R, update fam file with phenotype status from Emer
# added phenotypes to fam file in R

# remove samples which are +-3SD from mean in istats metrics
plink --bfile emer_ch_1 --remove istats_outliers_to_remove.txt --allow-extra-chr --allow-no-sex --make-bed --out emer_ch_2_no_outliers

# filter by variants not found in hapmap
# need to find a good way of comparing all variants to hapmap variants

# filter in plink for only biallelic SNVs (why??), call rate, LD
plink --bfile emer_ch_2_no_outliers --geno 0.02 --allow-extra-chr --allow-no-sex --make-bed --out emer_ch_3_snps98

plink --bfile emer_ch_3_snps98 --indep-pairwise 50 5 0.2 --allow-extra-chr --allow-no-sex --out ld_snps
plink --bfile emer_ch_3_snps98 --extract ld_snps.prune.in --allow-extra-chr --allow-no-sex --make-bed --out emer_ch_4_no_ld

# exclude samples with high missingness
plink --bfile emer_ch_4_no_ld --missing --allow-extra-chr --allow-no-sex --out individual_missing_rate
plink --bfile emer_ch_4_no_ld --mind 0.05 --allow-extra-chr --allow-no-sex --make-bed --out emer_ch_5_no_missing_ind

# check inbreeding het
plink --bfile emer_ch_5_no_missing_ind --het --allow-extra-chr --allow-no-sex --out het_statistics
# create hets to remove file in R
plink --bfile emer_ch_5_no_missing_ind --remove hets_to_remove.txt --allow-extra-chr --allow-no-sex --make-bed --out emer_ch_6_no_outlier_het

# remove hwe issue variants
#plink --bfile emer_ch_6_no_outlier_het --hardy --allow-extra-chr --allow-no-sex --out hwe_check
#plink --bfile emer_ch_6_no_outlier_het --hwe 0.000001 --allow-extra-chr --allow-no-sex --make-bed --out emer_ch_7_hwe_filt

# remove related samples
plink --bfile emer_ch_6_no_outlier_het --genome --allow-no-sex --allow-extra-chr --out relatives_check
#Manually analysed and created related_to_remove.txt, need to find good way of filtering
plink --bfile emer_ch_6_no_outlier_het --remove related_to_remove.txt --allow-no-sex --allow-extra-chr --make-bed --out emer_ch_7_unrelated

# sex check, impute sex
# Impute sex
plink --bfile emer_ch_7_unrelated --impute-sex --allow-no-sex --allow-extra-chr --make-bed --recode --out emer_ch_8_sex_checked
# check for dodgy sex
plink --bfile emer_ch_8_sex_checked --check-sex 0.25 0.75 --allow-no-sex --allow-extra-chr --out sex_check

# back to vcf
plink --bfile emer_ch_8_sex_checked --allow-no-sex --allow-extra-chr --recode vcf-iid --output-chr chrMT --out tmp4_post_plink

# population stratification
# try using peddy
bgzip -c tmp4_post_plink.vcf > tmp4_post_plink.vcf.gz
tabix -p vcf tmp4_post_plink.vcf.gz
python3.6 -m peddy -p 4 --sites hg38 --plot --prefix peddy_emer tmp4_post_plink.vcf.gz emer_ch_8_sex_checked.ped

# remove samples of non-European ethnicity
# get list of non-European samples in R
bcftools view -S european_to_keep.txt -O z -o tmp5_european.vcf.gz tmp4_post_plink.vcf.gz

# back to plink for pca
/data/kronos/NGS_Software/plink2/plink2 --vcf tmp5_european.vcf.gz --set-missing-var-ids '@:#\$r:\$a' --new-id-max-allele-len 1000 --make-bed --keep-allele-order --allow-extra-chr --out emer_ch_9_european
/data/kronos/NGS_Software/plink2/plink2 --bfile emer_ch_9_european --pca 20 --allow-extra-chr --keep-allele-order --out pca
# remove outliers from this too?
plink --bfile emer_ch_9_european --remove pca_samples_to_remove.txt --allow-extra-chr --make-bed --out emer_ch_10_pca_cleaned

# back to vcf
plink --bfile emer_ch_10_pca_cleaned --allow-extra-chr --recode vcf-iid --output-chr chrMT --out tmp6_post_pca

bcftools query -l tmp6_post_pca.vcf > samples_final_set.txt
bcftools view -S samples_final_set.txt -O z -o tmp7_final_samples.vcf.gz tmp1_no_multiallelics.vcf.gz

# split vcfs into SNPs and indels
# SNPs
/data/kronos/NGS_Software/bcftools-1.9/bcftools filter -i 'TYPE="snp" & QD > 2 & FS < 60 & MQ > 40 & MQRankSum > -12.5 & ReadPosRankSum > -8 & INFO/DP > 10 & GQ > 20 & F_MISSING < 0.05 & ExcessHet < 20 & InbreedingCoeff > -0.8' -O z -o tmp8_snps.vcf.gz tmp7_final_samples.vcf.gz
# indels
/data/kronos/NGS_Software/bcftools-1.9/bcftools filter -i 'TYPE="indel" & QD > 2 & FS < 200 & ReadPosRankSum > -20 & INFO/DP > 10 & GQ > 20 & InbreedingCoeff > -0.8' -O z -o tmp8_indels.vcf.gz tmp7_final_samples.vcf.gz
tabix -p vcf tmp8_snps.vcf.gz
tabix -p vcf tmp8_indels.vcf.gz

# re-merge with bcftools concat -
bcftools concat -a -O z -o tmp9_merge_filt.vcf.gz tmp8_snps.vcf.gz tmp8_indels.vcf.gz

# alternate method of filtering by variant type
# vcftools --gzvcf tmp6_final_samples.vcf.gz --remove-indels --recode --recode-INFO-all --out test_no_indels.vcf
# vcftools --gzvcf tmp6_final_samples.vcf.gz --keep-only-indels --recode --recode-INFO-all --out test_only_indels.vcf

# set missing genotypes to individual genotypes of low quality
vcftools --gzvcf tmp9_merge_filt.vcf.gz --minGQ 20 --minDP 10 --out tmp10_set_missing --recode --recode-INFO-all

# split vcf by case/control
# use the previous samples_final_set.txt in R and split by case/control
bcftools view -S case_samples_final_set.txt -O z -o tmp11_case_only.vcf.gz tmp10_set_missing.recode.vcf
bcftools view -S control_samples_final_set.txt -O z -o tmp11_control_only.vcf.gz tmp10_set_missing.recode.vcf

# recalculate mean read depths with vcftools --site-depth
vcftools --gzvcf tmp11_case_only.vcf.gz --site-depth --out case_new_depth
vcftools --gzvcf tmp11_control_only.vcf.gz --site-depth --out control_new_depth

# create list of variant locations to be included from depth files
# include from tmp10 vcf any variants in either output file from site depth filtering
bcftools view -T read_depth_pass_to_include.txt -O z -o tmp12_var_filt.vcf.gz tmp10_set_missing.recode.vcf
tabix -p vcf tmp12_var_filt.vcf.gz

# try bcftools stats
bcftools stats -s - tmp12_var_filt.vcf.gz > post_QC_bcftools_stats.txt

# get stats for case vs control
bcftools view -S case_samples_final_set.txt -O z -o tmp13_case_only.vcf.gz tmp12_var_filt.vcf.gz
bcftools view -S control_samples_final_set.txt -O z -o tmp13_control_only.vcf.gz tmp12_var_filt.vcf.gz
tabix -p vcf tmp13_case_only.vcf.gz
tabix -p vcf tmp13_control_only.vcf.gz
bcftools stats -s - tmp13_case_only.vcf.gz > post_QC_case_bcftools_stats.txt
bcftools stats -s - tmp13_control_only.vcf.gz > post_QC_control_bcftools_stats.txt

# pseq
gzip -cd tmp12_var_filt.vcf.gz | sed 's/VCFv4\.2/VCFv4.1/' | bgzip -c > tmp13_for_post_pseq.vcf.gz
tabix -p vcf tmp13_for_post_pseq.vcf.gz
# get individual stats
/data/kronos/NGS_Software/plinkseq-0.10/pseq tmp13_for_post_pseq.vcf.gz i-stats > post_QC_istats.txt
#get variant stats
/data/kronos/NGS_Software/plinkseq-0.10/pseq tmp13_for_post_pseq.vcf.gz v-stats > post_QC_vstats.txt

# filter for MAF > 0.005 for single variant testing
bcftools filter -i 'AF > 0.005' -O v -o tmp14_MAF_greater_0-005.vcf tmp12_var_filt.vcf.gz
# filter for MAF < 0.005 for burden rare variant testing
bcftools filter -i 'AF < 0.005' -O v -o tmp14_MAF_less_0-005.vcf tmp12_var_filt.vcf.gz

# split by syn/non-syn
#grep -f protein_truncating_variant_so_terms.txt tmp14_MAF_less_0-005.vcf > tmp15_rare_protein_truncating_no_head.vcf
#grep -f missense_variant_so_terms.txt tmp14_MAF_less_0-005.vcf > tmp15_rare_missense_no_head.vcf
#grep -f synonymous_variant_so_terms.txt tmp14_MAF_less_0-005.vcf > tmp15_rare_synonymous_no_head.vcf
cat protein_truncating_variant_so_terms.txt missense_variant_so_terms.txt > non_synonymous_variant_so_terms.txt
tee >(grep -F -f non_synonymous_variant_so_terms.txt    > tmp15_rare_nonsyn_noheader.vcf) < tmp14_MAF_less_0-005.vcf | \
      grep -F -f non_synonymous_variant_so_terms.txt -v > tmp15_rare_syn.vcf
#grep -f non_synonymous_variant_so_terms.txt tmp14_MAF_less_0-005.vcf > tmp15_rare_nonsyn_no_head.vcf
bcftools view -h tmp14_MAF_less_0-005.vcf > tmp15_head_replace.vcf
#cat tmp15_head_replace.vcf tmp15_rare_protein_truncating_no_head.vcf > tmp15_rare_protein_truncating.vcf
#cat tmp15_head_replace.vcf tmp15_rare_missense_no_head.vcf > tmp15_rare_missense.vcf
#cat tmp15_head_replace.vcf tmp15_rare_synonymous_no_head.vcf > tmp15_rare_synonymous.vcf
cat tmp15_head_replace.vcf tmp15_rare_nonsyn_noheader.vcf > tmp15_rare_nonsyn.vcf
#bgzip -c tmp15_rare_protein_truncating.vcf > tmp15_rare_protein_truncating.vcf.gz
#bgzip -c tmp15_rare_missense.vcf > tmp15_rare_missense.vcf.gz
bgzip -c tmp15_rare_syn.vcf > tmp15_rare_syn.vcf.gz
bgzip -c tmp15_rare_nonsyn.vcf > tmp15_rare_nonsyn.vcf.gz
#tabix -p vcf tmp15_rare_protein_truncating.vcf.gz
#tabix -p vcf tmp15_rare_missense.vcf.gz
tabix -p vcf tmp15_rare_syn.vcf.gz
tabix -p vcf tmp15_rare_nonsyn.vcf.gz
# remove remaining monomorphics
fill-an-ac tmp15_rare_nonsyn.vcf > tmp16_rare_nonsyn.vcf
fill-an-ac tmp15_rare_syn.vcf > tmp16_rare_syn.vcf
bcftools view -c1.0:minor -O v -o rare_nonsyn.vcf tmp16_rare_nonsyn.vcf
bcftools view -c1.0:minor -O v -o rare_syn.vcf tmp16_rare_syn.vcf

mkdir rvtests
cd rvtests

sed '1s/^/fid\tiid\tfatid\tmatid\tsex\tpheno\n/' ../emer_ch_1.fam > rvtests.fam

#/data/kronos/Genetics_Software/rvtests-2.1.0/executable/rvtest \
#  --inVcf /array/psivakumar/Emer_CH/burden/tmp15_rare_syn.vcf.gz \
#  --pheno /array/psivakumar/Emer_CH/burden/rvtests/rvtests.fam \
#  --pheno-name pheno \
#  --out /array/psivakumar/Emer_CH/burden/rvtests/rvtest_rare_synonymous_variants \
#  --geneFile /data/kronos/NGS_Reference/refFlat/refFlat_hg38.txt.gz \
#  --burden exactCMC

/data/kronos/Genetics_Software/rvtests-2.1.0/executable/rvtest \
  --inVcf /array/psivakumar/Emer_CH/burden/tmp15_rare_syn.vcf.gz \
  --pheno /array/psivakumar/Emer_CH/burden/rvtests/rvtests.fam \
  --pheno-name pheno \
  --out /array/psivakumar/Emer_CH/burden/rvtests/rvtest_rare_syn_variants \
  --setFile /array/psivakumar/Emer_CH/burden/rvtests/hg38_genes.set \
  --burden exactCMC

/data/kronos/Genetics_Software/rvtests-2.1.0/executable/rvtest \
  --inVcf /array/psivakumar/Emer_CH/burden/tmp15_rare_nonsyn.vcf.gz \
  --pheno /array/psivakumar/Emer_CH/burden/rvtests/rvtests.fam \
  --pheno-name pheno \
  --out /array/psivakumar/Emer_CH/burden/rvtests/rvtest_rare_nonsyn_variants \
  --setFile /array/psivakumar/Emer_CH/burden/rvtests/hg38_genes.set \
  --burden exactCMC

# Bonferroni correction in R
