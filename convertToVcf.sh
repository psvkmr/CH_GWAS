#!/bin/bash
for i in `seq 1 22`
do 
	plink --bfile all_24_hhMissing-updated-chr${i} --keep-allele-order --recode vcf-iid --out for_imputation_chr${i}
	bgzip for_imputation_chr${i}.vcf
	tabix -p vcf for_imputation_chr${i}.vcf.gz
done
