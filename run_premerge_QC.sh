#!/bin/bash
# $1 is input, $2 is output prefix
/data/kronos/NGS_Software/plink2/plink2 --bfile $1 --set-all-var-ids '@:#$r:$a' --new-id-max-allele-len 1000 --make-bed --out $2_1_newIds
# need to understand what to do with inconsistent hh calls, so leaving X chr out for now
# plink --bfile $2_1_newIds --chr 1-22 --make-bed --out $2_2_chr1_22
/data/kronos/NGS_Software/plink2/plink2 --bfile $2_1_newIds --rm-dup force-first --make-bed --out $2_2_rmDupIds
plink --bfile $2_2_rmDupIds --biallelic-only strict --make-bed --recode --out $2_3_biallelic
plink --bfile $2_3_biallelic --freq --out alleleFreqs_$2
