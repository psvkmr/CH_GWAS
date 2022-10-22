RES=/array/psivakumar/Emer_CH/gwas/uk_analysis/SAIGE_results/imputed_all_results_uk_only.csv
NEW_RES=/array/psivakumar/Emer_CH/gwas/uk_analysis/SAIGE_results/imputed_all_results_uk_only_fuma
scp $RES $NEW_RES
sed -i '1s/p.value/p-value/' $NEW_RES
sed -i '1s/CHR/chromosome/' $NEW_RES
sed -i '1s/POS/position/' $NEW_RES
sed -i '1s/BETA/Beta/' $NEW_RES
sed -i '1s/Allele2/effect_allele/' $NEW_RES
sed -i '1s/Allele1/non_effect_allele/' $NEW_RES
# make sure OFS tab is real tab space
awk 'BEGIN {FS=",";OFS=""} {print $1,$2,$4,$5,$10,$11,$13}' $NEW_RES > $NEW_RES.txt
gzip $NEW_RES.txt
