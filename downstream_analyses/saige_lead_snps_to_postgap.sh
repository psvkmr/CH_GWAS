#!/bin/bash

# lead snps filtered from SAIGE results in R on kronos
RES=/hades/psivakumar/Emer_CH/uk_swedish_alt/lead_snps.csv
NEW_RES=/hades/psivakumar/Emer_CH/uk_swedish_alt/lead_snps
scp $RES $NEW_RES
sed -i '1s/p.value/p-value/' $NEW_RES
sed -i '1s/CHR/chromosome/' $NEW_RES
sed -i '1s/POS/base_pair_location/' $NEW_RES
sed -i '1s/BETA/beta/' $NEW_RES
sed -i '1s/Allele2/effect_allele/' $NEW_RES
sed -i '1s/Allele1/non_effect_allele/' $NEW_RES
sed -i '1s/\<ID\>/variant_id/' $NEW_RES
sed -i '1s/SE/standard_error/' $NEW_RES
sed -i '1s/\<AF_Allele2\>/effect_allele_frequency/' $NEW_RES
# make sure OFS tab is real tab space
awk 'BEGIN {FS=",";OFS=""} {print $1,$2,$4,$5,$7,$10,$11,$13,$20}' $NEW_RES > $NEW_RES_for_postgap.tsv
rm $NEW_RES

/hades/Software/Genetics_Software/postgap/anaconda/bin/python2.7 /hades/Software/Genetics_Software/postgap/POSTGAP.py --database_dir /hades/Software/Genetics_Software/postgap/databases_dir/databases/ --summary_stats $NEW_RES_for_postgap.tsv > $NEW_RES_postgap_res.txt
