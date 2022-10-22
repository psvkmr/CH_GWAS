RES=/array/psivakumar/Emer_CH/gwas/uk_analysis/SAIGE_results/imputed_all_results_uk_only.csv
NEW_RES=/array/psivakumar/Emer_CH/gwas/uk_analysis/SAIGE_results/imputed_all_results_uk_only_rsIDs
awk -F"," '{print $3}' $RES > $NEW_RES
sed -i '1d' $NEW_RES
awk -F":" 'BEGIN {OFS="       ";} {print $1,$2,$3,$4}' $NEW_RES > $NEW_RES.tmp
sed -i 's/       /       \.      /2' $NEW_RES.tmp
/data/kronos/NGS_Software/VEP_95/ensembl-vep/variant_recoder --input_file $NEW_RES.tmp --grch37 --pretty > $NEW_RES.json

# test run
head -n 10 ${NEW_RES}.tmp > ${NEW_RES}_test.tmp
/data/kronos/NGS_Software/VEP_95/ensembl-vep/variant_recoder --input_file ${NEW_RES}_test.tmp --grch37 --pretty > ${NEW_RES}_test.json

# parallelise
mkdir imputed_res_by_chr
awk -F"       " '{print>"imputed_res_by_chr/"$1}' $NEW_RES.tmp
cd imputed_res_by_chr
for i in `seq 1 22`;
do
  /data/kronos/NGS_Software/VEP_95/ensembl-vep/variant_recoder --input_file ${i} --grch37 --pretty > ${i}.json 2> ${i}.log &
done

for i in `seq 3 4`;
do
  /data/kronos/NGS_Software/VEP_95/ensembl-vep/variant_recoder --input_file ${i} --grch37 --pretty > ${i}.json 2> ${i}.log &
done

/data/kronos/NGS_Software/VEP_95/ensembl-vep/variant_recoder --input_file 5 --grch37 --pretty > 5.json 2> 5.log &

# try curl
# vcf to hgvsg format
awk -F" " '{print $1":g."$2$4">"$5}' ${NEW_RES}_test.tmp > ${NEW_RES}_test_hgvsg

curl 'http://grch37.rest.ensembl.org/variant_recoder/homo_sapiens' -H 'Content-type:application/json' -H 'Accept:application/json' -X POST -d '{ "ids" : ["1:g.833223C>T", "1:g.833302C>T", "1:g.833641T>C", "1:g.833824T>C", "1:g.833927T>C", "1:g.834198T>C", "1:g.834832G>C", "1:g.834928A>G", "1:g.834999G>A", "1:g.835499A>G" ] }'

# never finished, found R method to get rsids in saige_to_ldsc.R
