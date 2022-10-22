FILES=/array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/*ldsc.txt
NAMES=$(ls ${FILES} | cut -d"/" -f 8 | awk -F 'ukbb_' '{print $2}' | awk -F '_sumstats_ldsc.txt' '{print $1}')
LENGTH=$(($(echo ${NAMES} | wc -w) -1))
ARR_F=(${FILES})
ARR_N=(${NAMES})

for i in `seq 0 ${LENGTH}`
do
  python2.7 /data/kronos/NGS_Software/ldsc/munge_sumstats.py --out /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/${ARR_N[${i}]} --merge-alleles /data/kronos/NGS_Software/ldsc/w_hm3.noMHC.snplist/w_hm3.noMHC.snplist --sumstats ${ARR_F[${i}]}
done

#add 1 to main trait output
python2.7 /data/kronos/NGS_Software/ldsc/munge_sumstats.py --out /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/1.ch --merge-alleles /data/kronos/NGS_Software/ldsc/w_hm3.noMHC.snplist/w_hm3.noMHC.snplist --sumstats /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/imputed_all_results_uk_ldhub.txt

SUMSTATS=$(ls /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/*.sumstats.gz | cut -d "/" -f 8)
SUMSTAT_NAMES=$(echo $SUMSTATS | tr ' ' ',')
python2.7 /data/kronos/NGS_Software/ldsc/ldsc.py --rg $SUMSTAT_NAMES --ref-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --w-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --out ch_all_ukbb
# to get heritability on liability scale, need population and sumstat prevalences of all diseases being tested
python2.7 /data/kronos/NGS_Software/ldsc/ldsc.py --rg $SUMSTAT_NAMES --ref-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --w-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --samp-prev 0.14,0.0028,0.057,0.0083,0.029,0.0011 --pop-prev 0.001,0.03,0.05,0.007,0.15,0.00025 --out ch_ukbb_all_liability_scale_h2


#python2.7 /data/kronos/NGS_Software/ldsc/ldsc.py --rg ch.sumstats.gz,anxiety.sumstats.gz,bipolar.sumstats.gz,depression.sumstats.gz,epilepsy.sumstats.gz,head.sumstats.gz,headache.sumstats.gz,migraine.sumstats.gz,nervous.sumstats.gz,schizophrenia.sumstats.gz --ref-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --w-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --out ch_ukbb_all
# to get heritability on liability scale, need population and sumstat prevalences of all diseases being tested
#python2.7 /data/kronos/NGS_Software/ldsc/ldsc.py --rg ch.sumstats.gz,bipolar.sumstats.gz,depression.sumstats.gz,epilepsy.sumstats.gz,migraine.sumstats.gz,schizophrenia.sumstats.gz --ref-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --w-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --samp-prev 0.14,0.0028,0.057,0.0083,0.029,0.0011 --pop-prev 0.001,0.03,0.05,0.007,0.15,0.00025 --out ch_ukbb_all_liability_scale_h2

#python2.7 /data/kronos/NGS_Software/ldsc/munge_sumstats.py --out /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/ch --merge-alleles /data/kronos/NGS_Software/ldsc/w_hm3.noMHC.snplist/w_hm3.noMHC.snplist --sumstats /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/imputed_all_results_uk_ldhub.txt
#python2.7 /data/kronos/NGS_Software/ldsc/munge_sumstats.py --out /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/ukbb_mig --merge-alleles /data/kronos/NGS_Software/ldsc/w_hm3.noMHC.snplist/w_hm3.noMHC.snplist --sumstats /array/psivakumar/Emer_CH/gwas/uk_analysis/ldsc/ukbb_all_sumstats_ldhub.txt
#python2.7 /data/kronos/NGS_Software/ldsc/ldsc.py --rg ch.sumstats.gz,ukbb_mig.sumstats.gz --ref-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --w-ld-chr /data/kronos/NGS_Software/ldsc/eur_w_ld_chr/ --out ch_ukbb_mig
