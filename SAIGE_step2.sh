#!/bin/bash
saige_step2(){
  Rscript /hades/Software/NGS_Software/SAIGE/extdata/step2_SPAtests.R \
          --vcfFile=/hades/psivakumar/Emer_CH/gwas/chr$1.imputed.poly.vcf.gz \
          --vcfFileIndex=/hades/psivakumar/Emer_CH/gwas/chr$1.imputed.poly.vcf.gz.tbi \
          --vcfField=DS \
          --chrom=$1 \
          --minMAC=1 \
          --sampleFile=/hades/psivakumar/Emer_CH/gwas/saige.sample \
          --GMMATmodelFile=/hades/psivakumar/Emer_CH/gwas/SAIGE_step1.rda \
          --varianceRatioFile=/hades/psivakumar/Emer_CH/gwas/SAIGE_step1.varianceRatio.txt \
          --SAIGEOutputFile=/hades/psivakumar/Emer_CH/gwas/SAIGE_step2_chr$1_output.txt \
          --numLinesOutput=2 \
          --IsOutputAFinCaseCtrl=TRUE
}
for i in `seq 1 22`;
do
    saige_step2 ${i} &> /hades/psivakumar/Emer_CH/gwas/SAIGE_step2_chr${i}.log &
done
