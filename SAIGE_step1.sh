logfile_step1=.Phenotype_step1.log { Rscript 
/hades/Software/NGS_Software/SAIGE/extdata/step1_fitNULLGLMM.R \
        --plinkFile=/hades/psivakumar/Emer_CH/gwas/saige_all_26 \
        --phenoFile=/hades/psivakumar/Emer_CH/gwas/saige.pheno \
        --phenoCol=Phenotype \
        --covarColList=PC1,PC2,PC3,PC4 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --skipModelFitting=FALSE \
        --outputPrefix=/hades/psivakumar/Emer_CH/gwas/SAIGE_step1 \
        --nThreads=2
} 2>&1 | tee $logfile_step1
