#!/bin/bash

for i in `seq 1 22`; 
do 
	bcftools view -q1.0:major -S ID_order_chrX.txt chr${i}.dose.vcf.gz -O v | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' >> imputed_monomorphic.txt
done
