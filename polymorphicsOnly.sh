#!/bin/bash

for i in `seq 1 22`; 
do 
	bcftools view -c1.0:minor -S ID_order_chrX.txt chr${i}.dose.vcf.gz -O z > chr${i}.imputed.poly.vcf.gz 2> err/chr${i}.imputed.poly.err &
done
