#!/bin/bash
# $1 is input file, $2 is output prefix
plink --bfile $1 --extract commonSnpsToExtract.txt --make-bed --recode --out $2_4_commonSnps
