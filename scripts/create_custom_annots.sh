#!/bin/bash
​
conda activate ldsc
​
chr=$(seq 22)
​
bfile=$1
outDir=$2
prefix=$3
​
mkdir $outDir
​
for c in $chr
do
    python2 ../make_annot.py \
	   --bed-file $bfile \
	   --bimfile 1000G_plinkfiles/1000G.mac5eur."$c".bim \
	   --annot-file $outDir/$prefix."$c".annot.gz
​
    python2 ../ldsc.py \
           --l2 \
           --bfile 1000G_plinkfiles/1000G.mac5eur."$c" \
	   --ld-wind-cm 1 \
           --annot $outDir/$prefix."$c".annot.gz \
           --thin-annot \
           --out $outDir/$prefix."$c" \
           --print-snps hapmap3_snps/hm."$c".snp
done
