#/bin/bash

conda activate ldsc

python2 ldsc/ldsc.py
        --h2 BRCA/BRCA_OncoArray.sumstats.gz
        --ref-ld-chr BRCA/brca-drivers/brca.,all-drivers/all.,ATAC/BRCA-ocr-annots/brca-ocr.,ATAC/immune-noBreast-annots/merged.,h3k27ac/merged-annots/h3k27ac.,partition/baseline/baseline.
        --w-ld-chr partition/weights_hm3_no_hla/weights.
        --frqfile-chr partition/1000G_frq/1000G.mac5eur.
        --out BRCA/BRCA_all -
        --overlap-annot
        --print-coefficients
