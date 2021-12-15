#!/bin/bash

scdrs munge-gs\
    --zscore-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.zscore_file \
    --out-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.gs \
    --weight zscore\
    --n-min 5\
    --n-max 10