#!/bin/bash
    
# # munge-gs --zscore-file
# scdrs munge-gs\
#     --zscore-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.zscore_file \
#     --out-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.from_zscore.gs \
#     --weight zscore\
#     --n-min 5\
#     --n-max 10

# # munge-gs --pval-file
# scdrs munge-gs\
#     --pval-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.pval_file \
#     --out-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.from_pval.gs \
#     --weight zscore\
#     --n-min 5\
#     --n-max 10

# scDRS CLI examples
H5AD_FILE=../scdrs/data/toydata_mouse.h5ad
COV_FILE=../scdrs/data/toydata_mouse.cov
GS_FILE_MOUSE=../scdrs/data/toydata_mouse.gs
GS_FILE_HUMAN=../scdrs/data/toydata_human.gs
OUT_FOLDER=../scdrs/data/res
SCORE_FILE=${OUT_FOLDER}/@.full_score.gz

# scdrs compute-score (for mouse gene set)
scdrs compute-score \
    --h5ad-file $H5AD_FILE\
    --h5ad-species mouse\
    --cov-file $COV_FILE\
    --gs-file $GS_FILE_MOUSE\
    --gs-species mouse\
    --ctrl-match-opt mean_var\
    --flag-filter-data False\
    --flag-raw_count False\
    --n-ctrl 20\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl_norm-score True\
    --out-folder $OUT_FOLDER

# scdrs compute-score (for human gene set)
scdrs compute-score \
    --h5ad-file $H5AD_FILE\
    --h5ad-species mouse\
    --cov-file $COV_FILE\
    --gs-file $GS_FILE_HUMAN\
    --gs-species human\
    --ctrl-match-opt mean_var\
    --flag-filter-data False\
    --flag-raw_count False\
    --n-ctrl 20\
    --flag-return-ctrl-raw-score False\
    --flag-return-ctrl_norm-score True\
    --out-folder $OUT_FOLDER

# perform-downstream --group-analysis
scdrs perform-downstream \
    --h5ad_file $H5AD_FILE \
    --score-file $SCORE_FILE \
    --group-analysis cell_type \
    --corr-analysis causal_variable,non_causal_variable,covariate \
    --gene-analysis \
    --out-folder $OUT_FOLDER \
    --flag-filter-data False \
    --flag-raw-count False \
    --knn-n-neighbors 15 \
    --knn-n-pcs 20
