#!/bin/bash

# scDRS CLI examples
H5AD_FILE=../scdrs/data/toydata_mouse.h5ad
COV_FILE=../scdrs/data/toydata_mouse.cov
GS_FILE_MOUSE=../scdrs/data/toydata_mouse.gs
GS_FILE_HUMAN=../scdrs/data/toydata_human.gs
OUT_FOLDER=../scdrs/data/res.from_scdrs_cli
SCORE_FILE=${OUT_FOLDER}/@.full_score.gz
    
# munge-gs
scdrs munge-gs\
    --zscore-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.zscore_file \
    --out-file /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/scdrs/data/gwas_gene.gs \
    --weight zscore\
    --n-min 5\
    --n-max 10

# # scdrs compute-score
# scdrs compute-score \
#     --h5ad-file $H5AD_FILE\
#     --h5ad-species mouse\
#     --cov-file $COV_FILE\
#     --gs-file $GS_FILE_MOUSE\
#     --gs-species mouse\
#     --ctrl-match-opt mean_var\
#     --flag-filter-data False\
#     --flag-raw_count False\
#     --n-ctrl 20\
#     --flag-return-ctrl-raw-score False\
#     --flag-return-ctrl_norm-score True\
#     --out-folder $OUT_FOLDER
    
# scdrs compute-score \
#     --h5ad-file $H5AD_FILE\
#     --h5ad-species mouse\
#     --cov-file $COV_FILE\
#     --gs-file $GS_FILE_HUMAN\
#     --gs-species human\
#     --ctrl-match-opt mean_var\
#     --flag-filter-data False\
#     --flag-raw_count False\
#     --n-ctrl 20\
#     --flag-return-ctrl-raw-score False\
#     --flag-return-ctrl_norm-score True\
#     --out-folder $OUT_FOLDER

# # perform-downstream --group-analysis
# scdrs perform-downstream \
#     --h5ad_file $H5AD_FILE \
#     --score-file $SCORE_FILE \
#     --group-analysis cell_type \
#     --corr-analysis causal_variable,non_causal_variable,covariate \
#     --gene-analysis \
#     --out-folder $OUT_FOLDER \
#     --flag-filter-data False \
#     --flag-raw-count False

# # perform-downstream --corr-analysis
# scdrs downstream \
#     --h5ad_file $H5AD_FILE \
#     --score-file $SCORE_FILE \
#     --corr-analysis causal_variable,non_causal_variable,covariate \
#     --out-folder cli/ \
#     --filter-data False \
#     --raw-count False

# # perform-downstream --gene-analysis
# scdrs downstream \
#     --h5ad_file $H5AD_FILE \
#     --score-file $SCORE_FILE \
#     --gene-analysis \
#     --out-folder cli/ \
#     --filter-data False \
#     --raw-count False
