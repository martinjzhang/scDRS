H5AD_FILE=../scdrs/data/toydata_mouse.h5ad
OUT_FOLDER=./

# # existing way
# SCORE_FILE=../scdrs/data/res/@.full_score.gz
# python3 ../compute_downstream.py \
#     --h5ad_file $H5AD_FILE \
#     --score_file $SCORE_FILE \
#     --cell_type cell_type \
#     --cell_variable causal_variable,non_causal_variable,covariate \
#     --flag_gene True \
#     --flag_filter False \
#     --flag_raw_count False \
#     --out_folder $OUT_FOLDER

# # group analysis
# scdrs downstream \
#     --h5ad_file $H5AD_FILE \
#     --score-file $SCORE_FILE \
#     --group-analysis cell_type \
#     --out-folder cli/ \
#     --filter-data False \
#     --raw-count False

# # correlation analysis
# scdrs downstream \
#     --h5ad_file $H5AD_FILE \
#     --score-file $SCORE_FILE \
#     --corr-analysis causal_variable,non_causal_variable,covariate \
#     --out-folder cli/ \
#     --filter-data False \
#     --raw-count False

# # gene analysis
# scdrs downstream \
#     --h5ad_file $H5AD_FILE \
#     --score-file $SCORE_FILE \
#     --gene-analysis \
#     --out-folder cli/ \
#     --filter-data False \
#     --raw-count False
