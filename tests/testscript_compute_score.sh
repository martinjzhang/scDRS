#!/bin/bash

H5AD_FILE=../scdrs/data/toydata_mouse.h5ad
COV_FILE=../scdrs/data/toydata_mouse.cov
GS_FILE=../scdrs/data/toydata_mouse.gs
OUT_FOLDER=../scdrs/data

python3 ../compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE\
    --gs_species mouse\
    --ctrl_match_opt mean_var\
    --flag_filter False\
    --flag_raw_count False\
    --n_ctrl 20\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score False\
    --out_folder $OUT_FOLDER
