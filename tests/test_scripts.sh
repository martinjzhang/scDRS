#!/bin/bash

H5AD_FILE=../scdrs/data/toydata_mouse.h5ad
COV_FILE=../scdrs/data/toydata_mouse.cov
GS_FILE_MOUSE=../scdrs/data/toydata_mouse.gs
GS_FILE_HUMAN=../scdrs/data/toydata_human.gs
OUT_FOLDER=../scdrs/data/res

python3 ../compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE_MOUSE\
    --gs_species mouse\
    --ctrl_match_opt mean_var\
    --flag_filter False\
    --flag_raw_count False\
    --n_ctrl 20\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER
    
python3 ../compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE_HUMAN\
    --gs_species human\
    --ctrl_match_opt mean_var\
    --flag_filter False\
    --flag_raw_count False\
    --n_ctrl 20\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER

SCORE_FILE=../scdrs/data/res/@.full_score.gz
python3 ../compute_downstream.py \
    --h5ad_file $H5AD_FILE\
    --score_file $SCORE_FILE\
    --cell_type cell_type\
    --cell_variable causal_variable,non_causal_variable,covariate\
    --flag_gene True\
    --flag_filter False\
    --flag_raw_count False\
    --out_folder $OUT_FOLDER
