#!/bin/bash
#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-06:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=32000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0         # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.err # Standard error

H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/simulation_data/single_cell_data/tms_facs_5ct.h5ad
COV_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/tabula_muris_senis/tms_facs.cov
GS_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/gs_file/magma_10kb_top1000_zscore.75_traits.batch/batch0.gs
OUT_FOLDER=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/score_file/debug.adj_prop

scdrs compute-score \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE\
    --gs_species human\
    --adj-prop cell_ontology_class\
    --flag_filter_data True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score False\
    --out_folder $OUT_FOLDER