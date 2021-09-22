#!/bin/bash
#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-06:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0         # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.err # Standard error

# Change this file to your local directory for scDRS_data_release_092121
DATA_PATH=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data_release_092121

BATCH_NUM=0
H5AD_FILE=${DATA_PATH}/tabula_muris_senis/tabula-muris-senis-facs-official-raw-obj.h5ad
COV_FILE=${DATA_PATH}/tabula_muris_senis/tms_facs.cov
GS_FILE=${DATA_PATH}/gs_file/magma_10kb_1000.74_traits.gs.batch/magma_10kb_1000.batch$BATCH_NUM.gs
OUT_FOLDER=${DATA_PATH}/score_file/score.tms_facs_with_cov.magma_10kb_1000.reproduce

# Change this to your local scDRS compute_score.py file
python3 /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE\
    --gs_species human\
    --flag_filter True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER
