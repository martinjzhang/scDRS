#!/bin/bash
#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=60000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./job_info/job_%A_%a.out # Standard output
#SBATCH -e ./job_info/job_%A_%a.err # Standard error

BATCH_NUM=$SLURM_ARRAY_TASK_ID
GS_PREFIX=$1

# BATCH_NUM=0

# # TMS FACS
H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/tabula_muris_senis/tabula-muris-senis-facs-official-raw-obj.h5ad
COV_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/tabula_muris_senis/tms_facs.cov

# # TMS droplet
# H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tabula-muris-senis-droplet-official-raw-obj.h5ad
# COV_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tms_droplet.cov

GS_FILE=./geneset/${GS_PREFIX}.gs.batch/batch${BATCH_NUM}.gs

OUT_FOLDER=./score_file/${GS_PREFIX}
mkdir -p ${OUT_FOLDER}

python3 /n/home11/mjzhang/gwas_informed_scRNAseq/scTRS/compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --cov_file $COV_FILE\
    --gs_file $GS_FILE\
    --gs_species human\
    --flag_filter True\
    --flag_raw_count True\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score False\
    --out_folder $OUT_FOLDER
