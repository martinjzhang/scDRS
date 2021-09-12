#!/bin/bash
#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-05:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=80000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-25          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./job_out.zeisel_2018 # Standard output
#SBATCH -e ./job_out.zeisel_2018 # Standard error

BATCH_NUM=$SLURM_ARRAY_TASK_ID
H5AD_FILE=../00_prepare_dataset/processed/zeisel_2018.raw.h5ad
COV_FILE=../00_prepare_dataset/processed/zeisel_2018.cov.tsv

GS_FILE=gs_file/mouse.gs.batch/batch${BATCH_NUM}.gs

OUT_FOLDER=./score_file/zeisel_2018

mkdir -p ${OUT_FOLDER}

python3 /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --gs_file $GS_FILE\
    --gs_species mouse\
    --cov_file ${COV_FILE} \
    --flag_filter True\
    --flag_raw_count True\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER
