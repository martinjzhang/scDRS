#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-03:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=24000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./job.out # Standard output
#SBATCH -e ./job.err # Standard error

BATCH_NUM=$SLURM_ARRAY_TASK_ID
METHOD=$1
H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tabula-muris-senis-facs-official-raw-obj.h5ad
GS_FILE=out/gs_file/${METHOD}.top500.batch${BATCH_NUM}.gs
OUT_FOLDER=out/score_file/score.facs.${METHOD}.top500

mkdir -p ${OUT_FOLDER}

python /n/home11/mjzhang/gwas_informed_scRNAseq/scTRS/compute_score.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --gs_file $GS_FILE\
    --gs_species human\
    --out_folder $OUT_FOLDER\
    --flag_filter True \
    --flag_raw_count True 
