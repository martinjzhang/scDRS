#!/bin/bash
#SBATCH -n 4                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=60000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-74          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./job_info/job_%A_%a.out # Standard output
#SBATCH -e ./job_info/job_%A_%a.err # Standard error


BATCH_NUM=$SLURM_ARRAY_TASK_ID
GS_PREFIX=$1
WEIGHT_OPT=$2


# # TMS FACS
H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/tabula_muris_senis/tabula-muris-senis-facs-official-raw-obj.h5ad
COV_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/tabula_muris_senis/tms_facs.cov


GS_FILE=./tms_facs/geneset/${GS_PREFIX}.gs.batch/batch${BATCH_NUM}.gs

OUT_FOLDER=./tms_facs/score_file/${GS_PREFIX}.$WEIGHT_OPT
mkdir -p ${OUT_FOLDER}

python3 /n/holystore01/LABS/price_lab/Users/khou/scDRS-revision/compute_score.py \
    --h5ad_file $H5AD_FILE \
    --h5ad_species mouse \
    --cov_file $COV_FILE \
    --weight_opt $WEIGHT_OPT \
    --gs_file $GS_FILE \
    --gs_species human\
    --flag_filter True\
    --flag_raw_count True\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER
