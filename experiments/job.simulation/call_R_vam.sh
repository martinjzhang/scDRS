#!/bin/bash
#SBATCH -n 2                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-6:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=32000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=4,8,12,16           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scTRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scTRS/experiments/job_info/job_%A_%a.err # Standard error

I_LINE=$SLURM_ARRAY_TASK_ID
# I_LINE=3
DNAME=tms_facs.ncell_10k

DATA_PATH=/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/simulation_data
GS_NAME=$(head -n $I_LINE ${DATA_PATH}/simu_list.txt | tail -n 1)
GS_FILE=${DATA_PATH}/gs_file/${GS_NAME}.gs
SC_FILE=${DATA_PATH}/single_cell_data/${DNAME}.tsv
OUT_FILE=${DATA_PATH}/score_file/result_vam/${DNAME}.${GS_NAME}.tsv

Rscript R_vam.R $SC_FILE $GS_FILE $OUT_FILE