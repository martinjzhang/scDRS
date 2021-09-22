#!/bin/bash
#SBATCH -n 2                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-6:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=32000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=1,2,6,10,11,19          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.out # Standard output
#SBATCH -e /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/experiments/job_info/job_%A_%a.err # Standard error

I_LINE=$SLURM_ARRAY_TASK_ID
# I_LINE=5

DATA_PATH=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/simulation_data
PERTURB_NAME=$(head -n $I_LINE ${DATA_PATH}/perturb_list.txt | tail -n 1)

GS_FILE=${DATA_PATH}/gs_file/all_ngene1000.gs
H5AD_FILE=${DATA_PATH}/single_cell_data/tms_facs.ncell_10k.h5ad
PERTURB_FILE=${DATA_PATH}/gs_file/${PERTURB_NAME}
OUT_FOLDER=${DATA_PATH}/score_file/${PERTURB_NAME}

[ -d ${OUT_FOLDER} ] || mkdir ${OUT_FOLDER}

python3 /n/home11/mjzhang/gwas_informed_scRNAseq/scDRS/compute_score_on_perturbed_data.py \
    --h5ad_file $H5AD_FILE\
    --h5ad_species mouse\
    --perturb_file $PERTURB_FILE\
    --gs_file $GS_FILE\
    --gs_species mouse\
    --flag_filter True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score False\
    --out_folder $OUT_FOLDER