#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-03:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=20gb           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --array=0-74          # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ./job_info/job_%A_%a.out # Standard output
#SBATCH -e ./job_info/job_%A_%a.err # Standard error

# SLURM_ARRAY_TASK_ID=1 bash compute_score.sh  
BATCH_NUM=$SLURM_ARRAY_TASK_ID

# GWAS gene set prefix
GS_PREFIX=$1
# scRNA-seq weight option
WEIGHT_OPT=$2

COV_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/tabula_muris_senis/tms_facs.cov
GS_FILE=geneset/${GS_PREFIX}.gs.batch/batch${BATCH_NUM}.gs

for rep_i in $(seq 0 19); do
    echo "rep${rep_i}"
    H5AD_FILE=/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/simulation_data/single_cell_data/tms_facs.ncell_10k_rep${rep_i}.h5ad
    OUT_FOLDER=~/holyscratch/continuous_weights_experiment/score_file/rep${rep_i}/${GS_PREFIX}.${WEIGHT_OPT}
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
done
