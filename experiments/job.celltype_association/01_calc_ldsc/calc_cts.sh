#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-02:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=16G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_err  # File to which STDERR will be written, %j inserts jobid

source activate ldsc

dset=$1
ref_dir="/n/holystore01/LABS/price_lab/Lab/ldsc/reference_files/1000G_EUR_Phase3/"
ldsc_dir=ldsc
sumstats_dir="/n/holystore01/LABS/price_lab/Lab/ldsc/sumstats_formatted"

trait=$(sed "${SLURM_ARRAY_TASK_ID}q;d" out/trait_list.txt)
mkdir -p out/${dset}/cts_result
python ${ldsc_dir}/ldsc.py \
    --h2-cts ${sumstats_dir}/${trait}.sumstats \
    --ref-ld-chr ${ref_dir}/baseline_v1.2/baseline. \
    --out out/${dset}/cts_result/${trait} \
    --ref-ld-chr-cts out/${dset}/ldsc.ldcts \
    --w-ld-chr ${ref_dir}/weights/weights.hm3_noMHC.
