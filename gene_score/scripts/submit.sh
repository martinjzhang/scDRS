#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-00:20          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=8G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_err  # File to which STDERR will be written, %j inserts jobid

export PATH=/n/holystore01/LABS/price_lab/Users/khou/miniconda3/bin:$PATH
export PYTHONNOUSERSITE=True

trait=$1

python utils.py process_gene_max_abs_z --raw_dir None --processed_dir out/processed/gwas_max_abs_z --trait ${trait}