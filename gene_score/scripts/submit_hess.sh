#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-02:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=16G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_err  # File to which STDERR will be written, %j inserts jobid

export PATH=/n/holystore01/LABS/price_lab/Users/khou/miniconda3/bin:$PATH
export PYTHONNOUSERSITE=True

trait=$1

python hess.py hess_cli --bfile_template "/n/holystore01/LABS/price_lab/Lab/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.{}" \ 
                        --sumstats_path /n/holystore01/LABS/price_lab/Lab/ldsc/sumstats_formatted/${trait}.sumstats \
                        --gene_loc_path /n/home12/khou/scTRS/gene_sets/out/raw/gwas_derived/NCBI37.3.gene.loc \
                        --out out/${trait}.csv