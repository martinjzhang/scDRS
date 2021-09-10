#!/bin/bash
#SBATCH -n 1                # Number of cores (-n)
#SBATCH -N 1                # Ensure that all cores are on one Node (-N)
#SBATCH -t 0-04:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=16G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o job_out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e job_err  # File to which STDERR will be written, %j inserts jobid

source activate ldsc

dset=$1
prefix=$(sed "${SLURM_ARRAY_TASK_ID}q;d" out/${dset}/ct_list.txt)

kg_dir="/n/holystore01/LABS/price_lab/Lab/ldsc/reference_files/1000G_EUR_Phase3/"
ldsc_dir=ldsc

for chr_i in $(seq 1 22); do
    # make annotation
    python ${ldsc_dir}/make_annot.py \
        --gene-set-file out/${dset}/ldscore/${prefix}.geneset \
        --gene-coord-file out/gene_loc/gene_loc.txt \
        --windowsize 100000 \
        --bimfile ${kg_dir}/plink_files/1000G.EUR.QC.${chr_i}.bim \
        --annot-file out/${dset}/ldscore/${prefix}.${chr_i}.annot.gz

    # calculate LD score
    python ${ldsc_dir}/ldsc.py \
        --l2 \
        --bfile ${kg_dir}/plink_files/1000G.EUR.QC.${chr_i} \
        --ld-wind-cm 1 \
        --annot out/${dset}/ldscore/${prefix}.${chr_i}.annot.gz \
        --thin-annot \
        --out out/${dset}/ldscore/${prefix}.${chr_i} \
        --print-snps ${kg_dir}/list.txt
done