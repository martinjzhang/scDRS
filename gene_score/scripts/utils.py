from urllib.request import urlopen
import pandas as pd
from bs4 import BeautifulSoup as bs
from os.path import join
import time
import fire
from tqdm import tqdm
import os
import numpy as np
import glob
import submitit

from statsmodels.stats.multitest import multipletests
import scipy.stats

def process_magma(raw_dir, processed_dir):
    
    df = pd.read_csv(join(raw_dir, 'MAGMA_v108_GENE_10_ZSTAT.txt'), sep='\t')
    for trait in df.columns:
        trait_df = df[[trait]].copy()
        trait_df.columns = ['z_score']
        trait_df = trait_df.dropna()
        trait_df['p_val'] = 1 - scipy.stats.norm.cdf(trait_df['z_score'].values)
        trait_df['fdr'] = multipletests(trait_df['p_val'].values, method='fdr_bh')[1]
        trait_df = pd.DataFrame({'GENE': trait_df.index, 'FDR': trait_df.fdr})
        trait_df.to_csv(join(processed_dir, f'{trait}.csv'), index=False)


def process_gwas_maxabsz(sumstats_path, out_path, window_size, gene_loc_path, reference_template):
    # read reference
    reference_list = []
    for chr_i in range(1, 23):
        reference_list.append(pd.read_csv(reference_template.format(chr_i), sep='\t', usecols=[0, 1, 3], header=None, names=['CHR', 'SNP', 'BP']))
    reference = pd.concat(reference_list).reset_index(drop=True)

    # read gene location
    gene_loc = pd.read_csv(gene_loc_path, delim_whitespace=True, header=None, usecols=[1,2,3,5], names=['CHR', 'START', 'STOP', 'GENE'])
    gene_loc = gene_loc[gene_loc['CHR'].isin(np.arange(1, 23).astype(str))]
    gene_loc['CHR'] = gene_loc['CHR'].astype(int)

    def extract_gene_info(sumstats, gene, window_size):
        gene_index = (sumstats.CHR == gene.CHR) & (gene.START - window_size <= sumstats.BP) & (sumstats.BP < gene.STOP + window_size)
        gene_sumstats = sumstats[gene_index]
        return gene_sumstats['Z'].abs().max()

    
    sumstats = pd.read_csv(sumstats_path, sep='\t')
    sumstats = pd.merge(sumstats, reference, on='SNP')

    gene_loc['MAX_ABS_Z'] = gene_loc.apply(lambda gene : extract_gene_info(sumstats, gene, window_size), axis=1)

    gene_loc[['GENE', 'MAX_ABS_Z']].to_csv(out_path, na_rep='NaN', index=False)

def process_gwas_maxabsz_cli(raw_dir, processed_dir, window_size=10_000, 
                             sumstats_dir="/n/holystore01/LABS/price_lab/Lab/ldsc/sumstats_formatted",
                             reference_template = "/n/holystore01/LABS/price_lab/Lab/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.{}.bim"):
    sumstats_path_list = glob.glob(join(sumstats_dir, "*.sumstats"))
    trait_list = [path.split('/')[-1][:-9] for path in sumstats_path_list]
    out_path_list = [join(processed_dir, f"{trait}.csv") for trait in trait_list]
    gene_loc_path = join(raw_dir, "NCBI37.3.gene.loc")
    
    # submit the jobs
    executor = submitit.AutoExecutor(folder="~/submitit_log/")
    executor.update_parameters(timeout_min=10, mem_gb=8, slurm_partition="serial_requeue")
    wrapper = lambda param: process_gwas_maxabsz(sumstats_path=param[0], out_path=param[1],
                                                 window_size=window_size, gene_loc_path=gene_loc_path, reference_template = reference_template)
    param_list = [(sumstats_path, out_path) for (sumstats_path, out_path) in zip(sumstats_path_list, out_path_list) if not os.path.exists(out_path)]
    print(param_list)
    jobs = executor.map_array(wrapper, param_list)



def process_hess(raw_dir, processed_dir):
    
    trait_list = [f[:-4] for f in os.listdir(raw_dir) if f.endswith('.csv')]

    for trait in trait_list:
        df = pd.read_csv(join(raw_dir, f'{trait}.csv'))
        df = df[['GENE', 'ESTIMATE']]
        df.to_csv(join(processed_dir, f'{trait}.csv'), index=False, na_rep='NaN')

    
if __name__ == '__main__':
    fire.Fire()
