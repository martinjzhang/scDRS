import fire
import pandas as pd
from os.path import join
import numpy as np
from scipy import linalg
import fire
import glob
from os.path import join
from pysnptools.snpreader import Bed
from tqdm import tqdm

def impute_std(geno):
    """
    impute the mean and then standardize
    geno: (num_indivs x num_snps) numpy array
    """
    mean = np.nanmean(geno, axis=0)
    nanidx = np.where(np.isnan(geno))
    geno[nanidx] = mean[nanidx[1]]
    std = np.std(geno, axis=0)
    std_geno = (geno - mean) / std
    return std_geno

def compute_corr(geno):
    num_indivs, num_snps = geno.shape
    std_geno = impute_std(geno)
    return np.dot(std_geno.T, std_geno) / num_indivs

def hess_estimator(ld, sumstats, max_num_eig, min_eigval=1.0):
    ld_w, ld_v = linalg.eigh(ld)
    idx = ld_w.argsort()[::-1]
    ld_w = ld_w[idx]
    ld_v = ld_v[:,idx]
    
    beta_gwas = sumstats.Z.values / np.sqrt(sumstats.N.values)
    all_prj = np.zeros_like(ld_w)
    for i in range(len(ld_w)):
        prj = np.dot(beta_gwas.T, ld_v[:, i])
        all_prj[i] = prj ** 2
    
    eig_num = min(max_num_eig, np.where(ld_w > min_eigval)[0].size)
    
    quad_form = np.sum(all_prj[0:eig_num] / ld_w[0 : eig_num])
    
    indiv_num = sumstats.N.mean()
    est = (indiv_num * quad_form - eig_num) / (indiv_num - eig_num)
    var = ((indiv_num / (indiv_num - eig_num)) ** 2) * \
        (2 * eig_num * (1 - est) / indiv_num + 4 * est) * \
        (1 - est) / indiv_num
    return est, var

def hess_cli(bfile_template, sumstats_path, gene_loc_path, out, chr_i=None, max_num_eig=50, min_eigval=0.5):
    raw_sumstats = pd.read_csv(sumstats_path, sep='\t')
    raw_gene_loc = pd.read_csv(gene_loc_path, delim_whitespace=True, header=None, usecols=[1,2,3,5], names=['CHR', 'START', 'STOP', 'GENE'])
    raw_gene_loc = raw_gene_loc[raw_gene_loc['CHR'].isin(np.arange(1, 23).astype(str))]
    raw_gene_loc['CHR'] = raw_gene_loc['CHR'].astype(int)
    
    if chr_i is None:
        chr_i_list = np.arange(1, 23)
    else:
        chr_i_list = [chr_i]
    df_estimates = {'GENE': [], 'ESTIMATE': [], 'VAR': []}
    for chr_i in chr_i_list:
        bim = pd.read_csv(bfile_template.format(chr_i) + '.bim', sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
        geno = Bed(bfile_template.format(chr_i), count_A1=True)
        sumstats = pd.merge(bim, raw_sumstats, on='SNP', suffixes=['', '_sumstats'])

        noflip_rows = (sumstats.A1 == sumstats.A1_sumstats)
        flip_rows = (sumstats.A1 == sumstats.A2_sumstats)

        sumstats.loc[flip_rows, 'Z'] = (-1) * sumstats.loc[flip_rows, 'Z']
        sumstats = sumstats.loc[np.logical_or(noflip_rows, flip_rows), :].drop(columns=['A1_sumstats', 'A2_sumstats'])

        gene_loc = raw_gene_loc.loc[raw_gene_loc.CHR == chr_i].reset_index(drop=True)
    
        for gene_i, gene in gene_loc.iterrows():
            gene_snp_index = (gene.START <= sumstats.BP) & (sumstats.BP < gene.STOP)
            gene_sumstats = sumstats.loc[gene_snp_index, :]
            gene_geno = geno[:, np.isin(geno.sid, gene_sumstats.SNP)].read().val
            gene_ld = compute_corr(gene_geno)
            est, var = hess_estimator(gene_ld, gene_sumstats, max_num_eig=50)
            df_estimates['GENE'].append(gene.GENE)
            df_estimates['ESTIMATE'].append(est)
            df_estimates['VAR'].append(var)
    df_estimates = pd.DataFrame(df_estimates)
    df_estimates.to_csv(out, index=False)

if __name__ == '__main__':
    fire.Fire()