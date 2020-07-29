import scanpy as sc
import numpy as np
import scipy as sp
from statsmodels.stats.multitest import multipletests
import pandas as pd

def score_cell(data, 
               gene_list, 
               suffix='',
               flag_correct_background=False,
               verbose=True,
               copy=False):
    
    """score cells based on the geneset 
    
    Args:
        data (AnnData): AnnData object
            adata.X should contain size-normalized log1p transformed count data
        gene_list (list): gene list 
        suffix (str): 'trs_'+suffix+['', '_z', '_p', '_bhp'] would be the name
        flag_correct_background (bool):
            If normalize for background mean and std. If True, normalize by 
            score = (score - mean)/std
        tissue (str): 'all' or one of the facs or droplet tissues
        

    Returns:
        adata (AnnData): Combined data for FACS and droplet
    """
    
    adata = data.copy() if copy else data
    
    gene_list_overlap = list(set(adata.var_names) & set(gene_list))
    if verbose:
        print('# score_cell: %d/%d gene_list genes also in adata'
              %(len(gene_list), len(gene_list_overlap)))
        print('# score_cell: suffix=%s, flag_correct_background=%s'
              %(suffix, flag_correct_background))
        
    trs_name = 'trs_%s'%suffix
    if trs_name in adata.obs.columns:
        print('# score_cell: overwrite original %s in adata.obs.columns'
              %trs_name)
        
    adata.obs[trs_name] = adata[:, gene_list_overlap].X.mean(axis=1)
    
    if flag_correct_background:
        v_mean,v_var = get_sparse_var(adata.X, axis=1)
        v_std = np.sqrt(v_var)
        adata.obs[trs_name] = (adata.obs[trs_name] - v_mean) / v_std * \
                                np.sqrt(len(gene_list_overlap))
        
    # Add z_score, p_value, and fdr 
    temp_v = adata.obs[trs_name].values
    adata.obs['%s_z'%trs_name] = (temp_v - temp_v.mean())/ temp_v.std()
    adata.obs['%s_p'%trs_name] = 1 - sp.stats.norm.cdf(adata.obs['%s_z'%trs_name].values)
    adata.obs['%s_bhp'%trs_name] = multipletests(adata.obs['%s_p'%trs_name].values,
                                                 method='fdr_bh')[1]
    
    return adata if copy else None


def score_cell_kangcheng_072920(data, 
               gene_list, 
               suffix='',
               flag_correct_background=False,
               flag_specific_expressed=False,
               verbose=True,
               copy=False):
    
    """score cells based on the geneset 
    
    Args:
        data (AnnData): AnnData object
            adata.X should contain size-normalized log1p transformed count data
        gene_list (list): gene list 
        suffix (str): 'trs_'+suffix+['', '_z', '_p', '_bhp'] would be the name
        flag_correct_background (bool):
            If normalize for background mean and std per_cell. If True, normalize by 
            score = (score - mean)/std, where mean and std is calculated within each cell
        flag_specific_expressed (bool):
            Whether transform gene expression to identify specific expressed genes.
            If True, for each gene, normalize score = (score - mean) / std, where mean and 
            std is calculated across the cells when calculating the TRS score, 
        tissue (str): 'all' or one of the facs or droplet tissues
        

    Returns:
        adata (AnnData): Combined data for FACS and droplet
    """
    
    adata = data.copy() if copy else data
    
    gene_list_overlap = list(set(adata.var_names) & set(gene_list))
    if verbose:
        print('# score_cell: %d/%d gene_list genes also in adata'
              %(len(gene_list), len(gene_list_overlap)))
        print('# score_cell: suffix=%s, flag_correct_background=%s, flag_specific_expressed=%s'
              %(suffix, flag_correct_background, flag_specific_expressed))
        
    trs_name = 'trs_%s'%suffix
    if trs_name in adata.obs.columns:
        print('# score_cell: overwrite original %s in adata.obs.columns'
              %trs_name)
        
    adata.obs[trs_name] = adata[:, gene_list_overlap].X.mean(axis=1)
    
    if flag_correct_background:
        cell_mean,cell_var = get_sparse_var(adata.X, axis=1)
        cell_std = np.sqrt(cell_var)
        # reshape to (1, #cells) vector
        cell_mean = cell_mean[:, np.newaxis]
        cell_std = cell_std[:, np.newaxis]
    
    gwas_mat = adata[:, gene_list_overlap].X
    if flag_correct_background:
        # normalize for each cell
        gwas_mat = (gwas_mat - cell_mean) / cell_std
        
    if flag_specific_expressed:
        # normalize for each gene
        gene_mean, gene_std = np.mean(gwas_mat, axis=0), np.std(gwas_mat, axis=0)
        gwas_mat = (gwas_mat - gene_mean) / gene_std
    
    adata.obs[trs_name] = gwas_mat.mean(axis=1)
        
    # Add z_score, p_value, and fdr 
    temp_v = adata.obs[trs_name].values
    adata.obs['%s_z'%trs_name] = (temp_v - temp_v.mean())/ temp_v.std()
    adata.obs['%s_p'%trs_name] = 1 - sp.stats.norm.cdf(adata.obs['%s_z'%trs_name].values)
    adata.obs['%s_bhp'%trs_name] = multipletests(adata.obs['%s_p'%trs_name].values,
                                                 method='fdr_bh')[1]
    
    return adata if copy else None

def get_sparse_var(sparse_X, axis=0):
    
    v_mean = sparse_X.mean(axis=axis)
    v_mean = np.array(v_mean).reshape([-1])
    v_var = sparse_X.power(2).mean(axis=axis)
    v_var = np.array(v_var).reshape([-1])
    v_var = v_var - v_mean**2
    
    
    return v_mean,v_var
    
    