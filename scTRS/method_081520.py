import scanpy as sc
import numpy as np
import scipy as sp
from statsmodels.stats.multitest import multipletests
import pandas as pd
import time

def score_cell(data, 
               gene_list, 
               suffix='',
               flag_correct_background=False,
               flag_nullgene=False,
               flag_permute_cell=False,
               random_seed=0,
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
    
    np.random.seed(random_seed)
    
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
    
    # Compute TRS
    adata.obs[trs_name] = adata[:, gene_list_overlap].X.mean(axis=1)
    if flag_nullgene:
        
        # A random set
        ind_select = np.random.permutation(adata.shape[1])[:len(gene_list_overlap)]
        gene_list_null = list(adata.var_names[ind_select])
        adata.obs['%s_nullgene'%trs_name] = adata[:, gene_list_null].X.mean(axis=1)
        
        # A random set with matching mean expression 
        temp_df = pd.DataFrame(index=adata.var_names)
        temp_df['gene'] = temp_df.index
        temp_df['mean'] = np.array(adata.X.mean(axis=0)).reshape([-1])
        temp_df['mean_bin'] = pd.cut(temp_df['mean'], bins=100)
        temp_df_bin = temp_df.groupby('mean_bin').agg({'gene':set})
        temp_df_bin = temp_df_bin.loc[~temp_df_bin['gene'].isna()]
        
        gene_set_overlap = set(gene_list_overlap)
        gene_list_null_me = []
        for bin_ in temp_df_bin.index:
            temp_n_gene = len(temp_df_bin.loc[bin_,'gene']&gene_set_overlap)
            temp_v_gene = np.array(list(temp_df_bin.loc[bin_,'gene']-gene_set_overlap))
            if (temp_n_gene>0) & (temp_v_gene.shape[0]>0):
                gene_list_null_me += list(temp_v_gene[np.random.permutation(temp_v_gene.shape[0])[0:temp_n_gene]])
        if verbose:
            print('# score_cell: %d trait genes with mean_exp=%0.3f'
                   %(len(gene_set_overlap), temp_df.loc[gene_set_overlap,'mean'].values.mean()))
            print('# score_cell: %d null_me genes with mean_exp=%0.3f'
                   %(len(gene_list_null_me), temp_df.loc[gene_list_null_me,'mean'].values.mean()))
        adata.obs['%s_nullgeneme'%trs_name] = adata[:, gene_list_null_me].X.mean(axis=1)
        
    if flag_permute_cell:
        temp_X = adata[:, gene_list_overlap].X.copy()
        temp_X = sp.sparse.csc_matrix(temp_X)
        for i_col in np.arange(temp_X.shape[1]):
            ind_perm = np.random.permutation(temp_X.shape[0])
            temp_X[:,i_col] = temp_X[ind_perm,i_col]
        adata.obs['%s_permcell'%trs_name] = temp_X.mean(axis=1)
    
    # Cell background correction
    if flag_correct_background:
        v_mean,v_var = get_sparse_var(adata.X, axis=1)
        v_std = np.sqrt(v_var)
        adata.obs[trs_name] = (adata.obs[trs_name] - v_mean) / v_std * \
                                np.sqrt(len(gene_list_overlap))
        if flag_nullgene:
            adata.obs['%s_nullgene'%trs_name] = (adata.obs['%s_nullgene'%trs_name] - v_mean) / v_std * \
                                np.sqrt(len(gene_list_null))
            adata.obs['%s_nullgeneme'%trs_name] = (adata.obs['%s_nullgeneme'%trs_name] - v_mean) / v_std * \
                                np.sqrt(len(gene_list_null_me))
        if flag_permute_cell:
            adata.obs['%s_permcell'%trs_name] = (adata.obs['%s_permcell'%trs_name] - v_mean) / v_std * \
                                np.sqrt(len(gene_list_overlap))
        
    # Add z_score, p_value, and fdr 
    temp_v = adata.obs[trs_name].values
    adata.obs['%s_z'%trs_name] = (temp_v - temp_v.mean())/ temp_v.std()
    adata.obs['%s_p'%trs_name] = 1 - sp.stats.norm.cdf(adata.obs['%s_z'%trs_name].values)
    adata.obs['%s_bhp'%trs_name] = multipletests(adata.obs['%s_p'%trs_name].values,
                                                 method='fdr_bh')[1]
    
    
    return adata if copy else None

# def score_cell(data, 
#                gene_list, 
#                suffix='',
#                flag_correct_background=False,
#                verbose=True,
#                copy=False):
#     
#     """score cells based on the geneset 
#     
#     Args:
#         data (AnnData): AnnData object
#             adata.X should contain size-normalized log1p transformed count data
#         gene_list (list): gene list 
#         suffix (str): 'trs_'+suffix+['', '_z', '_p', '_bhp'] would be the name
#         flag_correct_background (bool):
#             If normalize for background mean and std. If True, normalize by 
#             score = (score - mean)/std
#         tissue (str): 'all' or one of the facs or droplet tissues
#         
# 
#     Returns:
#         adata (AnnData): Combined data for FACS and droplet
#     """
#     
#     adata = data.copy() if copy else data
#     
#     gene_list_overlap = list(set(adata.var_names) & set(gene_list))
#     if verbose:
#         print('# score_cell: %d/%d gene_list genes also in adata'
#               %(len(gene_list), len(gene_list_overlap)))
#         print('# score_cell: suffix=%s, flag_correct_background=%s'
#               %(suffix, flag_correct_background))
#         
#     trs_name = 'trs_%s'%suffix
#     if trs_name in adata.obs.columns:
#         print('# score_cell: overwrite original %s in adata.obs.columns'
#               %trs_name)
#         
#     adata.obs[trs_name] = adata[:, gene_list_overlap].X.mean(axis=1)
#     
#     if flag_correct_background:
#         v_mean,v_var = get_sparse_var(adata.X, axis=1)
#         v_std = np.sqrt(v_var)
#         adata.obs[trs_name] = (adata.obs[trs_name] - v_mean) / v_std * \
#                                 np.sqrt(len(gene_list_overlap))
#         
#     # Add z_score, p_value, and fdr 
#     temp_v = adata.obs[trs_name].values
#     adata.obs['%s_z'%trs_name] = (temp_v - temp_v.mean())/ temp_v.std()
#     adata.obs['%s_p'%trs_name] = 1 - sp.stats.norm.cdf(adata.obs['%s_z'%trs_name].values)
#     adata.obs['%s_bhp'%trs_name] = multipletests(adata.obs['%s_p'%trs_name].values,
#                                                  method='fdr_bh')[1]
#     
#     return adata if copy else None


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

def gearys_c(adata, val_obs, prefix, stratify_obs=None, copy=False):
    """
    Interface of computing Geary's C statistics
    
    Args:
        adata: Anndata object
        val_obs: the obs name to calculate this statistics
        prefix: the name will be `prefix`_gearys_C
        stratify_obs: Calculate the statistics using `stratify_obs` obs column,
            must be a categorical variable
    """
    
    adata = adata.copy() if copy else adata
    if stratify_obs is not None:
        assert adata.obs[stratify_obs].dtype.name == 'category', \
            "`stratify_obs` must correspond to a Categorical column"
        categories = adata.obs[stratify_obs].unique()
        all_c_stats = np.zeros(adata.shape[0])
        for cat in categories:
            s_index = adata.obs[stratify_obs] == cat
            all_c_stats[s_index] = _gearys_c(adata[s_index], adata[s_index].obs[val_obs])
        
    else:
        all_c_stats = _gearys_c(adata, adata.obs[val_obs])
    
    gearys_C_name = prefix + '_gearys_C'
    if gearys_C_name in adata.obs.columns:
        print('# gearys_c: overwrite original %s in adata.obs.columns'
              %gearys_C_name)
    adata.obs[gearys_C_name] = all_c_stats
    # adata.obs[gearys_C_name] = adata.obs[gearys_C_name].astype('category')

    return adata if copy else None

    
def _gearys_c(adata, vals):
    """Compute Geary's C statistics for an AnnData
    Adopted from https://github.com/ivirshup/scanpy/blob/metrics/scanpy/metrics/_gearys_c.py
    
    C =
        \frac{
            (N - 1)\sum_{i,j} w_{i,j} (x_i - x_j)^2
        }{
            2W \sum_i (x_i - \bar{x})^2
        }
    
    Args:
        adata (AnnData): AnnData object
            adata.obsp["Connectivities] should contain the connectivity graph, 
            with shape `(n_obs, n_obs)`
        vals (Array-like):
            Values to calculate Geary's C for. If one dimensional, should have
            shape `(n_obs,)`. 
    Returns:
        C: the Geary's C statistics
    """
    graph = adata.obsp["connectivities"]
    assert graph.shape[0] == graph.shape[1]
    graph_data = graph.data.astype(np.float_, copy=False)
    assert graph.shape[0] == vals.shape[0]
    assert(np.ndim(vals) == 1)
    
    W = graph_data.sum()
    N = len(graph.indptr) - 1
    vals_bar = vals.mean()
    vals = vals.astype(np.float_)
    
    # numerators
    total = 0.0
    for i in range(N):
        s = slice(graph.indptr[i], graph.indptr[i + 1])
        # indices of corresponding neighbors
        i_indices = graph.indices[s]
        # corresponding connecting weights
        i_data = graph_data[s]
        total += np.sum(i_data * ((vals[i] - vals[i_indices]) ** 2))

    numer = (N - 1) * total
    denom = 2 * W * ((vals - vals_bar) ** 2).sum()
    C = numer / denom
        
    return C

def get_sparse_var(sparse_X, axis=0):
    
    v_mean = sparse_X.mean(axis=axis)
    v_mean = np.array(v_mean).reshape([-1])
    v_var = sparse_X.power(2).mean(axis=axis)
    v_var = np.array(v_var).reshape([-1])
    v_var = v_var - v_mean**2
    
    
    return v_mean,v_var
    
    