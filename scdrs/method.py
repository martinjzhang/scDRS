import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
from skmisc.loess import loess

def score_cell(data, 
               gene_list, 
               gene_weight=None,
               ctrl_match_key='mean',
               n_ctrl=500,
               n_genebin=200,
               weight_opt='vs',
               copy=False,
               return_ctrl_raw_score=False,
               return_ctrl_norm_score=False,
               random_seed=0,
               verbose=False,
               save_intermediate=None):
    
    """Score cells based on the trait gene set 
    
    Args
    ----
        data (n_cell, n_gene) : AnnData
            data.X should contain size-normalized log1p transformed count data
        gene_list (n_trait_gene) : list
            Trait gene list       
        gene_weight (n_trait_gene) : list/np.ndarray
            Gene weights for genes in the gene_list. 
            If gene_weight=None, the weigts are set to be one. 
        ctrl_match_key : str
            The quantity for matching control and trait genes. 
            ctrl_match_key should appear in data.obs.columns
        n_ctrl : int 
            Number of control genesets
        n_genebin : int
            Number of bins for dividing genes by ctrl_match_key
        weight_opt : str 
            Option for computing the raw score
            'uniform': average over the genes in the gene_list
            'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct)
            'inv_std': weighted average with weights equal to 1/std
            'od': overdispersion score
        copy : bool
            If to make copy of the AnnData object to avoid writing on the orignal data
        return_raw_ctrl_score : bool
            If to return control scores  
        return_norm_ctrl_score : bool
            If to return control scores  
        random_seed : int
            Random seed
        verbose : bool
            If to output messages 
        save_intermediate : str
            File path prefix for saving intermediate results

    Returns
    -------
        df_res (n_cell, n_key) : pd.DataFrame (dtype=np.float32)
            Columns: 
            0. raw_score
            1. norm_score: scores after cell-wise and trait-wise normalization
            2. mc_pval: p-values computed using only the control scores from the same cell 
            3. pval
            4. nlog10_pval: -log10(pval). Needed in case the single precision (np.float32) gives inaccurate p-values
            5. zscore: one-side z-score converted from pval
            
    TODO
    -------
        1. 
    """
    
    np.random.seed(random_seed)
    adata = data.copy() if copy else data
    n_cell,n_gene = adata.shape
    
    # Pre-compute statistics
    var_set = {'mean','var','var_tech'}
    obs_set = {'mean','var'}
    if (len(var_set-set(adata.var.columns))>0) | (len(obs_set-set(adata.obs.columns))>0):
        if verbose: 
            print('# score_cell: recompute statistics using method.compute_stats')
        compute_stats(adata)
    
    # Check options
    if ctrl_match_key not in adata.var.columns:
        raise ValueError('# score_cell: %s not in data.var.columns'%ctrl_match_key)
    weight_opt_list = ['uniform', 'vs', 'inv_std', 'od']
    if weight_opt not in weight_opt_list:
        raise ValueError('# score_cell: weight_opt not in [%s]'
                         %', '.join([str(x) for x in weight_opt_list]))
        
    if verbose:
        print('# score_cell: ctrl_match_key=%s, n_ctrl=%d, n_genebin=%d, weight_opt=%s'
              %(ctrl_match_key, n_ctrl, n_genebin, weight_opt))
        
    # Gene information
    df_gene = pd.DataFrame(index = adata.var_names, data={'gene':adata.var_names,
                                                          'mean':adata.var['mean'].values,
                                                          'var':adata.var['var'].values})
    if ctrl_match_key not in df_gene.columns:
        df_gene[ctrl_match_key] = adata.var[ctrl_match_key].values
    df_gene.drop_duplicates(subset='gene', inplace=True)
            
    # Update gene_list and gene_weight
    n_gene_old = len(list(gene_list))
    if gene_weight is None:
        gene_weight = [1]*n_gene_old
    dic_weight = {x:y for x,y in zip(gene_list, gene_weight)}
    gene_list = sorted(list(set(gene_list) & set(df_gene['gene'].values)))
    gene_weight = [dic_weight[x] for x in gene_list]
    
    if verbose: 
        print('# score_cell: %-15s %-15s %-20s'
              %('trait gene set,', '%d/%d genes,'%(len(gene_list),n_gene_old),
                'mean=%0.2e'%df_gene.loc[gene_list, 'mean'].mean()))
    
    # Select control gene sets
    dic_ctrl_list,dic_ctrl_weight = _select_ctrl_geneset(df_gene, gene_list, gene_weight,
                                                         ctrl_match_key, n_ctrl, n_genebin, random_seed)  

    # Compute raw scores        
    v_raw_score,v_score_weight = _compute_raw_score(adata, gene_list, gene_weight, weight_opt)
    
    mat_ctrl_raw_score,mat_ctrl_weight = np.zeros([n_cell,n_ctrl]),np.zeros([len(gene_list),n_ctrl])
    for i_ctrl in range(n_ctrl): 
        mat_ctrl_raw_score[:,i_ctrl],mat_ctrl_weight[:,i_ctrl] = \
            _compute_raw_score(adata, dic_ctrl_list[i_ctrl], dic_ctrl_weight[i_ctrl], weight_opt)
                
    # Compute normalized scores 
    v_var_ratio_c2t = np.ones(n_ctrl) 
    if weight_opt in ['uniform', 'vs', 'inv_std']: 
        # For raw scores compuated as weighted average. estimate variance ratio assuming independence
        for i_ctrl in range(n_ctrl):
            v_var_ratio_c2t[i_ctrl] = (df_gene.loc[dic_ctrl_list[i_ctrl], 'var'] * 
                                       mat_ctrl_weight[:,i_ctrl]**2).sum()
        v_var_ratio_c2t /= (df_gene.loc[gene_list, 'var']*v_score_weight**2).sum()
    
    v_norm_score,mat_ctrl_norm_score = _correct_background(v_raw_score, mat_ctrl_raw_score, v_var_ratio_c2t,
                                                           save_intermediate=save_intermediate)
    
    # Get p-values 
    mc_p = (1+(mat_ctrl_norm_score.T>=v_norm_score).sum(axis=0))/(1+n_ctrl)
    pooled_p = _get_p_from_empi_null(v_norm_score, mat_ctrl_norm_score.flatten())  
    nlog10_pooled_p = -np.log10(pooled_p)
    pooled_z = -sp.stats.norm.ppf(pooled_p).clip(min=-10,max=10)
    
    # Return result
    dic_res = {'raw_score':v_raw_score, 'norm_score':v_norm_score, 'mc_pval':mc_p,
               'pval':pooled_p, 'nlog10_pval':nlog10_pooled_p, 'zscore':pooled_z}
    if return_ctrl_raw_score:
        for i in range(n_ctrl):
            dic_res['ctrl_raw_score_%d'%i] = mat_ctrl_raw_score[:,i]
    if return_ctrl_norm_score:
        for i in range(n_ctrl):
            dic_res['ctrl_norm_score_%d'%i] = mat_ctrl_norm_score[:,i]
    df_res = pd.DataFrame(index=adata.obs.index,  data=dic_res, dtype=np.float32)
    return df_res


def _select_ctrl_geneset(input_df_gene, gene_list, gene_weight, 
                         ctrl_match_key, n_ctrl, n_genebin, random_seed):
    
    """Subroutine for score_cell, select control genesets 
    
    Args
    ----
        input_df_gene (adata.shape[1], n_statistic) : pd.DataFrame
            Gene-wise statistics
        gene_list (n_trait_gene) : list
            Trait gene list       
        gene_weight (n_trait_gene) : list/np.ndarray
            Gene weights for genes in the gene_list. 
        ctrl_match_key : str
            The quantity for matching control and trait genes. 
            ctrl_match_key should appear in input_df_gene.columns
        n_ctrl : int 
            Number of control gene sets
        n_genebin : int
            Number of bins for dividing genes by ctrl_match_key
        random_seed : int
            Random seed  
    Returns
    -------
        dic_ctrl_list : dictionary of lists
            dic_ctrl_list[i]: the i-th control gene list
        dic_ctrl_weight : dictionary of lists
            dic_ctrl_weight[i]: weights for the i-th control gene list
            
    """
    
    np.random.seed(random_seed)    
    df_gene = input_df_gene.copy()
    trait_gene_set = set(gene_list)
    dic_weight = {x:y for x,y in zip(gene_list,gene_weight)}
    
    # Divide genes into equal-sized bins based on ctrl_match_key
    if df_gene[ctrl_match_key].unique().shape[0] < df_gene.shape[0]/10:
        df_gene_bin = df_gene.groupby(ctrl_match_key).agg({'gene':set})    
    else:
        df_gene['qbin'] = pd.qcut(df_gene[ctrl_match_key], q=n_genebin, 
                                  labels=False, duplicates='drop')
        df_gene_bin = df_gene.groupby('qbin').agg({'gene':set})        
    
    # Find ctrl_match_key matched control genes
    dic_ctrl_list = {x:[] for x in range(n_ctrl)}
    dic_ctrl_weight = {x:[] for x in range(n_ctrl)}
    for bin_ in df_gene_bin.index:
        bin_gene = sorted(df_gene_bin.loc[bin_, 'gene'])
        bin_trait_gene = sorted(df_gene_bin.loc[bin_,'gene'] & trait_gene_set)
        if len(bin_trait_gene)>0:
            for i_list in np.arange(n_ctrl):
                dic_ctrl_list[i_list].extend(np.random.choice(bin_gene, size=len(bin_trait_gene),
                                                              replace=False))
                dic_ctrl_weight[i_list].extend([dic_weight[x] for x in bin_trait_gene])
                
    return dic_ctrl_list,dic_ctrl_weight


def _compute_raw_score(adata, gene_list, gene_weight, weight_opt):
    """Compute raw score
        v_score_weight = gene_weight * {uniform/vs/inv_std}
    
    Args
    ----
        adata (n_cell, n_gene) : AnnData
            adata.X should contain size-normalized log1p transformed count data
        gene_list (n_trait_gene) : list
            Trait gene list
        gene_weight (n_trait_gene) : list/np.ndarray
            Gene weights for genes in the gene_list 
        weight_opt : str 
            Option for computing the raw score
            - 'uniform': average over the genes in the gene_list
            - 'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct)
            - 'inv_std': weighted average with weights equal to 1/std
            - 'od': overdispersion score
            
    Returns
    -------
        v_raw_score (n_cell,) : np.ndarray
            Raw score
        v_score_weight (n_trait_gene,) : np.ndarray
            Gene weights score
    """
    
    # Compute raw score (overdispersion)
    if weight_opt=='od':
        return _compute_overdispersion_score(adata, gene_list, gene_weight)
    
    # Compute raw score (weighted average)
    if weight_opt=='uniform':
        v_score_weight = np.ones(len(gene_list))
    if weight_opt=='vs':
        v_score_weight = 1 / np.sqrt(adata.var.loc[gene_list,'var_tech'].values + 1e-2)
    if weight_opt=='inv_std':
        v_score_weight = 1 / np.sqrt(adata.var.loc[gene_list,'var'].values + 1e-2)
        
    if gene_weight is not None:
        v_score_weight = v_score_weight*np.array(gene_weight)
    v_score_weight = v_score_weight/v_score_weight.sum()
    v_raw_score = adata[:, gene_list].X.dot(v_score_weight).reshape([-1])  
        
    return v_raw_score,v_score_weight


def _compute_overdispersion_score(adata, gene_list, gene_weight):
    """Compute overdispersion score
    
        Let w_g_raw = gene_weight / \sigma_{tech,g}^2
        Let w_g = w_g_raw / \sum_g w_g_raw
        s_c = \sum_g w_g * [(X_cg - \mu_g)^2 - \sigma_{tech,g}^2]  
    Args
    ----
        adata (n_cell, n_gene) : AnnData
            adata.X should contain size-normalized log1p transformed count data
        gene_list (n_trait_gene) : list
            Trait gene list
        gene_weight (n_trait_gene) : list/np.ndarray
            Gene weights for genes in the gene_list 
            
    Returns
    -------
        v_raw_score (n_cell,) : np.ndarray
            Raw score
        v_score_weight (n_trait_gene,) : np.ndarray
            Gene weights score
    """
    
    v_mean = adata.var.loc[gene_list,'mean'].values
    v_var_tech = adata.var.loc[gene_list,'var_tech'].values
    
    v_w = 1 / (v_var_tech + 1e-2)
    if gene_weight is not None:
        v_w = v_w*np.array(gene_weight)
    v_w = v_w / v_w.sum()
    
    # Compute overdispersion score 
    if sp.sparse.issparse(adata.X):
        v_raw_score = adata[:, gene_list].X.power(2).dot(v_w).reshape([-1]) # Quadratic term 
    else:
        v_raw_score = (adata[:, gene_list].X**2).dot(v_w).reshape([-1]) # Quadratic term 
    v_raw_score = v_raw_score - adata[:, gene_list].X.dot(2*v_w*v_mean).reshape([-1]) # Linear term
    v_raw_score = v_raw_score + (v_w * (v_mean**2-v_var_tech)).sum() # Constant term
    
    return v_raw_score,np.ones(len(gene_list))


def _correct_background(v_raw_score, 
                        mat_ctrl_raw_score, 
                        v_var_ratio_c2t, 
                        save_intermediate=None):
    """Cell-wise and gene-wise background correction
    
    Args
    ----
        v_raw_score (n_cell,n_ctrl) : np.ndarray
            Trait raw score
        mat_ctrl_raw_score (n_cell,n_ctrl) : np.ndarray
            Control raw scores 
        v_var_ratio_c2t (n_cell) : np.ndarray
            Variance ratio between control scores and disease score
        save_intermediate : str
            File path prefix for saving intermediate results 
    Returns
    -------
        v_norm_score (n_cell,n_ctrl) : np.ndarray
            Trait normalized score
        mat_ctrl_norm_score (n_cell,n_ctrl) : np.ndarray
            Control normalized scores  
    """
    
    if save_intermediate is not None:
        np.savetxt(save_intermediate+'.raw_score.tsv.gz', 
                   v_raw_score, fmt='%.9e', delimiter='\t')
        np.savetxt(save_intermediate+'.ctrl_raw_score.tsv.gz', 
                   mat_ctrl_raw_score, fmt='%.9e', delimiter='\t')
    
    # Calibrate gene-sets (mean 0 and same ind. var)
    ind_zero_score = (v_raw_score==0)
    ind_zero_ctrl_score = (mat_ctrl_raw_score==0)
    
    v_raw_score = v_raw_score - v_raw_score.mean()
    mat_ctrl_raw_score = mat_ctrl_raw_score - mat_ctrl_raw_score.mean(axis=0)
    mat_ctrl_raw_score = mat_ctrl_raw_score/np.sqrt(v_var_ratio_c2t)
    if save_intermediate is not None:
        np.savetxt(save_intermediate+'.raw_score.1st_gs_alignment.tsv.gz', 
                   v_raw_score, fmt='%.9e', delimiter='\t')
        np.savetxt(save_intermediate+'.ctrl_raw_score.1st_gs_alignment.tsv.gz', 
                   mat_ctrl_raw_score, fmt='%.9e', delimiter='\t')
    
    # Cell-wise correction
    v_mean = mat_ctrl_raw_score.mean(axis=1)
    v_std = mat_ctrl_raw_score.std(axis=1)
    v_norm_score = v_raw_score.copy()
    v_norm_score = (v_norm_score - v_mean) / v_std
    mat_ctrl_norm_score = ((mat_ctrl_raw_score.T - v_mean) / v_std).T
    if save_intermediate is not None:
        np.savetxt(save_intermediate+'.raw_score.cellwise_standardization.tsv.gz', 
                   v_norm_score, fmt='%.9e', delimiter='\t')
        np.savetxt(save_intermediate+'.ctrl_raw_score.cellwise_standardization.tsv.gz', 
                   mat_ctrl_norm_score, fmt='%.9e', delimiter='\t')
            
    # Gene-set-wise correction
    v_norm_score = v_norm_score - v_norm_score.mean()
    mat_ctrl_norm_score = mat_ctrl_norm_score - mat_ctrl_norm_score.mean(axis=0)
    if save_intermediate is not None:
        np.savetxt(save_intermediate+'.raw_score.2nd_gs_alignment.tsv.gz', 
                   v_norm_score, fmt='%.9e', delimiter='\t')
        np.savetxt(save_intermediate+'.ctrl_raw_score.2nd_gs_alignment.tsv.gz', 
                   mat_ctrl_norm_score, fmt='%.9e', delimiter='\t')
    
    # Set cells with raw_score=0 to the minimum norm_score value
    norm_score_min = min(v_norm_score.min(), mat_ctrl_norm_score.min())
    v_norm_score[ind_zero_score] = norm_score_min-1e-3
    mat_ctrl_norm_score[ind_zero_ctrl_score] = norm_score_min
    if save_intermediate is not None:
        np.savetxt(save_intermediate+'.raw_score.final.tsv.gz', 
                   v_norm_score, fmt='%.9e', delimiter='\t')
        np.savetxt(save_intermediate+'.ctrl_raw_score.final.tsv.gz', 
                   mat_ctrl_norm_score, fmt='%.9e', delimiter='\t')
    
    return v_norm_score, mat_ctrl_norm_score


def compute_stats(data, n_mean_bin=20, n_var_bin=20, copy=False):
    """
    Precompute mean for each gene and mean&var for each cell
    """
    # Gene-wise statistics
    adata = data.copy() if copy else data
    adata.var['mean'],adata.var['var'] = _get_sparse_var(adata.X, axis=0)
    
    # Get the mean and var for the size-factor-normalized counts
    # It is highly correlated to the non-size-factor-normalized counts
    if sp.sparse.issparse(adata.X): # sp sparse matrix
        temp_X = adata.X.copy().expm1() # exp(X)-1 to get ct matrix from logct
    else:
        temp_X = np.expm1(adata.X) # numpy ndarray
    adata.var['ct_mean'],adata.var['ct_var'] = _get_sparse_var(temp_X, axis=0)
    del temp_X
    
    # Borrowed from scanpy _highly_variable_genes_seurat_v3
    not_const = adata.var['ct_var'].values>0
    estimat_var = np.zeros(adata.shape[1], dtype=np.float64)
    y = np.log10(adata.var['ct_var'].values[not_const])
    x = np.log10(adata.var['ct_mean'].values[not_const])
    model = loess(x, y, span=0.3, degree=2)
    model.fit()
    estimat_var[not_const] = model.outputs.fitted_values
    adata.var['ct_var_tech'] = 10**estimat_var
    # Recipe from Frost Nucleic Acids Research 2020
    adata.var['var_tech'] = adata.var['var']*adata.var['ct_var_tech']/adata.var['ct_var']
    adata.var.loc[adata.var['var_tech'].isna(),'var_tech'] = 0
    
    # Add n_mean_bin*n_var_bin mean_var bins
    n_bin_max = np.floor(np.sqrt(adata.shape[1]/10)).astype(int)
    if (n_mean_bin>n_bin_max) | (n_var_bin>n_bin_max):
        n_mean_bin,n_var_bin = n_bin_max,n_bin_max
        print('Too few genes for 20*20 bins, setting n_mean_bin=n_var_bin=%d'%(n_bin_max))
    v_mean_bin = pd.qcut(adata.var['mean'], n_mean_bin, labels=False, duplicates='drop')
    adata.var['mean_var'] = ''
    for bin_ in set(v_mean_bin):
        ind_select = (v_mean_bin==bin_)
        v_var_bin = pd.qcut(adata.var.loc[ind_select, 'var'], n_var_bin, labels=False, duplicates='drop')
        adata.var.loc[ind_select, 'mean_var'] = ['%s.%s'%(x,y) for x,y in zip(v_mean_bin[ind_select],v_var_bin)]
    
    # Cell-wise statistics
    adata.obs['mean'],adata.obs['var'] = _get_sparse_var(adata.X, axis=1)
    
    return adata if copy else None


def _get_sparse_var(sparse_X, axis=0):
    """
    Compute mean and var of a sparse matrix. 
    """   
    
    if sp.sparse.issparse(sparse_X):
        v_mean = sparse_X.mean(axis=axis)
        v_mean = np.array(v_mean).reshape([-1])
        v_var = sparse_X.power(2).mean(axis=axis)
        v_var = np.array(v_var).reshape([-1])
        v_var = v_var - v_mean**2
    else:
        v_mean = np.mean(sparse_X, axis=axis)
        v_var = np.var(sparse_X, axis=axis)
    
    return v_mean,v_var


def _get_p_from_empi_null(v_t,v_t_null):
    """Compute p-value from empirical null
    For score T and a set of null score T_1,...T_N, the p-value is 
        
        p= [1 + \Sigma_{i=1}^N 1_{ (T_i \geq T) }] / (1+N)
        
    If T, T1, ..., T_N are i.i.d. variables following a null distritbuion, 
    then p is super-uniform. 
    
    The naive algorithm is N^2. Here we provide an O(N log N) algorithm to 
    compute the p-value for each of the N elements in v_t
    
    Args
    ----
        v_t (M,): np.ndarray
            The observed score.
        v_t_null (N,): np.ndarray
            The null score.

    Returns
    -------
        v_p: (M,): np.ndarray
            P-value for each element in v_t
    """
    
    v_t = np.array(v_t)
    v_t_null = np.array(v_t_null)
    
    v_t_null = np.sort(v_t_null)    
    v_pos = np.searchsorted(v_t_null, v_t, side='left')
    v_p = (v_t_null.shape[0]-v_pos+1)/(v_t_null.shape[0]+1)
    return v_p


##############################################################################
######################### Code for comparison methods ########################
##############################################################################

def score_cell_vision(adata, gene_list):
    
    """Score cells based on the trait gene set 
    
    Args
    ----
        data (n_cell, n_gene) : AnnData
            data.X should contain size-normalized log1p transformed count data
        gene_list (n_trait_gene) : list
            Trait gene list  

    Returns
    -------
        df_res (n_cell, n_key) : pd.DataFrame (dtype=np.float32)
            Columns: 
            1. score: Vision signature score
            2. pval: p-value computed from the Vision score
    """
    
    gene_list = sorted(set(gene_list) & set(adata.var_names))
    v_mean,v_var = _get_sparse_var(adata.X, axis=1)
    
    v_score = adata[:, gene_list].X.mean(axis=1)
    v_score = np.array(v_score).reshape([-1])
    v_score = (v_score - v_mean) / np.sqrt(v_var/len(gene_list))
    v_p = 1-sp.stats.norm.cdf(v_score)
    
    v_norm_score = (v_score - v_score.mean())/v_score.std()
    v_norm_p = 1-sp.stats.norm.cdf(v_norm_score)
    
    dic_res = {'score':v_score, 'pval':v_p, 'norm_score':v_norm_score, 'norm_pval':v_norm_p}
    df_res = pd.DataFrame(index=adata.obs.index,  data=dic_res, dtype=np.float32)
    return df_res


def score_cell_scanpy(adata, gene_list):
    
    """Score cells based on the trait gene set 
    
    Args
    ----
        data (n_cell, n_gene) : AnnData
            data.X should contain size-normalized log1p transformed count data
        gene_list (n_trait_gene) : list
            Trait gene list  

    Returns
    -------
        df_res (n_cell, n_key) : pd.DataFrame (dtype=np.float32)
            Columns: 
            1. score: Vision signature score
            2. pval: p-value computed from the Vision score
    """
    
    gene_list = sorted(set(gene_list) & set(adata.var_names))
    sc.tl.score_genes(adata, gene_list=gene_list)
    
    v_score = adata.obs['score']
    v_score_z = (v_score - v_score.mean()) / np.sqrt(v_score.var())
    v_p = 1-sp.stats.norm.cdf(v_score_z)
    
    dic_res = {'score':v_score, 'score_z':v_score_z, 'pval':v_p}
    df_res = pd.DataFrame(index=adata.obs.index,  data=dic_res, dtype=np.float32)
    return df_res

##############################################################################
######################## Code for downstream analysis ########################
##############################################################################

def reg_out(mat_Y, mat_X):
    """Regress mat_X out of mat_Y
    
    Args
    ----
        mat_Y (n_sample, n_response) : np.ndarray
            Response variable
        mat_X (n_sample, n_covariates) : np.ndarray
            Covariates
            
    Returns
    -------
        mat_Y_resid (n_sample, n_response) : np.ndarray
            Response variable residual
    """
    
    mat_X = np.array(mat_X)
    if len(mat_X.shape)==1:
        mat_X = mat_X.reshape([-1,1]) 
    mat_Y = np.array(mat_Y)
    if len(mat_Y.shape)==1:
        mat_Y = mat_Y.reshape([-1,1]) 
    
    n_sample = mat_Y.shape[0]
    mat_xtx = np.dot(mat_X.T, mat_X)/n_sample
    mat_xty = np.dot(mat_X.T, mat_Y)/n_sample
    mat_coef = np.linalg.solve(mat_xtx, mat_xty)
    mat_Y_resid = mat_Y - mat_X.dot(mat_coef)
    
    if mat_Y_resid.shape[1]==1:
        mat_Y_resid = mat_Y_resid.reshape([-1])
    
    return mat_Y_resid


def correlate_gene(data,
                   trs_name='trs_ez',
                   suffix='',
                   corr_opt='pearson',
                   cov_list=None,
                   copy=False):
    
    """Compute the correlation between gene expressions and TRS
    
    Args
    ----
        data (n_cell, n_gene) : AnnData
            adata.X should contain size-normalized log1p transformed count data
        trs_name : str
            The variable to correlate gene expression with. Should be one column in data.obs.
        suffix : str
            The name of the added gene-wise correlation would be 'trs_corr'+suffix.
        corr_opt : str
            Option for computing the correlation
            'pearson': Pearson's correlation
            'spearman': Spearman's correlation
        cov_list : list of str
            Covariates to control for.
            The covariates are first centered and then regressed out from 
                both trs_name and the gene expression before computing the correlation.
            Elements in cov_list should be present in data.obs.columns
        copy : bool
            If to make copy of the AnnData object
            
    Returns
    -------
        adata (AnnData): 
            Add the columns 'scdrs_corr'+suffix to data.var
    """
    
    adata = data.copy() if copy else data
    
    # Check options
    corr_opt_list = ['pearson', 'spearman']
    if corr_opt not in corr_opt_list:
        raise ValueError('# compute_scdrs_corr: corr_opt not in [%s]'
                         %', '.join([str(x) for x in corr_opt_list]))
    if trs_name not in adata.obs.columns:
        raise ValueError('# compute_scdrs_corr: %s not in data.obs.columns'%trs_name)
    if cov_list is not None:
        temp_list = list(set(cov_list) - set(adata.obs.columns))
        if len(temp_list)>0:
            raise ValueError('# compute_scdrs_corr: covariates %s not in data.obs.columns'
                             %','.join(temp_list))
    
    # Get data 
    mat_X = data.X.toarray()
    v_trs = data.obs[trs_name].values.copy()
    
    # Regress out covariates
    if cov_list is not None:
        mat_cov = adata.obs[cov_list].values.copy()        
        mat_cov = mat_cov - mat_cov.mean(axis=0)
        v_trs = reg_out(v_trs, mat_cov)
        mat_X = reg_out(mat_X, mat_cov)
    
    # Compute correlation
    if corr_opt=='pearson':
        v_corr = _pearson_corr(mat_X, v_trs)
               
    if corr_opt=='spearman':
        v_corr = _spearman_corr(mat_X, v_trs)
    
    adata.var['scdrs_corr'+suffix] = v_corr
        
    return adata if copy else None 
    
    
def _pearson_corr(mat_X, mat_Y):
    """Pearson's correlation between every columns in mat_X and mat_Y
    
    Args
    ----
        mat_X (N,M1): np.ndarray
        
        mat_Y (N,M2): np.ndarray

    Returns
    -------
        mat_corr: (M1,M2): np.ndarray
            Correlation matrix
    """
    
    # If sparse, use _pearson_corr_sparse
    if sp.sparse.issparse(mat_X) | sp.sparse.issparse(mat_Y):
        return _pearson_corr_sparse(mat_X, mat_Y)
    
    # Reshape 
    if len(mat_X.shape)==1:
        mat_X = mat_X.reshape([-1,1])
    if len(mat_Y.shape)==1:
        mat_Y = mat_Y.reshape([-1,1])
        
    mat_X = (mat_X-mat_X.mean(axis=0))/mat_X.std(axis=0).clip(min=1e-8)
    mat_Y = (mat_Y-mat_Y.mean(axis=0))/mat_Y.std(axis=0).clip(min=1e-8)
    mat_corr = mat_X.T.dot(mat_Y)/mat_X.shape[0]
    mat_corr = np.array(mat_corr, dtype=np.float32)
    
    if (mat_X.shape[1]==1) | (mat_Y.shape[1]==1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr
    

def _pearson_corr_sparse(mat_X, mat_Y):
    """Pearson's correlation between every columns in mat_X and mat_Y (sparse matrix)
    
    Args
    ----
        mat_X (N,M1): sp.sparse
        
        mat_Y (N,M2): sp.sparse

    Returns
    -------
        mat_corr: (M1,M2): np.ndarray
            Correlation matrix
    """

    # Reshape 
    if len(mat_X.shape)==1:
        mat_X = mat_X.reshape([-1,1])
    if len(mat_Y.shape)==1:
        mat_Y = mat_Y.reshape([-1,1])
        
    # Convert to sparse matrix if not already sparse 
    if sp.sparse.issparse(mat_X) is False:
        mat_X = sp.sparse.csr_matrix(mat_X)
    if sp.sparse.issparse(mat_Y) is False:
        mat_Y = sp.sparse.csr_matrix(mat_Y)
        
    # Compute v_mean,v_var
    v_X_mean,v_X_var = _get_sparse_var(mat_X, axis=0)
    v_X_sd = np.sqrt(v_X_var).clip(min=1e-8)
    v_Y_mean,v_Y_var = _get_sparse_var(mat_Y, axis=0)
    v_Y_sd = np.sqrt(v_Y_var).clip(min=1e-8)
    
    mat_corr = mat_X.T.dot(mat_Y)/mat_X.shape[0]
    mat_corr = mat_corr - v_X_mean.reshape([-1,1]).dot(v_Y_mean.reshape([1,-1]))
    mat_corr = mat_corr / v_X_sd.reshape([-1,1]).dot(v_Y_sd.reshape([1,-1]))
    mat_corr = np.array(mat_corr, dtype=np.float32)
    
    if (mat_X.shape[1]==1) | (mat_Y.shape[1]==1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr   
    
    
def _spearman_corr(mat_X, mat_Y):
    """Spearman's correlation between every columns in mat_X and mat_Y
    
    Args
    ----
        mat_X (N,M1): np.ndarray
        
        mat_Y (N,M2): np.ndarray

    Returns
    -------
        mat_corr (M1,M2): np.ndarray
            Correlation matrix
    """
    
    # Reshape 
    if len(mat_X.shape)==1:
        mat_X = mat_X.reshape([-1,1])
    if len(mat_Y.shape)==1:
        mat_Y = mat_Y.reshape([-1,1])
        
    mat_X = _get_rank(mat_X, axis=0)
    mat_Y = _get_rank(mat_Y, axis=0)

    mat_X = (mat_X-mat_X.mean(axis=0))/mat_X.std(axis=0).clip(min=1e-8)
    mat_Y = (mat_Y-mat_Y.mean(axis=0))/mat_Y.std(axis=0).clip(min=1e-8)
    mat_corr = mat_X.T.dot(mat_Y)/mat_X.shape[0]
    
    if (mat_X.shape[1]==1) | (mat_Y.shape[1]==1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr
    

def _get_rank(mat_X, axis=0):
    """Spearman's correlation between every columns in mat_X and mat_Y
    
    Args
    ----
        mat_X (N,M): np.ndarray
        axis: int
            axis=0: column-wise rank (across rows)
            axis=1: row-wise rank (across columns)
    Returns
    -------
        mat_rank  (N,M): np.ndarray
            Rank matrix
    """
    
    if axis==0:
        mat_X = np.argsort(mat_X, axis=0)
        mat_rank = np.empty_like(mat_X)
        temp_v = np.arange(mat_X.shape[0])
        for i_col in range(mat_X.shape[1]):
            mat_rank[mat_X[:,i_col], i_col] = temp_v
        
    if axis==1:
        mat_X = np.argsort(mat_X, axis=1)
        mat_rank = np.empty_like(mat_X)
        temp_v = np.arange(mat_X.shape[1])
        for i_row in range(mat_X.shape[0]):
            mat_rank[i_row, mat_X[i_row,:]] = temp_v
    
    return mat_rank
    
    
def compute_gene_contrib(data,
                         gene_list,
                         random_seed=0,
                         copy=False,
                         verbose=False):
    
    """Find the genes 
    
    Args
    ----
        data (n_cell, n_gene) : AnnData
            adata.X should contain size-normalized log1p transformed count data
        gene_list (n_trait_gene) : list
            Trait gene list
        copy : bool
            If to make copy of the AnnData object
            
    Returns
    -------
        adata (AnnData): 
            Add the columns 'trs_corr'+suffix to data.var
    """
    
    np.random.seed(random_seed)
    adata = data.copy() if copy else data
    
    # Pre-compute statistics
    if 'var_tech' not in adata.var.columns:
        if verbose: 
            print('# score_cell: recompute statistics using method.compute_stats')
        compute_stats(adata)
        
    # Compute contribution for each gene 
    gene_list = sorted(set(gene_list) & set(adata.var_names))
    v_score_weight = 1 / np.sqrt(adata.var.loc[gene_list,'var_tech'].values + 1e-2)
    df_contrib = pd.DataFrame(index=adata.obs_names, columns=gene_list,
                              data=adata[:,gene_list].X.toarray())
    df_contrib = df_contrib * v_score_weight
    
    return df_contrib

##############################################################################
################################## Old code ##################################
##############################################################################

# Can we delete this one?
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


# Can we delete this one?
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