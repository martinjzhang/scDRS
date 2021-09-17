import pytest 
import numpy as np
import scipy as sp
import pandas as pd
import os
from anndata import read_h5ad
import scdrs


"""
    score_cell
"""
def test_score_cell():
    
    # Load toy data 
    DATA_PATH=scdrs.__path__[0]
    H5AD_FILE=os.path.join(DATA_PATH,'data/toydata_mouse.h5ad')
    COV_FILE=os.path.join(DATA_PATH,'data/toydata_mouse.cov')
    GS_FILE=os.path.join(DATA_PATH,'data/toydata_mouse.gs')
    assert os.path.exists(H5AD_FILE), "built-in data toydata_mouse.h5ad missing"
    assert os.path.exists(COV_FILE), "built-in data toydata_mouse.cov missing"
    assert os.path.exists(GS_FILE), "built-in data toydata_mouse.gs missing"
    
    # Load built-in data 
    adata = read_h5ad(H5AD_FILE)
    
    df_cov = pd.read_csv(COV_FILE, sep='\t', index_col=0)
    cov_list = list(df_cov.columns)
    adata.obs.drop([x for x in cov_list if x in adata.obs.columns], axis=1, inplace=True)
    adata.obs = adata.obs.join(df_cov)
    adata.obs.fillna(adata.obs[cov_list].mean(), inplace=True)
    adata.var['mean'] = adata.X.mean(axis=0).T
    if sp.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X -= adata.var['mean'].values
    adata.X = scdrs.method.reg_out(adata.X, adata.obs[cov_list].values)
    adata.X += adata.var['mean']
    
    df_gs = pd.read_csv(GS_FILE, sep='\t')
    df_gs.index = df_gs['TRAIT']
    
    # Compute score 
    scdrs.method.compute_stats(adata)
    gene_list = df_gs.loc['toydata_gs_mouse','GENESET'].split(',')
    gene_list = sorted(set(gene_list) & set(adata.var_names))            
    df_res = scdrs.method.score_cell(adata, gene_list, ctrl_match_key='mean_var', n_ctrl=20, weight_opt='vs',
                                     return_ctrl_raw_score=False, return_ctrl_norm_score=False, verbose=False)
    # Check raw_score
    v_val = df_res['raw_score'].values[[0,4,9]]
    v_ref = np.array([3.4417892, 3.3382401, 3.2723658])
    err_msg = 'raw_score of cells 0,4,9 does not match reference: actual=%s expect=%s'%\
        (','.join(['%0.4f'%x for x in v_val]), ','.join(['%0.4f'%x for x in v_ref]))
    assert np.absolute(v_val-v_ref).mean()<1e-3, err_msg
    
    # Check raw_score
    v_val = df_res['raw_score'].values[[0,4,9]]
    v_ref = np.array([3.4417892, 3.3382401, 3.2723658])
    err_msg = 'raw_score of cells 0,4,9 does not match reference: actual=%s expect=%s'%\
        (','.join(['%0.4f'%x for x in v_val]), ','.join(['%0.4f'%x for x in v_ref]))
    assert np.absolute(v_val-v_ref).mean()<1e-3, err_msg
    
    # Check norm_score
    v_val = df_res['norm_score'].values[[0,4,9]]
    v_ref = np.array([7.084109, 6.1206293, -2.7951033])
    err_msg = 'norm_score of cells 0,4,9 does not match reference: actual=%s expect=%s'%\
        (','.join(['%0.4f'%x for x in v_val]), ','.join(['%0.4f'%x for x in v_ref]))
    assert np.absolute(v_val-v_ref).mean()<1e-3, err_msg
        
    return