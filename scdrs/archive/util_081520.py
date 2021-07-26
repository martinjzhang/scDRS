import numpy as np
import scipy as sp
import pandas as pd

def test_sctrs():
    print('# test_sctrs')
    return 

def compute_div(v_x, v_y, opt='kl', n_bin=100):
    """compute divergence between two sets of samples 
    
    Args:
        v_x (list or np.adarray): sample 1
        v_y (list or np.adarray): sample 2
        opt (str): divergence measure, one of ['kl', 'l1']
        n_bin=100 (int): number of bins

    Returns:
        div: divergence between two sets of samples 
    """
    
    if opt not in ['kl', 'l1']:
        raise ValueError('# compute_divergence: %s not in [kl, l1]'%opt)
        
    v_x = np.array(v_x, dtype=float)
    v_y = np.array(v_y, dtype=float)
    
    min_ = min(v_x.min(), v_y.min())
    max_ = max(v_x.max(), v_y.max())
    bins_ = np.linspace(min_, max_, n_bin+1)
    
    p_x = np.histogram(v_x, bins=bins_)[0]+1/(bins_.shape[0]-1)
    p_y = np.histogram(v_y, bins=bins_)[0]+1/(bins_.shape[0]-1)
    p_x = p_x/p_x.sum()
    p_y = p_y/p_y.sum()
    
    if opt=='kl':
        div = (p_x*np.log(p_x/p_y)).sum()
        
    if opt=='l1':
        div = np.absolute(p_x-p_y).sum()
    
    return div 