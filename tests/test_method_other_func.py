import pytest 
import numpy as np
import scipy as sp

import scdrs


"""
    _select_ctrl_geneset
"""
def test_select_ctrl_geneset():
    return


"""
    _compute_raw_score
"""
def test_compute_raw_score():
    return


"""
    _correct_background
"""
def test_correct_background():
    return


"""
    compute_stats
"""
def test_compute_stats():
    return


"""
    _get_sparse_var
"""
def test_get_sparse_var():
    return


"""
    _get_p_from_empi_null
"""
def test_get_p_from_empi_null():
    v_t = [0,1]
    v_t_null = [0.5, 0.6, 0.7]
    v_p = scdrs.method._get_p_from_empi_null(v_t,v_t_null)
    assert np.absolute(v_p[0]-1)<0.001, "First MC-p should be 1"
    assert np.absolute(v_p[1]-0.25)<0.001, "Second MC-p should be 1/(3+1)=0.25"
    return


"""
    reg_out
"""
def test_reg_out():
    return


"""
    correlate_gene
"""
def test_correlate_gene():
    return


"""
    _pearson_corr
"""
def test_pearson_corr():
    return


"""
    _spearman_corr
"""
def test_spearman_corr():
    return


"""
    _get_rank
"""
def test_get_rank():
    return


"""
    compute_gene_contrib
"""
def test_compute_gene_contrib():
    return
