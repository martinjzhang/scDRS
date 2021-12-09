import pytest
import numpy as np
import scipy as sp
from scipy import sparse

import scdrs


def test_preprocess():
    """
    Test scdrs.pp.preprocess
    """
    return


def test_compute_stats():
    """
    Test scdrs.pp.compute_stats
    """
    return


def test_reg_out():
    """
    Test scdrs.pp.reg_out
    """
    
    mat_list = [np.array([1, 2, 3]), 
                np.array([[1, 2], [3, 4], [5, 6]]), 
                sparse.csr_matrix([1, 2, 3]).T,
                sparse.csr_matrix([[1, 2], [3, 4], [5, 6]])]
    for mat_X in mat_list:
        for mat_Y in mat_list:
            mat_Y_resid = scdrs.pp.reg_out(mat_Y, mat_X)
            
            # Compare results with np.linalg.lstsq
            if sparse.issparse(mat_X):
                mat_X = mat_X.toarray()
            if len(mat_X.shape)==1:
                mat_X = mat_X.reshape([-1,1])
            if sparse.issparse(mat_Y):
                mat_Y = mat_Y.toarray()
            if len(mat_Y.shape)==1:
                mat_Y = mat_Y.reshape([-1,1])
            mat_beta,_,_,_ = np.linalg.lstsq(mat_X, mat_Y, rcond=-1)
            mat_Y_resid_true = mat_Y - mat_X.dot(mat_beta)
            
            print('mat_X', mat_X)
            print('mat_Y', mat_Y)
            err_msg = 'avg_abs_dif=%0.2e'%np.absolute(mat_Y_resid- mat_Y_resid_true).mean()
            assert np.allclose(
                        mat_Y_resid, mat_Y_resid_true, rtol=1e-5, equal_nan=True
                    ), err_msg

    return


def test_get_mean_var():
    """
    Test scdrs.pp._get_mean_var
    """

    mat_test = np.array([[1, 2], [3, 4]], dtype=float)

    for opt in ["dense", "sparse"]:
        for axis_ in [0, 1]:
            if opt == "sparse":
                mat_ = sparse.csr_matrix(mat_test)
            else:
                mat_ = mat_test

            v_mean, v_var = scdrs.pp._get_mean_var(mat_, axis=axis_)

            v_mean_true = np.mean(mat_test, axis=axis_)
            v_var_true = np.var(mat_test, axis=axis_)

            err_msg = (
                "mode=%s, axis=%d, avg_abs_mean_dif=%0.2e, avg_abs_var_dif=%0.2e"
                % (
                    opt,
                    axis_,
                    np.absolute(v_mean - v_mean_true).mean(),
                    np.absolute(v_var - v_var_true).mean(),
                )
            )
            assert np.allclose(
                v_mean, v_mean_true, rtol=1e-5, equal_nan=True
            ) & np.allclose(v_var, v_var_true, rtol=1e-5, equal_nan=True), err_msg

    return


def test_get_mean_var_implicit_cov_corr():
    """
    Test scdrs.pp._get_mean_var_implicit_cov_corr
    """
    return
