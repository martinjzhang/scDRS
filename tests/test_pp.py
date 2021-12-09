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
    
    mat_Y = np.array()
    mat_X = np.array()
    
    
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
