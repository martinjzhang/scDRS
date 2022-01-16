import pytest
import numpy as np
import scipy as sp
import pandas as pd
from scipy import sparse
from .test_method_score_cell_main import load_toy_data
import scdrs


def test_preprocess_sparse_cov():
    """
    Test scdrs.pp.preprocess: sparse + cov (implicit-covariate-correction mode)
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    scdrs.pp.preprocess(adata, cov=df_cov)

    assert "SCDRS_PARAM" in adata.uns, "'SCDRS_PARAM' should be in adata.uns"
    KEY_LIST = [
        "FLAG_SPARSE",
        "FLAG_COV",
        "COV_MAT",
        "COV_BETA",
        "COV_GENE_MEAN",
        "GENE_STATS",
        "CELL_STATS",
    ]
    for key_ in KEY_LIST:
        assert key_ in adata.uns["SCDRS_PARAM"], (
            "'%s' should be in adata.uns['SCDRS_PARAM']" % key_
        )
    assert adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"] is True, "FLAG_SPARSE=True"
    assert adata.uns["SCDRS_PARAM"]["FLAG_COV"] is True, "FLAG_COV=True"

    return


def test_preprocess_sparse_nocov():
    """
    Test scdrs.pp.preprocess: sparse + nocov (normal mode)
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    scdrs.pp.preprocess(adata, cov=None)

    assert "SCDRS_PARAM" in adata.uns, "'SCDRS_PARAM' should be in adata.uns"
    KEY_LIST = [
        "FLAG_SPARSE",
        "FLAG_COV",
        "GENE_STATS",
        "CELL_STATS",
    ]
    for key_ in KEY_LIST:
        assert key_ in adata.uns["SCDRS_PARAM"], (
            "'%s' should be in adata.uns['SCDRS_PARAM']" % key_
        )
    assert adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"] is True, "FLAG_SPARSE=True"
    assert adata.uns["SCDRS_PARAM"]["FLAG_COV"] is False, "FLAG_COV=False"

    return


def test_preprocess_dense_cov():
    """
    Test scdrs.pp.preprocess: dense + cov (normal mode)
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    adata.X = adata.X.toarray()
    scdrs.pp.preprocess(adata, cov=df_cov)

    assert "SCDRS_PARAM" in adata.uns, "'SCDRS_PARAM' should be in adata.uns"
    KEY_LIST = [
        "FLAG_SPARSE",
        "FLAG_COV",
        "GENE_STATS",
        "CELL_STATS",
    ]
    for key_ in KEY_LIST:
        assert key_ in adata.uns["SCDRS_PARAM"], (
            "'%s' should be in adata.uns['SCDRS_PARAM']" % key_
        )
    assert adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"] is False, "FLAG_SPARSE=False"
    assert adata.uns["SCDRS_PARAM"]["FLAG_COV"] is True, "FLAG_COV=True"

    return


def test_preprocess_dense_nocov():
    """
    Test scdrs.pp.preprocess: dense + nocov (normal mode)
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    adata.X = adata.X.toarray()
    scdrs.pp.preprocess(adata, cov=None)

    assert "SCDRS_PARAM" in adata.uns, "'SCDRS_PARAM' should be in adata.uns"
    KEY_LIST = [
        "FLAG_SPARSE",
        "FLAG_COV",
        "GENE_STATS",
        "CELL_STATS",
    ]
    for key_ in KEY_LIST:
        assert key_ in adata.uns["SCDRS_PARAM"], (
            "'%s' should be in adata.uns['SCDRS_PARAM']" % key_
        )
    assert adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"] is False, "FLAG_SPARSE=False"
    assert adata.uns["SCDRS_PARAM"]["FLAG_COV"] is False, "FLAG_COV=False"

    return


def test_preprocess_consistency():
    """
    Test scdrs.pp.preprocess: consistency between sparse+cov and dense+cov
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata_sparse = adata[:, :50].copy()
    scdrs.pp.preprocess(adata_sparse, cov=df_cov)

    adata_dense = adata[:, :50].copy()
    adata_dense.X = adata_dense.X.toarray()
    scdrs.pp.preprocess(adata_dense, cov=df_cov)

    # COV
    mat_X_sparse = (
        adata_sparse.X.toarray()
        + adata_sparse.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
            adata_sparse.uns["SCDRS_PARAM"]["COV_BETA"].values.T
        )
        + adata_sparse.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].values
    )
    err_msg = "Covariate-corrected matrices should be same for sparse+cov and dense+cov"
    assert np.allclose(
        mat_X_sparse,
        adata_dense.X,
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    # GENE_STATS
    err_msg = "'GENE_STATS' should be same for sparse+cov and dense+cov"
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        adata_dense.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    # CELL_STATS
    err_msg = "'CELL_STATS' should be same for sparse+cov and dense+cov"
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        adata_dense.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    return


def test_preprocess_consistency_nocov():
    """
    Test scdrs.pp.preprocess: consistency between sparse+nocov and dense+nocov
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata_sparse = adata[:, :50].copy()
    scdrs.pp.preprocess(adata_sparse, cov=None)

    adata_dense = adata[:, :50].copy()
    adata_dense.X = adata_dense.X.toarray()
    scdrs.pp.preprocess(adata_dense, cov=None)

    # GENE_STATS
    err_msg = "'GENE_STATS' should be same for sparse+nocov and dense+nocov"
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        adata_dense.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    # CELL_STATS
    err_msg = "'CELL_STATS' should be same for sparse+nocov and dense+nocov"
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        adata_dense.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    return


def test_compute_stats():
    """
    Test scdrs.pp.compute_stats: normal mode

    TODO:
    -----
    1. The following columns are not tested: df_gene["var_tech"],
    df_gene["ct_var_tech"], df_gene["mean_var"]
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    df_gene, df_cell = scdrs.pp.compute_stats(adata)

    mat_X = adata.X.toarray()
    mat_X_ct = np.expm1(mat_X)

    # df_gene["mean"]
    v_mean_true = np.mean(mat_X, axis=0)
    err_msg = (
        "avg_abs_gene_mean_dif=%0.2e"
        % np.absolute(df_gene["mean"] - v_mean_true).mean()
    )
    assert np.allclose(df_gene["mean"], v_mean_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["var"]
    v_var_true = np.var(mat_X, axis=0)
    err_msg = (
        "avg_abs_gene_var_dif=%0.2e" % np.absolute(df_gene["var"] - v_var_true).mean()
    )
    assert np.allclose(df_gene["var"], v_var_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["ct_mean"]
    v_mean_true = np.mean(mat_X_ct, axis=0)
    err_msg = (
        "avg_abs_gene_ct_mean_dif=%0.2e"
        % np.absolute(df_gene["ct_mean"] - v_mean_true).mean()
    )
    assert np.allclose(
        df_gene["ct_mean"], v_mean_true, rtol=1e-4, equal_nan=True
    ), err_msg

    # df_gene["ct_var"]
    v_var_true = np.var(mat_X_ct, axis=0)
    err_msg = (
        "avg_abs_gene_ct_var_dif=%0.2e"
        % np.absolute(df_gene["ct_var"] - v_var_true).mean()
    )
    assert np.allclose(
        df_gene["ct_var"], v_var_true, rtol=1e-4, equal_nan=True
    ), err_msg

    # df_cell["mean"]
    v_mean_true = np.mean(mat_X, axis=1)
    err_msg = (
        "avg_abs_cell_mean_dif=%0.2e"
        % np.absolute(df_cell["mean"] - v_mean_true).mean()
    )
    assert np.allclose(df_cell["mean"], v_mean_true, rtol=1e-4, equal_nan=True), err_msg

    # df_cell["var"]
    v_var_true = np.var(mat_X, axis=1)
    err_msg = (
        "avg_abs_cell_var_dif=%0.2e" % np.absolute(df_cell["var"] - v_var_true).mean()
    )
    assert np.allclose(df_cell["var"], v_var_true, rtol=1e-4, equal_nan=True), err_msg

    return


def test_compute_stats_implicit_covariate_correction():
    """
    Test scdrs.pp.compute_stats: implicit covariate correction mode

    TODO:
    -----
    1. The following columns are not tested: df_gene["var_tech"],
    df_gene["ct_var_tech"], df_gene["mean_var"]
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    scdrs.pp.preprocess(adata, cov=df_cov)
    df_gene, df_cell = scdrs.pp.compute_stats(adata, implicit_cov_corr=True)

    # adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN
    mat_X = adata.X.toarray()
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
        adata.uns["SCDRS_PARAM"]["COV_BETA"].values.T
    )
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].values
    mat_X_ct = np.expm1(mat_X)

    # df_gene["mean"]
    v_mean_true = np.mean(mat_X, axis=0)
    err_msg = (
        "avg_abs_gene_mean_dif=%0.2e"
        % np.absolute(df_gene["mean"] - v_mean_true).mean()
    )
    assert np.allclose(df_gene["mean"], v_mean_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["var"]
    v_var_true = np.var(mat_X, axis=0)
    err_msg = (
        "avg_abs_gene_var_dif=%0.2e" % np.absolute(df_gene["var"] - v_var_true).mean()
    )
    assert np.allclose(df_gene["var"], v_var_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["ct_mean"]
    v_mean_true = np.mean(mat_X_ct, axis=0)
    err_msg = (
        "avg_abs_gene_ct_mean_dif=%0.2e"
        % np.absolute(df_gene["ct_mean"] - v_mean_true).mean()
    )
    assert np.allclose(
        df_gene["ct_mean"], v_mean_true, rtol=1e-4, equal_nan=True
    ), err_msg

    # df_gene["ct_var"]
    v_var_true = np.var(mat_X_ct, axis=0)
    err_msg = (
        "avg_abs_gene_ct_var_dif=%0.2e"
        % np.absolute(df_gene["ct_var"] - v_var_true).mean()
    )
    assert np.allclose(
        df_gene["ct_var"], v_var_true, rtol=1e-4, equal_nan=True
    ), err_msg

    # df_cell["mean"]
    v_mean_true = np.mean(mat_X, axis=1)
    err_msg = (
        "avg_abs_cell_mean_dif=%0.2e"
        % np.absolute(df_cell["mean"] - v_mean_true).mean()
    )
    assert np.allclose(df_cell["mean"], v_mean_true, rtol=1e-4, equal_nan=True), err_msg

    # df_cell["var"]
    v_var_true = np.var(mat_X, axis=1)
    err_msg = (
        "avg_abs_cell_var_dif=%0.2e" % np.absolute(df_cell["var"] - v_var_true).mean()
    )
    assert np.allclose(df_cell["var"], v_var_true, rtol=1e-4, equal_nan=True), err_msg

    return


def test_reg_out():
    """
    Test scdrs.pp.reg_out
    """

    mat_list = [
        np.array([1, 2, 3]),
        np.array([[1, 2], [3, 4], [5, 6]]),
        sparse.csr_matrix([1, 2, 3]).T,
        sparse.csr_matrix([[1, 2], [3, 4], [5, 6]]),
    ]
    for mat_X in mat_list:
        for mat_Y in mat_list:
            mat_Y_resid = scdrs.pp.reg_out(mat_Y, mat_X)

            # Compare results with np.linalg.lstsq
            if sparse.issparse(mat_X):
                mat_X = mat_X.toarray()
            if len(mat_X.shape) == 1:
                mat_X = mat_X.reshape([-1, 1])
            if sparse.issparse(mat_Y):
                mat_Y = mat_Y.toarray()
            if len(mat_Y.shape) == 1:
                mat_Y = mat_Y.reshape([-1, 1])
            mat_beta, _, _, _ = np.linalg.lstsq(mat_X, mat_Y, rcond=-1)
            mat_Y_resid_true = mat_Y - mat_X.dot(mat_beta)

            print("mat_X")
            print(mat_X)
            print("mat_Y")
            print(mat_Y)
            err_msg = (
                "avg_abs_dif=%0.2e" % np.absolute(mat_Y_resid - mat_Y_resid_true).mean()
            )
            assert np.allclose(
                mat_Y_resid, mat_Y_resid_true, rtol=1e-4, equal_nan=True
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
                v_mean, v_mean_true, rtol=1e-4, equal_nan=True
            ) & np.allclose(v_var, v_var_true, rtol=1e-4, equal_nan=True), err_msg

    return


def test_get_mean_var_implicit_cov_corr():
    """
    Test scdrs.pp._get_mean_var_implicit_cov_corr
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    np.random.seed(0)
    adata.uns["SCDRS_PARAM"] = {
        "COV_MAT": df_cov,
        "COV_BETA": pd.DataFrame(
            index=adata.var_names,
            columns=df_cov.columns,
            data=np.random.randn(adata.shape[1], df_cov.shape[1]),
        ),
        "COV_GENE_MEAN": pd.Series(
            data=np.random.randn(adata.shape[1]), index=adata.var_names
        ),
    }

    # adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN
    mat_X = adata.X.toarray()
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
        adata.uns["SCDRS_PARAM"]["COV_BETA"].values.T
    )
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].values

    for axis_ in [0, 1]:
        v_mean, v_var = scdrs.pp._get_mean_var_implicit_cov_corr(adata, axis=axis_)
        v_mean_true = np.mean(mat_X, axis=axis_)
        v_var_true = np.var(mat_X, axis=axis_)
        err_msg = "axis=%d, avg_abs_mean_dif=%0.2e, avg_abs_var_dif=%0.2e" % (
            axis_,
            np.absolute(v_mean - v_mean_true).mean(),
            np.absolute(v_var - v_var_true).mean(),
        )
        assert np.allclose(
            v_mean, v_mean_true, rtol=1e-4, equal_nan=True
        ) & np.allclose(v_var, v_var_true, rtol=1e-4, equal_nan=True), err_msg

    return