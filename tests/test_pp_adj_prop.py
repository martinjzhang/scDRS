import pytest
import numpy as np
import scipy as sp
import pandas as pd
from scipy import sparse
from .test_method_score_cell_main import load_toy_data
from anndata import AnnData
import scdrs


def load_toy_data_adj_prop():
    """
    Toy data for checking `adj_prop`
    """
    np.random.seed(0)
    n_cell, n_cell_ct1, n_gene = 1000, 100, 100

    df_obs = pd.DataFrame(
        index=["cell%d" % x for x in range(n_cell)],
        data={
            "cell_group": ["ct1"] * n_cell_ct1 + ["ct2"] * (n_cell - n_cell_ct1),
            "cov": np.random.randn(n_cell),
        },
    )
    mat_X = np.random.rand(n_cell, n_gene)
    mat_X[df_obs["cell_group"] == "ct1", :50] = (
        mat_X[df_obs["cell_group"] == "ct1", :50] + 0.5
    ) * 2
    mat_X[df_obs["cell_group"] == "ct2", 50:] = (
        mat_X[df_obs["cell_group"] == "ct2", 50:] + 0.5
    ) * 2

    adata = AnnData(X=sp.sparse.csr_matrix(mat_X), obs=df_obs)
    adata_balance = adata[: 2 * n_cell_ct1].copy()
    df_cov = df_obs[["cov"]].copy()

    return adata, adata_balance, df_cov


def test_preprocess_consistency_adj_prop():
    """
    Test scdrs.pp.preprocess: consistency between sparse+cov+adj_prop and dense+cov+adj_prop
    """

    adata, adata_balance, df_cov = load_toy_data_adj_prop()

    adata_sparse = adata.copy()
    scdrs.pp.preprocess(adata_sparse, cov=df_cov, adj_prop="cell_group")

    adata_dense = adata.copy()
    adata_dense.X = adata_dense.X.toarray()
    scdrs.pp.preprocess(adata_dense, cov=df_cov, adj_prop="cell_group")

    # COV
    mat_X_sparse = (
        adata_sparse.X.toarray()
        + adata_sparse.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
            adata_sparse.uns["SCDRS_PARAM"]["COV_BETA"].values.T
        )
        + adata_sparse.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].values
    )
    max_abs_dif = np.absolute(mat_X_sparse - adata_dense.X).max()
    err_msg = (
        "Covariate-corrected matrices different between `sparse+cov+adj_prop` and `dense+cov+adj_prop`, max_abs_dif=%0.3e"
        % max_abs_dif
    )
    assert np.allclose(
        mat_X_sparse,
        adata_dense.X,
        rtol=1e-4,
        atol=1e-5,
        equal_nan=True,
    ), err_msg

    # GENE_STATS
    max_abs_dif = np.absolute(
        adata_sparse.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float)
        - adata_dense.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float)
    ).max()
    err_msg = (
        "'GENE_STATS' different between `sparse+cov+adj_prop` and `dense+cov+adj_prop`, max_abs_dif=%0.3e"
        % max_abs_dif
    )
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        adata_dense.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    # CELL_STATS
    max_abs_dif = np.absolute(
        adata_sparse.uns["SCDRS_PARAM"]["CELL_STATS"].values
        - adata_dense.uns["SCDRS_PARAM"]["CELL_STATS"].values
    ).max()
    err_msg = (
        "'CELL_STATS' different between `sparse+cov+adj_prop` and `dense+cov+adj_prop`, max_abs_dif=%0.3e"
        % max_abs_dif
    )
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        adata_dense.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    return


def test_preprocess_consistency_adj_prop_nocov():
    """
    Test scdrs.pp.preprocess: consistency between sparse+nocov+adj_prop and dense+nocov+adj_prop
    """

    adata, adata_balance, df_cov = load_toy_data_adj_prop()

    adata_sparse = adata.copy()
    scdrs.pp.preprocess(adata_sparse, adj_prop="cell_group")

    adata_dense = adata.copy()
    adata_dense.X = adata_dense.X.toarray()
    scdrs.pp.preprocess(adata_dense, adj_prop="cell_group")

    # GENE_STATS
    max_abs_dif = np.absolute(
        adata_sparse.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float)
        - adata_dense.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float)
    ).max()
    err_msg = (
        "'GENE_STATS' different between `sparse+nocov+adj_prop` and `dense+nocov+adj_prop`, max_abs_dif=%0.3e"
        % max_abs_dif
    )
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        adata_dense.uns["SCDRS_PARAM"]["GENE_STATS"].values.astype(float),
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    # CELL_STATS
    max_abs_dif = np.absolute(
        adata_sparse.uns["SCDRS_PARAM"]["CELL_STATS"].values
        - adata_dense.uns["SCDRS_PARAM"]["CELL_STATS"].values
    ).max()
    err_msg = (
        "'CELL_STATS' different between `sparse+nocov+adj_prop` and `dense+nocov+adj_prop`, max_abs_dif=%0.3e"
        % max_abs_dif
    )
    assert np.allclose(
        adata_sparse.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        adata_dense.uns["SCDRS_PARAM"]["CELL_STATS"].values,
        rtol=1e-4,
        equal_nan=True,
    ), err_msg

    return


def test_preprocess_accuracy_adj_prop():
    """
    Test scdrs.pp.preprocess: accuracy of `adj_prop` for approximating balanced data
    """

    adata, adata_balance, df_cov = load_toy_data_adj_prop()

    scdrs.pp.preprocess(adata_balance, cov=df_cov)

    adata_adj = adata.copy()
    scdrs.pp.preprocess(adata_adj, cov=df_cov, adj_prop="cell_group")

    adata_unadj = adata.copy()
    scdrs.pp.preprocess(adata_unadj, cov=df_cov)

    for term in ["mean", "var", "var_tech"]:
        v_balance = adata_balance.uns["SCDRS_PARAM"]["GENE_STATS"][term]
        v_adj = adata_adj.uns["SCDRS_PARAM"]["GENE_STATS"][term]
        v_unadj = adata_unadj.uns["SCDRS_PARAM"]["GENE_STATS"][term]
        mae_adj = np.absolute(v_balance - v_adj).mean()
        mae_unadj = np.absolute(v_balance - v_unadj).mean()
        err_msg = (
            "%s, MAE relative to balanced data: adjusted=%0.3f, unadjusted=%0.3f, ratio=%0.3f"
            % (term, mae_adj, mae_unadj, mae_adj / mae_unadj)
        )
        assert mae_adj / mae_unadj < 0.1, err_msg

    v_balance = adata_balance.uns["SCDRS_PARAM"]["GENE_STATS"]["mean_var"]
    v_adj = adata_adj.uns["SCDRS_PARAM"]["GENE_STATS"]["mean_var"]
    v_unadj = adata_unadj.uns["SCDRS_PARAM"]["GENE_STATS"]["mean_var"]
    acc_adj = (v_balance == v_adj).mean()
    acc_unadj = (v_balance == v_unadj).mean()
    err_msg = (
        "`mean_var` match with balanced data: adjusted=%0.3f, unadjusted=%0.3f, ratio=%0.3f"
        % (acc_adj, acc_unadj, acc_adj / acc_unadj)
    )
    assert acc_adj / acc_unadj > 1.5, err_msg

    return


def test_preprocess_accuracy_adj_prop_nocov():
    """
    Test scdrs.pp.preprocess: accuracy of `adj_prop` for approximating balanced data (nocov)
    """

    adata, adata_balance, df_cov = load_toy_data_adj_prop()

    scdrs.pp.preprocess(adata_balance)

    adata_adj = adata.copy()
    scdrs.pp.preprocess(adata_adj, adj_prop="cell_group")

    adata_unadj = adata.copy()
    scdrs.pp.preprocess(adata_unadj)

    for term in ["mean", "var", "var_tech"]:
        v_balance = adata_balance.uns["SCDRS_PARAM"]["GENE_STATS"][term]
        v_adj = adata_adj.uns["SCDRS_PARAM"]["GENE_STATS"][term]
        v_unadj = adata_unadj.uns["SCDRS_PARAM"]["GENE_STATS"][term]
        mae_adj = np.absolute(v_balance - v_adj).mean()
        mae_unadj = np.absolute(v_balance - v_unadj).mean()
        err_msg = (
            "%s, MAE relative to balanced data: adjusted=%0.3f, unadjusted=%0.3f, ratio=%0.3f"
            % (term, mae_adj, mae_unadj, mae_adj / mae_unadj)
        )
        assert mae_adj / mae_unadj < 0.1, err_msg

    v_balance = adata_balance.uns["SCDRS_PARAM"]["GENE_STATS"]["mean_var"]
    v_adj = adata_adj.uns["SCDRS_PARAM"]["GENE_STATS"]["mean_var"]
    v_unadj = adata_unadj.uns["SCDRS_PARAM"]["GENE_STATS"]["mean_var"]
    acc_adj = (v_balance == v_adj).mean()
    acc_unadj = (v_balance == v_unadj).mean()
    err_msg = (
        "`mean_var` match with balanced data: adjusted=%0.3f, unadjusted=%0.3f, ratio=%0.3f"
        % (acc_adj, acc_unadj, acc_adj / acc_unadj)
    )
    assert acc_adj / acc_unadj > 1.5, err_msg

    return


def test_compute_stats_weighted():
    """
    Test scdrs.pp.compute_stats: normal mode + cell_weight

    TODO:
    -----
    1. The following columns are not tested: df_gene["var_tech"],
    df_gene["ct_var_tech"], df_gene["mean_var"]
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    np.random.seed(0)
    cell_weight = np.random.rand(adata.shape[0]) * 10
    df_gene, df_cell = scdrs.pp.compute_stats(adata, cell_weight=cell_weight)

    mat_X = adata.X.toarray()
    mat_X_ct = np.expm1(mat_X)

    # df_gene["mean"]
    v_mean_true = np.average(mat_X, axis=0, weights=cell_weight)
    err_msg = (
        "avg_abs_gene_mean_dif=%0.2e"
        % np.absolute(df_gene["mean"] - v_mean_true).mean()
    )
    assert np.allclose(df_gene["mean"], v_mean_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["var"]
    v_var_true = np.average(mat_X ** 2, axis=0, weights=cell_weight)
    v_var_true = v_var_true - v_mean_true ** 2
    err_msg = (
        "avg_abs_gene_var_dif=%0.2e" % np.absolute(df_gene["var"] - v_var_true).mean()
    )
    assert np.allclose(df_gene["var"], v_var_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["ct_mean"]
    v_mean_true = np.average(mat_X_ct, axis=0, weights=cell_weight)
    err_msg = (
        "avg_abs_gene_ct_mean_dif=%0.2e"
        % np.absolute(df_gene["ct_mean"] - v_mean_true).mean()
    )
    assert np.allclose(
        df_gene["ct_mean"], v_mean_true, rtol=1e-4, equal_nan=True
    ), err_msg

    # df_gene["ct_var"]
    v_var_true = np.average(mat_X_ct ** 2, axis=0, weights=cell_weight)
    v_var_true = v_var_true - v_mean_true ** 2
    err_msg = (
        "avg_abs_gene_ct_var_dif=%0.2e"
        % np.absolute(df_gene["ct_var"] - v_var_true).mean()
    )
    assert np.allclose(
        df_gene["ct_var"], v_var_true, rtol=1e-4, equal_nan=True
    ), err_msg

    return


def test_compute_stats_implicit_covariate_correction_weighted():
    """
    Test scdrs.pp.compute_stats: implicit covariate correction mode + cell_weight

    TODO:
    -----
    1. The following columns are not tested: df_gene["var_tech"],
    df_gene["ct_var_tech"], df_gene["mean_var"]
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :50].copy()
    np.random.seed(0)
    cell_weight = np.random.rand(adata.shape[0]) * 10
    scdrs.pp.preprocess(adata, cov=df_cov)
    df_gene, df_cell = scdrs.pp.compute_stats(
        adata, implicit_cov_corr=True, cell_weight=cell_weight
    )

    # adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN
    mat_X = adata.X.toarray()
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
        adata.uns["SCDRS_PARAM"]["COV_BETA"].values.T
    )
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].values
    mat_X_ct = np.expm1(mat_X)

    # df_gene["mean"]
    v_mean_true = np.average(mat_X, axis=0, weights=cell_weight)
    err_msg = (
        "avg_abs_gene_mean_dif=%0.2e"
        % np.absolute(df_gene["mean"] - v_mean_true).mean()
    )
    assert np.allclose(df_gene["mean"], v_mean_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["var"]
    v_var_true = np.average(mat_X ** 2, axis=0, weights=cell_weight)
    v_var_true = v_var_true - v_mean_true ** 2
    err_msg = (
        "avg_abs_gene_var_dif=%0.2e" % np.absolute(df_gene["var"] - v_var_true).mean()
    )
    assert np.allclose(df_gene["var"], v_var_true, rtol=1e-4, equal_nan=True), err_msg

    # df_gene["ct_mean"]
    v_mean_true = np.average(mat_X_ct, axis=0, weights=cell_weight)
    err_msg = (
        "avg_abs_gene_ct_mean_dif=%0.2e"
        % np.absolute(df_gene["ct_mean"] - v_mean_true).mean()
    )
    assert np.allclose(
        df_gene["ct_mean"], v_mean_true, rtol=1e-4, equal_nan=True
    ), err_msg

    # df_gene["ct_var"]
    v_var_true = np.average(mat_X_ct ** 2, axis=0, weights=cell_weight)
    v_var_true = v_var_true - v_mean_true ** 2
    err_msg = (
        "avg_abs_gene_ct_var_dif=%0.2e"
        % np.absolute(df_gene["ct_var"] - v_var_true).mean()
    )
    assert np.allclose(
        df_gene["ct_var"], v_var_true, rtol=1e-4, equal_nan=True
    ), err_msg

    return


def test_get_mean_var_weighted():
    """
    Test scdrs.pp._get_mean_var with weights
    """

    mat_test = np.array([[1, 2], [3, 4], [5, 6]], dtype=float)

    for opt in ["dense", "sparse"]:
        for axis_ in [0, 1]:
            if opt == "sparse":
                mat_ = sparse.csr_matrix(mat_test)
            else:
                mat_ = mat_test

            if axis_ == 0:
                weights = np.array([8.8, 2.1, 5.6])
            else:
                weights = np.array([5.42, 4.66])

            v_mean, v_var = scdrs.pp._get_mean_var(mat_, axis=axis_, weights=weights)

            v_mean_true = np.average(mat_test, axis=axis_, weights=weights)
            v_var_true = np.average(mat_test ** 2, axis=axis_, weights=weights)
            v_var_true = v_var_true - v_mean_true ** 2

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


def test_get_mean_var_implicit_cov_corr_weighted():
    """
    Test scdrs.pp._get_mean_var_implicit_cov_corr with weights
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
        np.random.seed(0)
        weights = np.random.rand(mat_X.shape[axis_]) * 10

        v_mean, v_var = scdrs.pp._get_mean_var_implicit_cov_corr(
            adata, axis=axis_, weights=weights
        )
        v_mean_true = np.average(mat_X, axis=axis_, weights=weights)
        v_var_true = np.average(mat_X ** 2, axis=axis_, weights=weights)
        v_var_true = v_var_true - v_mean_true ** 2

        err_msg = "axis=%d, avg_abs_mean_dif=%0.2e, avg_abs_var_dif=%0.2e" % (
            axis_,
            np.absolute(v_mean - v_mean_true).mean(),
            np.absolute(v_var - v_var_true).mean(),
        )
        assert np.allclose(
            v_mean, v_mean_true, rtol=1e-4, equal_nan=True
        ) & np.allclose(v_var, v_var_true, rtol=1e-4, equal_nan=True), err_msg
    return
