import pytest
import numpy as np
import scipy as sp
import pandas as pd
import os
from anndata import read_h5ad
import scdrs


def load_toy_data():
    """
    Load toy data
    """

    DATA_PATH = scdrs.__path__[0]

    H5AD_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.h5ad")
    COV_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.cov")
    GS_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.gs")
    REF_NOCOV_FILE = os.path.join(
        DATA_PATH, "data/toydata_gs_mouse.ref_Ctrl20_CovNone.score.gz"
    )
    REF_COV_FILE = os.path.join(
        DATA_PATH, "data/toydata_gs_mouse.ref_Ctrl20_CovConstCovariate.score.gz"
    )

    assert os.path.exists(H5AD_FILE), "built-in data toydata_mouse.h5ad missing"
    assert os.path.exists(COV_FILE), "built-in data toydata_mouse.cov missing"
    assert os.path.exists(GS_FILE), "built-in data toydata_mouse.gs missing"
    assert os.path.exists(
        REF_NOCOV_FILE
    ), "built-in data toydata_gs_mouse.ref_Ctrl20_CovNone.score.gz missing"
    assert os.path.exists(
        REF_COV_FILE
    ), "built-in data toydata_gs_mouse.ref_Ctrl20_CovConstCovariate.score.gz missing"

    adata = read_h5ad(H5AD_FILE)
    df_cov = pd.read_csv(COV_FILE, sep="\t", index_col=0)
    df_gs = pd.read_csv(GS_FILE, sep="\t")
    df_gs.index = df_gs["TRAIT"]
    dic_res_ref = {
        "REF_NOCOV": pd.read_csv(REF_NOCOV_FILE, sep="\t", index_col=0),
        "REF_COV": pd.read_csv(REF_COV_FILE, sep="\t", index_col=0),
    }

    return adata, df_cov, df_gs, dic_res_ref


def compare_score_file(df_res, df_res_ref):
    """
    Compare df_res
    """

    col_list = ["raw_score", "norm_score", "mc_pval", "pval"]
    for col in col_list:
        v_ = df_res[col].values
        v_ref = df_res_ref[col].values
        err_msg = "Inconsistent values: {}\n".format(col)
        err_msg += "|{:^15}|{:^15}|{:^15}|{:^15}|\n".format(
            "OBS", "REF", "DIF", "REL_DIF"
        )
        for i in range(v_.shape[0]):
            err_msg += "|{:^15.3e}|{:^15.3e}|{:^15.3e}|{:^15.3e}|\n".format(
                v_[i],
                v_ref[i],
                v_[i] - v_ref[i],
                np.absolute((v_[i] - v_ref[i]) / v_ref[i]),
            )
        assert np.allclose(v_, v_ref, rtol=1e-3, equal_nan=True), err_msg
    return None


def test_score_cell_nocov():
    """
    score_cell: dense + nocov
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    gene_list = df_gs.loc["toydata_gs_mouse", "GENESET"].split(",")
    gene_list = sorted(set(gene_list) & set(adata.var_names))

    scdrs.method.compute_stats(adata)
    df_res = scdrs.method.score_cell(
        adata,
        gene_list,
        ctrl_match_key="mean_var",
        n_ctrl=20,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=False,
        verbose=False,
    )
    compare_score_file(df_res, dic_res_ref["REF_NOCOV"])

    return None


def test_score_cell_dense_cov():
    """
    score_cell: dense + cov
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    gene_list = df_gs.loc["toydata_gs_mouse", "GENESET"].split(",")
    gene_list = sorted(set(gene_list) & set(adata.var_names))

    # Load cov
    cov_list = list(df_cov.columns)
    if len(set(df_cov.index) & set(adata.obs_names)) < 0.1 * adata.shape[0]:
        raise ValueError("--cov_file does not match the cells in --h5ad_file")
    adata.obs.drop(
        [x for x in cov_list if x in adata.obs.columns], axis=1, inplace=True
    )
    adata.obs = adata.obs.join(df_cov)
    adata.obs.fillna(adata.obs[cov_list].mean(), inplace=True)
    # fixit
    v_resid = scdrs.method.reg_out(np.ones(adata.shape[0]), adata.obs[cov_list])
    if ((v_resid ** 2).mean() > 0.05) and ("const" not in cov_list):
        adata.obs["const"] = 1
        cov_list = ["const"] + cov_list
    # fixit

    adata.var["mean"] = adata.X.mean(axis=0).T
    if sp.sparse.issparse(adata.X):
        adata.X = adata.X.toarray()
    adata.X -= adata.var["mean"].values
    adata.X = scdrs.method.reg_out(adata.X, adata.obs[cov_list].values)
    adata.X += adata.var["mean"]

    # Compute score
    scdrs.method.compute_stats(adata)
    df_res = scdrs.method.score_cell(
        adata,
        gene_list,
        ctrl_match_key="mean_var",
        n_ctrl=20,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=False,
        verbose=False,
    )
    compare_score_file(df_res, dic_res_ref["REF_COV"])

    return None


def test_score_cell_sparse_cov():
    """
    score_cell: sparse + cov
    """
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    gene_list = df_gs.loc["toydata_gs_mouse", "GENESET"].split(",")
    gene_list = sorted(set(gene_list) & set(adata.var_names))

    # Load cov
    cov_list = list(df_cov.columns)
    if len(set(df_cov.index) & set(adata.obs_names)) < 0.1 * adata.shape[0]:
        raise ValueError("--cov_file does not match the cells in --h5ad_file")
    adata.obs.drop(
        [x for x in cov_list if x in adata.obs.columns], axis=1, inplace=True
    )
    adata.obs = adata.obs.join(df_cov)
    adata.obs.fillna(adata.obs[cov_list].mean(), inplace=True)
    # fixit
    v_resid = scdrs.method.reg_out(np.ones(adata.shape[0]), adata.obs[cov_list])
    if ((v_resid ** 2).mean() > 0.05) and ("const" not in cov_list):
        adata.obs["const"] = 1
        cov_list = ["const"] + cov_list
    # fixit

    cov_values = adata.obs[cov_list].values
    scdrs.sparse_util.sparse_reg_out(adata, cov_values)

    # Compute score
    scdrs.sparse_util.sparse_compute_stats(adata)
    df_res = scdrs.method.score_cell(
        adata,
        gene_list,
        ctrl_match_key="mean_var",
        n_ctrl=20,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=False,
        verbose=False,
        sparse=True,
    )
    compare_score_file(df_res, dic_res_ref["REF_COV"])

    return None
