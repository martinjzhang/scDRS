import pytest
import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
import scanpy as sc
from scipy.stats import rankdata

import scdrs
from .test_method_score_cell_main import load_toy_data


def gearyc_naive(mat_W, v_x):
    """
    Naive implementation of Geary's C, used for unit tests.
    """
    N = v_x.shape[0]
    C = 0
    for i in range(N):
        for j in range(N):
            C += mat_W[i, j] * (v_x[i] - v_x[j]) ** 2
    C = C * (N - 1) / 2 / np.sum(mat_W) / np.sum((v_x - np.mean(v_x)) ** 2)
    return C


def test_downstream_group_analysis():
    """
    Test scdrs.method.downstream_group_analysis
    @Kangcheng: could you review? Specifically, I could only reproduce an
    approximate result for 'hetero_mcz' (rtol=0.05). Could you check if I
    missed something?
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    df_full_score = dic_res_ref["REF_COV_FULL"]
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40, adata.n_obs - 1))

    dic_res = scdrs.method.downstream_group_analysis(
        adata, df_full_score, ["cell_type"]
    )

    assert "cell_type" in dic_res, "Expect 'cell_type' in `dic_res`."

    df_res = dic_res["cell_type"].loc[["causal_cell", "non_causal_cell"]]

    # Check n_cell and n_ctrl
    err_msg = "'n_cell': expect=[10, 20], actual=%s" % str(list(df_res["n_cell"]))
    assert list(df_res["n_cell"]) == [10, 20], err_msg
    err_msg = "'n_ctrl': expect=[20, 20], actual=%s" % str(list(df_res["n_ctrl"]))
    assert list(df_res["n_ctrl"]) == [20, 20], err_msg

    temp_df = df_full_score.join(adata.obs[["cell_type"]])
    ctrl_list = [x for x in temp_df if x.startswith("ctrl_norm_score_")]

    # Check assoc_mcp and assoc_mcz (5% quantile)
    v_mcp_true = np.zeros(2)
    v_mcz_true = np.zeros(2)
    for i_ct, ct in enumerate(["causal_cell", "non_causal_cell"]):
        cell_list_group = temp_df["cell_type"] == ct
        q95 = np.quantile(temp_df.loc[cell_list_group, "norm_score"], 0.95, axis=0)
        v_q95_ctrl = np.quantile(temp_df.loc[cell_list_group, ctrl_list], 0.95, axis=0)
        v_mcp_true[i_ct] = ((v_q95_ctrl >= q95).sum() + 1) / (v_q95_ctrl.shape[0] + 1)
        v_mcz_true[i_ct] = (q95 - v_q95_ctrl.mean()) / v_q95_ctrl.std()
    err_msg = "'assoc_mcp': expect=%s, actual=%s" % (
        str(list(v_mcp_true)),
        str(list(df_res["assoc_mcp"])),
    )
    assert np.allclose(df_res["assoc_mcp"], v_mcp_true), err_msg
    err_msg = "'assoc_mcz': expect=%s, actual=%s" % (
        str(list(v_mcz_true)),
        str(list(df_res["assoc_mcz"])),
    )
    assert np.allclose(df_res["assoc_mcz"], v_mcz_true), err_msg

    # Check hetero_mcp and hetero_mcz
    v_mcp_true = np.zeros(2)
    v_mcz_true = np.zeros(2)
    for i_ct, ct in enumerate(["causal_cell", "non_causal_cell"]):
        cell_list_group = temp_df["cell_type"] == ct
        adata_group = adata[cell_list_group, :].copy()
        mat_W = adata_group.obsp["connectivities"].toarray()
        # GC for disease score
        v_score_group = temp_df.loc[cell_list_group, "norm_score"].values
        gc = gearyc_naive(mat_W, v_score_group)
        # GC for disease-distribution-normalized control scores
        v_gc_ctrl = np.zeros(len(ctrl_list))
        for icol, col in enumerate(ctrl_list):
            v_ctrl_score = temp_df.loc[cell_list_group, col].values
            v_ctrl_score = np.sort(v_score_group)[
                rankdata(v_ctrl_score, method="ordinal") - 1
            ]
            v_gc_ctrl[icol] = gearyc_naive(mat_W, v_ctrl_score)
        v_mcp_true[i_ct] = ((v_gc_ctrl <= gc).sum() + 1) / (v_gc_ctrl.shape[0] + 1)
        # add ddof=1 to be consistent with the test_gearysc original implementation
        v_mcz_true[i_ct] = -(gc - v_gc_ctrl.mean()) / v_gc_ctrl.std(ddof=1)
    err_msg = "'hetero_mcp': expect=%s, actual=%s" % (
        str(list(v_mcp_true)),
        str(list(df_res["hetero_mcp"])),
    )
    assert np.allclose(df_res["hetero_mcp"], v_mcp_true), err_msg
    err_msg = "'hetero_mcz': expect=%s, actual=%s" % (
        str(list(v_mcz_true)),
        str(list(df_res["hetero_mcz"])),
    )
    assert np.allclose(df_res["hetero_mcz"], v_mcz_true), err_msg

    # Check n_fdr_0.05, n_fdr_0.1, n_fdr_0.2
    for col in ["n_fdr_0.05", "n_fdr_0.1", "n_fdr_0.2"]:
        err_msg = "'%s': expect=[5, 0], actual=%s" % (col, str(list(df_res[col])))
        assert list(df_res[col]) == [5, 0], err_msg

    return


def test_gearys_c():
    """
    Test scdrs.method.gearys_c
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=min(40, adata.n_obs - 1))

    v_x = np.arange(adata.shape[0])
    gc = scdrs.method.gearys_c(adata, v_x)
    gc_true = gearyc_naive(adata.obsp["connectivities"].toarray(), v_x)
    err_msg = "Geary's C: expect=%s, actual=%s" % (gc_true, gc)
    assert np.allclose(gc, gc_true), err_msg
    return


def test_downstream_corr_analysis():
    """
    Test scdrs.method.downstream_corr_analysis
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    df_full_score = dic_res_ref["REF_COV_FULL"]

    df_res = scdrs.method.downstream_corr_analysis(
        adata, df_full_score, ["causal_variable"]
    )
    n_ctrl, corr_mcp, corr_mcz = df_res.loc["causal_variable"]

    control_list = [x for x in df_full_score.columns if x.startswith("ctrl_norm_score")]
    v_x = adata.obs["causal_variable"].values
    corr_ = np.corrcoef(v_x, df_full_score["norm_score"])[0, 1]
    v_corr_ = np.array([np.corrcoef(v_x, df_full_score[x])[0, 1] for x in control_list])
    corr_mcp_true = ((v_corr_ >= corr_).sum() + 1) / (v_corr_.shape[0] + 1)
    corr_mcz_true = (corr_ - v_corr_.mean()) / v_corr_.std()

    err_msg = "'n_ctrl': expect=20, actual=%s" % n_ctrl
    assert n_ctrl == 20, err_msg

    err_msg = "'corr_mcp': expect=%s, actual=%s" % (corr_mcp_true, corr_mcp)
    assert np.allclose(corr_mcp, corr_mcp_true), err_msg

    err_msg = "'corr_mcz': expect=%s, actual=%s" % (corr_mcz_true, corr_mcz)
    assert np.allclose(corr_mcz, corr_mcz_true), err_msg
    return


def test_downstream_gene_analysis():
    """
    Test scdrs.method.downstream_gene_analysis
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :20].copy()
    adata_dense = adata.copy()
    adata_dense.X = adata_dense.X.toarray()
    df_full_score = dic_res_ref["REF_COV_FULL"]

    # Gene-analysis with sparse and dense data
    df_res = scdrs.method.downstream_gene_analysis(adata, df_full_score)
    df_res_dense = scdrs.method.downstream_gene_analysis(adata_dense, df_full_score)

    # True values
    df_res_true = pd.DataFrame(index=adata.var_names)
    df_res_true["CORR"] = [
        np.corrcoef(adata[:, x].X.toarray().flatten(), df_full_score["norm_score"])[
            0, 1
        ]
        for x in df_res_true.index
    ]
    df_res_true.sort_values("CORR", ascending=False, inplace=True)
    df_res_true["RANK"] = np.arange(df_res_true.shape[0])

    # Check sparse
    err_msg = (
        "'CORR' (sparse `adata.X`): ave_abs_dif=%0.2e"
        % np.absolute(df_res["CORR"] - df_res_true["CORR"]).mean()
    )
    assert np.allclose(df_res["CORR"], df_res_true["CORR"]), err_msg
    err_msg = (
        "'RANK' (sparse `adata.X`): ave_abs_dif=%0.2e"
        % np.absolute(df_res["RANK"] - df_res_true["RANK"]).mean()
    )
    assert np.allclose(df_res["RANK"], df_res_true["RANK"]), err_msg

    # Check dense
    err_msg = (
        "'CORR' (dense `adata.X`): ave_abs_dif=%0.2e"
        % np.absolute(df_res_dense["CORR"] - df_res_true["CORR"]).mean()
    )
    assert np.allclose(df_res_dense["CORR"], df_res_true["CORR"]), err_msg
    err_msg = (
        "'RANK' (dense `adata.X`): ave_abs_dif=%0.2e"
        % np.absolute(df_res_dense["RANK"] - df_res_true["RANK"]).mean()
    )
    assert np.allclose(df_res_dense["RANK"], df_res_true["RANK"]), err_msg

    return


def test_pearson_corr():
    """
    Test scdrs.method._pearson_corr
    """

    np.random.seed(0)
    mat_X = np.random.randn(3, 2)
    mat_X_sparse = sparse.csr_matrix(mat_X)
    mat_Y = np.random.randn(3, 2)
    mat_Y_sparse = sparse.csr_matrix(mat_Y)

    mat_corr_true = np.zeros([2, 2])
    for i in range(2):
        for j in range(2):
            mat_corr_true[i, j] = np.corrcoef(mat_X[:, i], mat_Y[:, j])[0, 1]

    # Basic
    mat_corr = scdrs.method._pearson_corr(mat_X, mat_Y)
    err_msg = (
        "`mat_X` and `mat_Y`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true), err_msg

    # Sparse vs. dense
    mat_corr = scdrs.method._pearson_corr(mat_X_sparse, mat_Y_sparse)
    err_msg = (
        "`mat_X_sparse` and `mat_Y_sparse`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true), err_msg

    mat_corr = scdrs.method._pearson_corr(mat_X_sparse, mat_Y)
    err_msg = (
        "`mat_X_sparse` and `mat_Y`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true), err_msg

    mat_corr = scdrs.method._pearson_corr(mat_X, mat_Y_sparse)
    err_msg = (
        "`mat_X` and `mat_Y_sparse`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true), err_msg

    # 1D vs. nD
    mat_corr = scdrs.method._pearson_corr(mat_X[:, 0], mat_Y[:, 0])
    err_msg = (
        "`mat_X[:,0]` and `mat_Y[:,0]`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true[0, 0]).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true[0, 0]), err_msg

    mat_corr = scdrs.method._pearson_corr(mat_X[:, 0], mat_Y)
    err_msg = (
        "`mat_X[:,0]` and `mat_Y`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true[0, :]).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true[0, :]), err_msg

    mat_corr = scdrs.method._pearson_corr(mat_X, mat_Y[:, 0])
    err_msg = (
        "`mat_X` and `mat_Y[:,0]`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true[:, 0]).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true[:, 0]), err_msg

    mat_corr = scdrs.method._pearson_corr(mat_X[:, 0], mat_Y_sparse)
    err_msg = (
        "`mat_X[:,0]` and `mat_Y_sparse`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true[0, :]).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true[0, :]), err_msg

    mat_corr = scdrs.method._pearson_corr(mat_X, mat_Y_sparse[:, 0])
    err_msg = (
        "`mat_X` and `mat_Y_sparse[:, 0]`: ave_abs_dif=%0.2e"
        % np.absolute(mat_corr - mat_corr_true[:, 0]).mean()
    )
    assert np.allclose(mat_corr, mat_corr_true[:, 0]), err_msg

    return
