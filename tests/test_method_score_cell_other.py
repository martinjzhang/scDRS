import pytest
import numpy as np
import scipy as sp
import pandas as pd

import scdrs
from .test_method_score_cell_main import load_toy_data


def test_select_ctrl_geneset():
    """
    Test scdrs.method._select_ctrl_geneset
    """

    np.random.seed(0)

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()

    df_gene = pd.DataFrame(
        index=adata.var_names[:100], columns=["gene", "categorical", "continuous"]
    )
    df_gene["gene"] = df_gene.index
    df_gene["categorical"] = [1] * 50 + [2] * 50
    df_gene["continuous"] = np.random.rand(100)

    gene_list = list(df_gene.index[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.5, 3.8, 4.1, 5.2]

    for ctrl_match_key in ["categorical", "continuous"]:
        dic_ctrl_list, dic_ctrl_weight = scdrs.method._select_ctrl_geneset(
            df_gene, gene_list, gene_weight, ctrl_match_key, 1, 10, 0
        )

        ctrl_gene_list_sort = np.array(dic_ctrl_list[0])[np.argsort(dic_ctrl_weight[0])]
        ctrl_weight_sort = np.array(dic_ctrl_weight[0])[np.argsort(dic_ctrl_weight[0])]

        err_msg = "ctrl_match_key={}\n".format(ctrl_match_key)
        err_msg += "|{:^15}|{:^15}|{:^15}|{:^15}|{:^15}|{:^15}|\n".format(
            "GENE", ctrl_match_key, "WEIGHT", "CTRL_GENE", ctrl_match_key, "WEIGHT"
        )
        for i in range(len(gene_list)):
            err_msg += (
                "|{:^15}|{:^15.3f}|{:^15.3f}|{:^15}|{:^15.3f}|{:^15.3f}|\n".format(
                    gene_list[i],
                    df_gene.loc[gene_list[i], ctrl_match_key],
                    gene_weight[i],
                    ctrl_gene_list_sort[i],
                    df_gene.loc[ctrl_gene_list_sort[i], ctrl_match_key],
                    ctrl_weight_sort[i],
                )
            )
        assert (
            np.allclose(
                df_gene.loc[gene_list, ctrl_match_key],
                df_gene.loc[ctrl_gene_list_sort, ctrl_match_key],
                rtol=0,
                atol=0.1,
            )
            & np.allclose(gene_weight, ctrl_weight_sort)
        ), err_msg
    return


def test_compute_raw_score_dense_nocov_vs():
    """
    Test scdrs.method._compute_raw_score: dense+nocov+vs
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])

    adata.X = adata.X.toarray()
    scdrs.pp.preprocess(adata, cov=None)
    v_raw_score, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, np.ones(len(gene_list)), "vs"
    )

    v_score_weight_true = 1 / np.sqrt(
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values + 1e-2
    )
    v_score_weight_true = v_score_weight_true / v_score_weight_true.sum()
    v_raw_score_true = adata[:, gene_list].X.dot(v_score_weight_true)

    err_msg = "Dense+nocov+vs: avg_abs_score_dif=%0.2e, avg_abs_weight_dif=%0.2e" % (
        np.absolute(v_raw_score - v_raw_score_true).mean(),
        np.absolute(v_score_weight - v_score_weight_true).mean(),
    )
    assert np.allclose(v_raw_score, v_raw_score_true) & np.allclose(
        v_score_weight, v_score_weight_true
    ), err_msg

    return


def test_compute_raw_score_dense_nocov_vs_weight():
    """
    Test scdrs.method._compute_raw_score: dense+nocov+vs+weight
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.3, 1.8, 0.2, 3]

    adata.X = adata.X.toarray()
    scdrs.pp.preprocess(adata, cov=None)
    v_raw_score, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, gene_weight, "vs"
    )

    v_score_weight_true = 1 / np.sqrt(
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values + 1e-2
    )
    v_score_weight_true = v_score_weight_true * np.array(gene_weight)
    v_score_weight_true = v_score_weight_true / v_score_weight_true.sum()
    v_raw_score_true = adata[:, gene_list].X.dot(v_score_weight_true)

    err_msg = (
        "Dense+nocov+vs+weight: avg_abs_score_dif=%0.2e, avg_abs_weight_dif=%0.2e"
        % (
            np.absolute(v_raw_score - v_raw_score_true).mean(),
            np.absolute(v_score_weight - v_score_weight_true).mean(),
        )
    )
    assert np.allclose(v_raw_score, v_raw_score_true) & np.allclose(
        v_score_weight, v_score_weight_true
    ), err_msg

    return


def test_compute_raw_score_sparse_nocov_vs():
    """
    Test scdrs.method._compute_raw_score: sparse+nocov+vs
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])

    scdrs.pp.preprocess(adata, cov=None)
    v_raw_score, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, np.ones(len(gene_list)), "vs"
    )

    v_score_weight_true = 1 / np.sqrt(
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values + 1e-2
    )
    v_score_weight_true = v_score_weight_true / v_score_weight_true.sum()
    v_raw_score_true = adata[:, gene_list].X.dot(v_score_weight_true)

    err_msg = "Sparse+nocov+vs: avg_abs_score_dif=%0.2e, avg_abs_weight_dif=%0.2e" % (
        np.absolute(v_raw_score - v_raw_score_true).mean(),
        np.absolute(v_score_weight - v_score_weight_true).mean(),
    )
    assert np.allclose(v_raw_score, v_raw_score_true) & np.allclose(
        v_score_weight, v_score_weight_true
    ), err_msg

    return


def test_compute_raw_score_sparse_cov_vs():
    """
    Test scdrs.method._compute_raw_score: sparse+nocov+vs
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])

    scdrs.pp.preprocess(adata, cov=df_cov)
    v_raw_score, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, np.ones(len(gene_list)), "vs"
    )

    # adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN
    mat_X = adata[:, gene_list].X.toarray()
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
        adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list].values.T
    )
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

    v_score_weight_true = 1 / np.sqrt(
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values + 1e-2
    )
    v_score_weight_true = v_score_weight_true / v_score_weight_true.sum()
    v_raw_score_true = mat_X.dot(v_score_weight_true)

    err_msg = "Sparse+cov+vs: avg_abs_score_dif=%0.2e, avg_abs_weight_dif=%0.2e" % (
        np.absolute(v_raw_score - v_raw_score_true).mean(),
        np.absolute(v_score_weight - v_score_weight_true).mean(),
    )
    assert np.allclose(v_raw_score, v_raw_score_true) & np.allclose(
        v_score_weight, v_score_weight_true
    ), err_msg

    return


def test_compute_raw_score_sparse_cov_vs_weight():
    """
    Test scdrs.method._compute_raw_score: sparse+nocov+vs+weight
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.3, 1.8, 0.2, 3]

    scdrs.pp.preprocess(adata, cov=df_cov)
    v_raw_score, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, gene_weight, "vs"
    )

    # adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN
    mat_X = adata[:, gene_list].X.toarray()
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
        adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list].values.T
    )
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

    v_score_weight_true = 1 / np.sqrt(
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values + 1e-2
    )
    v_score_weight_true = v_score_weight_true * np.array(gene_weight)
    v_score_weight_true = v_score_weight_true / v_score_weight_true.sum()
    v_raw_score_true = mat_X.dot(v_score_weight_true)

    err_msg = (
        "Sparse+cov+vs+weight: avg_abs_score_dif=%0.2e, avg_abs_weight_dif=%0.2e"
        % (
            np.absolute(v_raw_score - v_raw_score_true).mean(),
            np.absolute(v_score_weight - v_score_weight_true).mean(),
        )
    )
    assert np.allclose(v_raw_score, v_raw_score_true) & np.allclose(
        v_score_weight, v_score_weight_true
    ), err_msg

    return


def test_compute_raw_score_sparse_cov_uniform():
    """
    Test scdrs.method._compute_raw_score: sparse+nocov+uniform
    """

    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])

    scdrs.pp.preprocess(adata, cov=df_cov)
    v_raw_score, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, np.ones(len(gene_list)), "uniform"
    )

    # adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN
    mat_X = adata[:, gene_list].X.toarray()
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_MAT"].values.dot(
        adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list].values.T
    )
    mat_X = mat_X + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

    v_score_weight_true = np.ones(len(gene_list)) * 0.2
    v_score_weight_true = v_score_weight_true / v_score_weight_true.sum()
    v_raw_score_true = mat_X.dot(v_score_weight_true)

    err_msg = (
        "Sparse+cov+uniform: avg_abs_score_dif=%0.2e, avg_abs_weight_dif=%0.2e"
        % (
            np.absolute(v_raw_score - v_raw_score_true).mean(),
            np.absolute(v_score_weight - v_score_weight_true).mean(),
        )
    )
    assert np.allclose(v_raw_score, v_raw_score_true) & np.allclose(
        v_score_weight, v_score_weight_true
    ), err_msg

    return


def test_compute_overdispersion_score_cov():
    """
    Test scdrs.method._compute_overdispersion_score: sparse+cov and dense+cov
    """

    # Sparse+cov
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.3, 1.8, 0.2, 3]

    scdrs.pp.preprocess(adata, cov=df_cov)
    v_raw_score_sparse, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, gene_weight, "od"
    )

    # Dense+cov
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.3, 1.8, 0.2, 3]

    adata.X = adata.X.toarray()
    scdrs.pp.preprocess(adata, cov=df_cov)
    v_raw_score_dense, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, gene_weight, "od"
    )

    # True val (from dense+cov)
    mat_X = adata[:, gene_list].X.toarray()

    v_mean = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "mean"].values
    v_var_tech = (
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values
    )
    v_w = np.array(gene_weight) / (v_var_tech + 1e-2)
    v_w = v_w / v_w.sum()

    v_raw_score_true = (((mat_X - v_mean) ** 2 - v_var_tech) * v_w).sum(axis=1)

    err_msg = "Sparse+cov+od: avg_abs_score_dif=%0.2e" % (
        np.absolute(v_raw_score_sparse - v_raw_score_true).mean(),
    )
    assert np.allclose(v_raw_score_sparse, v_raw_score_true, rtol=1e-4), err_msg

    err_msg = "Dense+cov+od: avg_abs_score_dif=%0.2e" % (
        np.absolute(v_raw_score_dense - v_raw_score_true).mean(),
    )
    assert np.allclose(v_raw_score_dense, v_raw_score_true, rtol=1e-4), err_msg

    return


def test_compute_overdispersion_score_nocov():
    """
    Test scdrs.method._compute_overdispersion_score: sparse+nocov and dense+nocov
    """

    # Sparse+nocov
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.3, 1.8, 0.2, 3]

    scdrs.pp.preprocess(adata, cov=None)
    v_raw_score_sparse, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, gene_weight, "od"
    )

    # Dense+nocov
    adata, df_cov, df_gs, dic_res_ref = load_toy_data()
    adata = adata[:, :100].copy()
    gene_list = list(adata.var_names[[0, 55, 27, 80, 2]])
    gene_weight = [1.1, 2.3, 1.8, 0.2, 3]

    adata.X = adata.X.toarray()
    scdrs.pp.preprocess(adata, cov=None)
    v_raw_score_dense, v_score_weight = scdrs.method._compute_raw_score(
        adata, gene_list, gene_weight, "od"
    )

    # True val (from dense+cov)
    mat_X = adata[:, gene_list].X.toarray()

    v_mean = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "mean"].values
    v_var_tech = (
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values
    )
    v_w = np.array(gene_weight) / (v_var_tech + 1e-2)
    v_w = v_w / v_w.sum()

    v_raw_score_true = (((mat_X - v_mean) ** 2 - v_var_tech) * v_w).sum(axis=1)

    err_msg = "Sparse+nocov+od: avg_abs_score_dif=%0.2e" % (
        np.absolute(v_raw_score_sparse - v_raw_score_true).mean(),
    )
    assert np.allclose(v_raw_score_sparse, v_raw_score_true, rtol=1e-4), err_msg

    err_msg = "Dense+nocov+od: avg_abs_score_dif=%0.2e" % (
        np.absolute(v_raw_score_dense - v_raw_score_true).mean(),
    )
    assert np.allclose(v_raw_score_dense, v_raw_score_true, rtol=1e-4), err_msg

    return


def test_correct_background_causal_cell():
    """
    Test scdrs.method._correct_background: check if it can find the causal cell
    """

    dic_norm_score = {}
    dic_ctrl_score = {}

    for config in [
        "no-perturb",
        "zero-mask",
        "cell-perturb",
        "gs-perturb",
        "all-perturb",
    ]:

        np.random.seed(0)
        v_raw_score = np.random.rand(500)
        v_raw_score[5] = 3  # Causal cell
        mat_ctrl_raw_score = np.random.rand(500, 100)
        ind_zero_ctrl_score = mat_ctrl_raw_score < 0.01

        # Cell-wise transformation
        if config in ["cell-perturb", "all-perturb"]:
            v_translate_cell = np.random.rand(500) * 1
            v_raw_score = v_raw_score + v_translate_cell
            mat_ctrl_raw_score = (mat_ctrl_raw_score.T + v_translate_cell).T

            v_scale_cell = np.random.rand(500) * 0.5 + 0.5
            v_raw_score = v_raw_score * v_scale_cell
            mat_ctrl_raw_score = (mat_ctrl_raw_score.T * v_scale_cell).T

        # Gene-gene-wise scaling
        v_var_ratio_c2t = np.ones(100)
        if config in ["gs-perturb", "all-perturb"]:
            mat_ctrl_raw_score += np.random.rand(100) * 1
            v_var_ratio_c2t = np.random.rand(100) * 0.5 + 0.5
            mat_ctrl_raw_score = mat_ctrl_raw_score * np.sqrt(v_var_ratio_c2t)

        if config in ["zero-mask", "all-perturb"]:
            mat_ctrl_raw_score[ind_zero_ctrl_score] = 0

        v_norm_score, mat_ctrl_norm_score = scdrs.method._correct_background(
            v_raw_score, mat_ctrl_raw_score, v_var_ratio_c2t, save_intermediate=None
        )

        dic_norm_score[config] = v_norm_score
        dic_ctrl_score[config] = mat_ctrl_norm_score

        # Check causal cell
        err_msg = "Config=%s, cell #5 should be causal (norm disease score=%0.2e)" % (
            config,
            v_norm_score[5],
        )
        assert v_norm_score[5] > 3, err_msg

    return


def test_get_p_from_empi_null():
    """
    Test scdrs.method._get_p_from_empi_null
    """
    v_t = [0, 1]
    v_t_null = [0.5, 0.6, 0.7]
    v_p = scdrs.method._get_p_from_empi_null(v_t, v_t_null)
    assert np.absolute(v_p[0] - 1) < 0.001, "First MC-p should be 1"
    assert np.absolute(v_p[1] - 0.25) < 0.001, "Second MC-p should be 1/(3+1)=0.25"
    return