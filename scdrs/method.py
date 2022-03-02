import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
from tqdm import tqdm
import anndata
from typing import List, Dict, Tuple
from statsmodels.stats.multitest import multipletests
import scdrs


def score_cell(
    data,
    gene_list,
    gene_weight=None,
    ctrl_match_key="mean_var",
    n_ctrl=1000,
    n_genebin=200,
    weight_opt="vs",
    copy=False,
    return_ctrl_raw_score=False,
    return_ctrl_norm_score=False,
    random_seed=0,
    verbose=False,
    save_intermediate=None,
):

    """Score cells based on the disease gene set.

    Preprocessing information `data.uns["SCDRS_PARAM"]` is required
    (run `scdrs.pp.preprocess` first).

    It operates in implicit-covariate-correction mode if both `FLAG_SPARSE`
    and `FLAG_COV` are `True`, where computations are based on the implicit
    covariate-corrected data 
    
        `CORRECTED_X = data.X + COV_MAT * COV_BETA + COV_GENE_MEAN`.

    It operates in normal mode otherwise, where computations are based on `data.X`,


    Parameters
    ----------
    data : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : array_like, default=None
        Gene weights of length n_disease_gene for genes in the gene_list.
        If gene_weight=None, the weights are set to be one.
    ctrl_match_key : str, default="mean_var"
        Gene-level statistic used for matching control and disease genes;
        should be in `data.uns["SCDRS_PARAM"]["GENE_STATS"]`.
    n_ctrl : int, default=1000
        Number of control gene sets.
    n_genebin : int, default=200
        Number of bins for dividing genes by ctrl_match_key if
        `data.uns["SCDRS_PARAM"]["GENE_STATS"][ctrl_match_key]` is a continuous variable.
    weight_opt : str, default="vs"
        Option for computing the raw score
        
        - 'uniform': average over the genes in the gene_list.
        - 'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct).
        - 'inv_std': weighted average with weights equal to 1/std.
        - 'od': overdispersion score.
        
    copy : bool, default=False
        If to make copy of the AnnData object to avoid writing on the orignal data.
    return_raw_ctrl_score : bool, default=False
        If to return raw control scores.
    return_norm_ctrl_score : bool, default=False
        If to return normalized control scores.
    random_seed : int, default=0
        Random seed.
    verbose : bool, default=False
        If to output messages.
    save_intermediate : str, default=None
        File path prefix for saving intermediate results.

    Returns
    -------    
    df_res : pandas.DataFrame (dtype=np.float32)
        scDRS results of shape (n_cell, n_key) with columns
        
        - raw_score: raw disease scores.
        - norm_score: normalized disease scores.
        - mc_pval: Monte Carlo p-values based on the normalized control scores of the same cell.
        - pval: scDRS individual cell-level disease-association p-values.
        - nlog10_pval: -log10(pval). Needed in case the single precision (np.float32) gives inaccurate p-values
        - zscore: one-side z-score converted from pval.
        - ctrl_raw_score_*: raw control scores.
        - ctrl_norm_score_*: normalized control scores.
    """

    np.random.seed(random_seed)
    adata = data.copy() if copy else data
    n_cell, n_gene = adata.shape

    # Check preprocessing information
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `scdrs.pp.preprocess` first"

    # Check GENE_STATS from adata.uns["SCDRS_PARAM"]
    assert (
        "GENE_STATS" in adata.uns["SCDRS_PARAM"]
    ), "adata.uns['SCDRS_PARAM']['GENE_STATS'] not found, run `scdrs.pp.preprocess` first"

    gene_stats_set_expect = {"mean", "var", "var_tech"}
    gene_stats_set = set(adata.uns["SCDRS_PARAM"]["GENE_STATS"])
    assert (
        len(gene_stats_set_expect - gene_stats_set) == 0
    ), "One of 'mean', 'var', 'var_tech' not found in adata.uns['SCDRS_PARAM']['GENE_STATS'], run `scdrs.pp.preprocess` first"

    # Check if ctrl_match_key is in GENE_STATS
    assert ctrl_match_key in adata.uns["SCDRS_PARAM"]["GENE_STATS"], (
        "ctrl_match_key=%s not found in adata.uns['SCDRS_PARAM']['GENE_STATS']"
        % ctrl_match_key
    )

    # Check if weight_opt is legal
    assert weight_opt in {"uniform", "vs", "inv_std", "od"}, (
        "weight_opt=%s is not one of {'uniform', 'vs', 'inv_std', 'od}'" % weight_opt
    )

    if verbose:
        msg = "# scdrs.method.score_cell summary:"
        msg += "\n    n_cell=%d, n_gene=%d," % (n_cell, n_gene)
        msg += "\n    n_disease_gene=%d," % len(gene_list)
        msg += "\n    n_ctrl=%d, n_genebin=%d," % (n_ctrl, n_genebin)
        msg += "\n    ctrl_match_key='%s'," % ctrl_match_key
        msg += "\n    weight_opt='%s'," % weight_opt
        msg += "\n    return_ctrl_raw_score=%s," % return_ctrl_raw_score
        msg += "\n    return_ctrl_norm_score=%s," % return_ctrl_norm_score
        msg += "\n    random_seed=%d, verbose=%s," % (random_seed, verbose)
        msg += "\n    save_intermediate=%s," % save_intermediate
        print(msg)

    # Load parameters
    flag_sparse = adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"]
    flag_cov = adata.uns["SCDRS_PARAM"]["FLAG_COV"]

    df_gene = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[adata.var_names].copy()
    df_gene["gene"] = df_gene.index
    df_gene.drop_duplicates(subset="gene", inplace=True)

    gene_list = list(gene_list)
    if gene_weight is not None:
        gene_weight = list(gene_weight)
    else:
        gene_weight = [1] * len(gene_list)

    # Overlap gene_list with df_gene["gene"]
    dic_gene_weight = {x: y for x, y in zip(gene_list, gene_weight)}
    gene_list = sorted(set(gene_list) & set(df_gene["gene"]))
    gene_weight = [dic_gene_weight[x] for x in gene_list]

    if verbose:
        print(
            "# scdrs.method.score_cell: use %d overlapping genes for scoring"
            % len(gene_list)
        )

    # Select control gene sets
    dic_ctrl_list, dic_ctrl_weight = _select_ctrl_geneset(
        df_gene, gene_list, gene_weight, ctrl_match_key, n_ctrl, n_genebin, random_seed
    )

    # Compute raw scores
    v_raw_score, v_score_weight = _compute_raw_score(
        adata, gene_list, gene_weight, weight_opt
    )

    mat_ctrl_raw_score = np.zeros([n_cell, n_ctrl])
    mat_ctrl_weight = np.zeros([len(gene_list), n_ctrl])
    for i_ctrl in tqdm(range(n_ctrl), desc="Computing control scores"):
        v_ctrl_raw_score, v_ctrl_weight = _compute_raw_score(
            adata, dic_ctrl_list[i_ctrl], dic_ctrl_weight[i_ctrl], weight_opt
        )
        mat_ctrl_raw_score[:, i_ctrl] = v_ctrl_raw_score
        mat_ctrl_weight[:, i_ctrl] = v_ctrl_weight

    # Compute normalized scores
    v_var_ratio_c2t = np.ones(n_ctrl)
    if (ctrl_match_key == "mean_var") & (weight_opt in ["uniform", "vs", "inv_std"]):
        # For mean_var matched control genes and raw scores computed as weighted average,
        # estimate variance ratio assuming independence.
        for i_ctrl in range(n_ctrl):
            v_var_ratio_c2t[i_ctrl] = (
                df_gene.loc[dic_ctrl_list[i_ctrl], "var"]
                * mat_ctrl_weight[:, i_ctrl] ** 2
            ).sum()
        v_var_ratio_c2t /= (df_gene.loc[gene_list, "var"] * v_score_weight ** 2).sum()

    v_norm_score, mat_ctrl_norm_score = _correct_background(
        v_raw_score,
        mat_ctrl_raw_score,
        v_var_ratio_c2t,
        save_intermediate=save_intermediate,
    )

    # Get p-values
    mc_p = (1 + (mat_ctrl_norm_score.T >= v_norm_score).sum(axis=0)) / (1 + n_ctrl)
    pooled_p = _get_p_from_empi_null(v_norm_score, mat_ctrl_norm_score.flatten())
    nlog10_pooled_p = -np.log10(pooled_p)
    pooled_z = -sp.stats.norm.ppf(pooled_p).clip(min=-10, max=10)

    # Return result
    dic_res = {
        "raw_score": v_raw_score,
        "norm_score": v_norm_score,
        "mc_pval": mc_p,
        "pval": pooled_p,
        "nlog10_pval": nlog10_pooled_p,
        "zscore": pooled_z,
    }
    if return_ctrl_raw_score:
        for i in range(n_ctrl):
            dic_res["ctrl_raw_score_%d" % i] = mat_ctrl_raw_score[:, i]
    if return_ctrl_norm_score:
        for i in range(n_ctrl):
            dic_res["ctrl_norm_score_%d" % i] = mat_ctrl_norm_score[:, i]
    df_res = pd.DataFrame(index=adata.obs.index, data=dic_res, dtype=np.float32)
    return df_res


def _select_ctrl_geneset(
    input_df_gene,
    gene_list,
    gene_weight,
    ctrl_match_key,
    n_ctrl,
    n_genebin,
    random_seed,
):

    """Subroutine for `scdrs.method.score_cell`. Select control gene sets that match
    the disease gene set by `ctrl_match_key`.

    It recognizes `ctrl_match_key` as categorical if the number of unique values is
    less than 10% of the total number of values, and otherwise continuous. For
    categorical `ctrl_match_key`, genes are matched within each category. For continuous
    `ctrl_match_key`, genes are divided into `n_genebin` bins and are matched within
    each bin. A matched control gene takes the same weight as the disease gene,

    Args
    ----
    input_df_gene : pd.DataFrame
        Gene-level statistics of shape (n_gene, n_stats).
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : list
        Gene weights of length n_disease_gene for genes in the gene_list.
    ctrl_match_key : str
        Gene-level statistic used for matching control and disease genes;
        should be in `input_df_gene`.
    n_ctrl : int
        Number of control gene sets.
    n_genebin : int
        Number of bins for dividing genes by ctrl_match_key if
        `input_df_gene[ctrl_match_key]` is a continuous variable.
    random_seed : int
        Random seed.

    Returns
    -------
    dic_ctrl_list : dict of lists
        dic_ctrl_list[i]: the i-th control gene list
    dic_ctrl_weight : dict of lists
        dic_ctrl_weight[i]: weights for the i-th control gene list

    """

    np.random.seed(random_seed)
    df_gene = input_df_gene.copy()
    if "gene" not in df_gene:
        df_gene["gene"] = df_gene.index
    disease_gene_set = set(gene_list)
    dic_gene_weight = {x: y for x, y in zip(gene_list, gene_weight)}

    # Divide genes into equal-sized bins based on ctrl_match_key
    if df_gene[ctrl_match_key].unique().shape[0] < df_gene.shape[0] / 10:
        df_gene_bin = df_gene.groupby(ctrl_match_key).agg({"gene": set})
    else:
        df_gene["qbin"] = pd.qcut(
            df_gene[ctrl_match_key], q=n_genebin, labels=False, duplicates="drop"
        )
        df_gene_bin = df_gene.groupby("qbin").agg({"gene": set})

    # Find ctrl_match_key matched control genes
    dic_ctrl_list = {x: [] for x in range(n_ctrl)}
    dic_ctrl_weight = {x: [] for x in range(n_ctrl)}
    for bin_ in df_gene_bin.index:
        bin_gene = sorted(df_gene_bin.loc[bin_, "gene"])
        bin_disease_gene = sorted(df_gene_bin.loc[bin_, "gene"] & disease_gene_set)
        if len(bin_disease_gene) > 0:
            for i_list in np.arange(n_ctrl):
                dic_ctrl_list[i_list].extend(
                    np.random.choice(
                        bin_gene, size=len(bin_disease_gene), replace=False
                    )
                )
                dic_ctrl_weight[i_list].extend(
                    [dic_gene_weight[x] for x in bin_disease_gene]
                )

    return dic_ctrl_list, dic_ctrl_weight


def _compute_raw_score(adata, gene_list, gene_weight, weight_opt):
    """Compute raw score
        v_score_weight = gene_weight * {uniform/vs/inv_std}
        `SCDRS_PARAM` is assumed to have been computed using `sparse_reg_out`

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : list
        Gene weights of length n_disease_gene for genes in the gene_list.
    weight_opt : str
        Option for computing the raw score
        - 'uniform': average over the genes in the gene_list.
        - 'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct).
        - 'inv_std': weighted average with weights equal to 1/std.
        - 'od': overdispersion score.

    Returns
    -------
    v_raw_score : np.ndarray
        Raw score of shape (n_cell,).
    v_score_weight : np.ndarray
        Gene weights of shape (n_disease_gene,).

    Notes
    -----

    """

    gene_list = list(gene_list)
    gene_weight = list(gene_weight)

    assert weight_opt in {"uniform", "vs", "inv_std", "od"}, (
        "weight_opt=%s is not one of {'uniform', 'vs', 'inv_std', 'od}'" % weight_opt
    )

    # Compute overdispersion score
    # (used only for benchmarking, do not support implicit covariate correction mode)
    if weight_opt == "od":
        return _compute_overdispersion_score(adata, gene_list, gene_weight)

    # Compute other weighted average scores
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `scdrs.pp.preprocess` first"

    df_gene = adata.uns["SCDRS_PARAM"]["GENE_STATS"]
    flag_sparse = adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"]
    flag_cov = adata.uns["SCDRS_PARAM"]["FLAG_COV"]

    if weight_opt == "uniform":
        v_score_weight = np.ones(len(gene_list))
    if weight_opt == "vs":
        v_score_weight = 1 / np.sqrt(df_gene.loc[gene_list, "var_tech"].values + 1e-2)
    if weight_opt == "inv_std":
        v_score_weight = 1 / np.sqrt(df_gene.loc[gene_list, "var"].values + 1e-2)

    if gene_weight is not None:
        v_score_weight = v_score_weight * np.array(gene_weight)
    v_score_weight = v_score_weight / v_score_weight.sum()

    if flag_sparse and flag_cov:
        # Implicit covariate correction mode
        cell_list = list(adata.obs_names)
        cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])
        cov_mat = adata.uns["SCDRS_PARAM"]["COV_MAT"].loc[cell_list, cov_list].values
        cov_beta = (
            adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values.T
        )
        gene_mean = adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

        # Compute v_raw_score = transformed_X @ v_score_weight
        # where transformed_X = adata.X + cov_mat @ cov_beta + gene_mean
        v_raw_score = (
            adata[:, gene_list].X.dot(v_score_weight)
            + cov_mat @ (cov_beta @ v_score_weight)
            + gene_mean @ v_score_weight
        ).flatten()
    else:
        # Normal mode
        v_raw_score = adata[:, gene_list].X.dot(v_score_weight).reshape([-1])

    return v_raw_score, v_score_weight


def _compute_overdispersion_score(adata, gene_list, gene_weight):
    """Compute overdispersion score

        Raw weight: w_g_raw = gene_weight / \sigma_{tech,g}^2

        Normalized weight: w_g = w_g_raw / \sum_g w_g_raw

        Overdispersion score: s_c = \sum_g w_g * [(X_cg - \mu_g)^2 - \sigma_{tech,g}^2]

    Args
    ----
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    gene_list : list
        Disease gene list of length n_disease_gene.
    gene_weight : list
        Gene weights of length n_disease_gene for genes in the gene_list.

    Returns
    -------
    v_raw_score : np.ndarray
        Raw score of shape (n_cell,).
    v_score_weight : np.ndarray
        Gene weights of shape (n_disease_gene,).
    """

    gene_list = list(gene_list)
    gene_weight = list(gene_weight)

    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `scdrs.pp.preprocess` first"

    # Mode check
    flag_sparse = adata.uns["SCDRS_PARAM"]["FLAG_SPARSE"]
    flag_cov = adata.uns["SCDRS_PARAM"]["FLAG_COV"]
    if flag_sparse and flag_cov:
        cell_list = list(adata.obs_names)
        cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])
        mat_X = (
            adata[:, gene_list].X.toarray()
            + adata.uns["SCDRS_PARAM"]["COV_MAT"]
            .loc[cell_list, cov_list]
            .values.dot(
                adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values.T
            )
            + adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values
        )
    else:
        mat_X = adata[:, gene_list].X

    v_mean = adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "mean"].values
    v_var_tech = (
        adata.uns["SCDRS_PARAM"]["GENE_STATS"].loc[gene_list, "var_tech"].values
    )

    v_w = 1 / (v_var_tech + 1e-2)
    if gene_weight is not None:
        v_w = v_w * np.array(gene_weight)
    v_w = v_w / v_w.sum()

    # Compute overdispersion score
    if sp.sparse.issparse(mat_X):
        v_raw_score = mat_X.power(2).dot(v_w).reshape([-1])  # Quadratic term
    else:
        v_raw_score = (mat_X ** 2).dot(v_w).reshape([-1])  # Quadratic term
    v_raw_score = v_raw_score - mat_X.dot(2 * v_w * v_mean).reshape([-1])  # Linear term
    v_raw_score = (
        v_raw_score + (v_w * (v_mean ** 2 - v_var_tech)).sum()
    )  # Constant term

    return v_raw_score, np.ones(len(gene_list))


def _correct_background(
    v_raw_score, mat_ctrl_raw_score, v_var_ratio_c2t, save_intermediate=None
):
    """Cell-wise and gene-wise background correction

    Args
    ----
    v_raw_score : np.ndarray
        Disease raw score of shape (n_cell,).
    mat_ctrl_raw_score : np.ndarray
        Disease raw control scores of shape (n_cell,n_ctrl).
    v_var_ratio_c2t : np.ndarray
        Ratio of independent variance between control scores and disease score,
        of shape (n_ctrl).
    save_intermediate : str
        File path prefix for saving intermediate results.

    Returns
    -------
    v_norm_score : np.ndarray
        Normalized disease score of shape (n_cell,)
    mat_ctrl_norm_score : np.ndarray
        Normalized control scores of shape (n_cell,n_ctrl).
    """

    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.tsv.gz",
            v_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.tsv.gz",
            mat_ctrl_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Zero-values are assigned the smallest values at the end
    ind_zero_score = v_raw_score == 0
    ind_zero_ctrl_score = mat_ctrl_raw_score == 0

    # First gene set alignment: mean 0 and same independent variance
    v_raw_score = v_raw_score - v_raw_score.mean()
    mat_ctrl_raw_score = mat_ctrl_raw_score - mat_ctrl_raw_score.mean(axis=0)
    mat_ctrl_raw_score = mat_ctrl_raw_score / np.sqrt(v_var_ratio_c2t)
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.1st_gs_alignment.tsv.gz",
            v_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.1st_gs_alignment.tsv.gz",
            mat_ctrl_raw_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Cell-wise standardization
    v_mean = mat_ctrl_raw_score.mean(axis=1)
    v_std = mat_ctrl_raw_score.std(axis=1)
    v_norm_score = v_raw_score.copy()
    v_norm_score = (v_norm_score - v_mean) / v_std
    mat_ctrl_norm_score = ((mat_ctrl_raw_score.T - v_mean) / v_std).T
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.cellwise_standardization.tsv.gz",
            v_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.cellwise_standardization.tsv.gz",
            mat_ctrl_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Second gene set alignment: mean 0
    v_norm_score = v_norm_score - v_norm_score.mean()
    mat_ctrl_norm_score = mat_ctrl_norm_score - mat_ctrl_norm_score.mean(axis=0)
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.2nd_gs_alignment.tsv.gz",
            v_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.2nd_gs_alignment.tsv.gz",
            mat_ctrl_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )

    # Set cells with raw_score=0 to the minimum norm_score value
    norm_score_min = min(v_norm_score.min(), mat_ctrl_norm_score.min())
    v_norm_score[ind_zero_score] = norm_score_min - 1e-3
    mat_ctrl_norm_score[ind_zero_ctrl_score] = norm_score_min
    if save_intermediate is not None:
        np.savetxt(
            save_intermediate + ".raw_score.final.tsv.gz",
            v_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )
        np.savetxt(
            save_intermediate + ".ctrl_raw_score.final.tsv.gz",
            mat_ctrl_norm_score,
            fmt="%.9e",
            delimiter="\t",
        )

    return v_norm_score, mat_ctrl_norm_score


def _get_p_from_empi_null(v_t, v_t_null):
    """Compute p-value from empirical null
    For score T and a set of null score T_1,...T_N, the p-value is

        p= [1 + \Sigma_{i=1}^N 1_{ (T_i \geq T) }] / (1+N)

    If T, T_1, ..., T_N are i.i.d. variables following a null distritbuion,
    then p is super-uniform.

    The naive algorithm is N^2. Here we provide an O(N log N) algorithm to
    compute the p-value for each of the N elements in v_t

    Args
    ----
    v_t : np.ndarray
        Observed score of shape (M,).
    v_t_null : np.ndarray
        Null scores of shape (N,).

    Returns
    -------
    v_p: : np.ndarray
        P-value for each element in v_t of shape (M,).
    """

    v_t = np.array(v_t)
    v_t_null = np.array(v_t_null)

    v_t_null = np.sort(v_t_null)
    v_pos = np.searchsorted(v_t_null, v_t, side="left")
    v_p = (v_t_null.shape[0] - v_pos + 1) / (v_t_null.shape[0] + 1)
    return v_p


##############################################################################
######################### Code for comparison methods ########################
##############################################################################
def score_cell_vision(adata, gene_list):

    """Score cells based on the trait gene set

    Args
    ----
    data (n_cell, n_gene) : AnnData
        data.X should contain size-normalized log1p transformed count data
    gene_list (n_disease_gene) : list
        Trait gene list

    Returns
    -------
    df_res (n_cell, n_key) : pd.DataFrame (dtype=np.float32)
        Columns:
        1. score: Vision signature score
        2. pval: p-value computed from the Vision score
    """

    gene_list = sorted(set(gene_list) & set(adata.var_names))
    v_mean, v_var = scdrs.pp._get_mean_var(adata.X, axis=1)

    v_score = adata[:, gene_list].X.mean(axis=1)
    v_score = np.array(v_score).reshape([-1])
    v_score = (v_score - v_mean) / np.sqrt(v_var / len(gene_list))
    v_p = 1 - sp.stats.norm.cdf(v_score)

    v_norm_score = (v_score - v_score.mean()) / v_score.std()
    v_norm_p = 1 - sp.stats.norm.cdf(v_norm_score)

    dic_res = {
        "score": v_score,
        "pval": v_p,
        "norm_score": v_norm_score,
        "norm_pval": v_norm_p,
    }
    df_res = pd.DataFrame(index=adata.obs.index, data=dic_res, dtype=np.float32)
    return df_res


def score_cell_scanpy(adata, gene_list):

    """Score cells based on the trait gene set

    Args
    ----
        data (n_cell, n_gene) : AnnData
            data.X should contain size-normalized log1p transformed count data
        gene_list (n_disease_gene) : list
            Trait gene list

    Returns
    -------
        df_res (n_cell, n_key) : pd.DataFrame (dtype=np.float32)
            Columns:
            1. score: Vision signature score
            2. pval: p-value computed from the Vision score
    """

    gene_list = sorted(set(gene_list) & set(adata.var_names))
    sc.tl.score_genes(adata, gene_list=gene_list)

    v_score = adata.obs["score"]
    v_score_z = (v_score - v_score.mean()) / np.sqrt(v_score.var())
    v_p = 1 - sp.stats.norm.cdf(v_score_z)

    dic_res = {"score": v_score, "score_z": v_score_z, "pval": v_p}
    df_res = pd.DataFrame(index=adata.obs.index, data=dic_res, dtype=np.float32)
    return df_res


##############################################################################
######################## Code for downstream analysis ########################
##############################################################################
def downstream_group_analysis(
    adata: anndata.AnnData,
    df_full_score: pd.DataFrame,
    group_cols: List[str],
    fdr_thresholds: List[float] = [0.05, 0.1, 0.2],
) -> Dict[str, pd.DataFrame]:
    """
    scDRS group-level analysis.

    For each annotation in `group_cols` and each group of cells in the annotation, compute:
    
    1. Proportion of FDR < 0.1 cells.
    2. Group-level trait association.
    3. Group-level heterogeneity.

    `connectivities` is expected in `adata.obsp` for the group-level heterogeneity analysis.
    Recommended parameters: `sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20)`.

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    df_full_score : pd.DataFrame
        scDRS `.full_score` file for a given trait.
    group_cols : list of str
        List of column names in adata.obs used to define cell groups.

    Returns
    -------    
    dict_df_res : Dict[str, pd.DataFrame]
        Group-level statistics (n_group, n_stats) keyed by the group names.
    """

    assert (
        "connectivities" in adata.obsp
    ), "Expect `connectivities` in `adata.obsp`; run `sc.pp.neighbors` first"

    assert (
        len(set(group_cols) - set(adata.obs)) == 0
    ), "Missing `group_cols` variables from `adata.obs.columns`."

    # Align cells between `adata` and `df_full_score`.
    cell_list = sorted(set(adata.obs_names) & set(df_full_score.index))
    control_list = [x for x in df_full_score.columns if x.startswith("ctrl_norm_score")]
    n_ctrl = len(control_list)
    df_reg = adata.obs.loc[cell_list, group_cols].copy()
    df_reg = df_reg.join(
        df_full_score.loc[cell_list, ["norm_score"] + control_list + ["pval"]]
    )

    # Group-level analysis; dict_df_res : group_col -> df_res
    dict_df_res = {}
    for group_col in group_cols:
        group_list = sorted(set(adata.obs[group_col]))
        res_cols = [
            "n_cell",
            "n_ctrl",
            "assoc_mcp",
            "assoc_mcz",
            "hetero_mcp",
            "hetero_mcz",
        ]
        for fdr_threshold in fdr_thresholds:
            res_cols.append(f"n_fdr_{fdr_threshold}")

        df_res = pd.DataFrame(index=group_list, columns=res_cols, dtype=np.float32)

        df_fdr = pd.DataFrame(
            {"fdr": multipletests(df_reg["pval"].values, method="fdr_bh")[1]},
            index=df_reg.index,
        )

        for group in group_list:
            group_cell_list = list(df_reg.index[df_reg[group_col] == group])
            # Basic info
            df_res.loc[group, ["n_cell", "n_ctrl"]] = [len(group_cell_list), n_ctrl]

            # Number of FDR < fdr_threshold cells in each group
            for fdr_threshold in fdr_thresholds:
                df_res.loc[group, f"n_fdr_{fdr_threshold}"] = (
                    df_fdr.loc[group_cell_list, "fdr"].values < fdr_threshold
                ).sum()

        # Association
        for group in group_list:
            group_cell_list = list(df_reg.index[df_reg[group_col] == group])
            score_q95 = np.quantile(df_reg.loc[group_cell_list, "norm_score"], 0.95)
            v_ctrl_score_q95 = np.quantile(
                df_reg.loc[group_cell_list, control_list], 0.95, axis=0
            )
            mc_p = ((v_ctrl_score_q95 >= score_q95).sum() + 1) / (
                v_ctrl_score_q95.shape[0] + 1
            )
            mc_z = (score_q95 - v_ctrl_score_q95.mean()) / v_ctrl_score_q95.std()
            df_res.loc[group, ["assoc_mcp", "assoc_mcz"]] = [mc_p, mc_z]

        # Heterogeneity
        df_rls = test_gearysc(
            adata[cell_list], df_reg.loc[cell_list, :], groupby=group_col
        )
        for ct in group_list:
            mc_p, mc_z = df_rls.loc[ct, ["pval", "zsc"]]
            df_res.loc[ct, ["hetero_mcp", "hetero_mcz"]] = [mc_p, mc_z]

        # write to dict for this group_col
        dict_df_res[group_col] = df_res

    return dict_df_res


def downstream_corr_analysis(
    adata: anndata.AnnData, df_full_score: pd.DataFrame, var_cols: List[str]
) -> pd.DataFrame:
    """
    scDRS cell-level correlation analysis.

    For a given individual cell-level annotation (e.g., T cell effectorness gradient),
    assess association between disease and the individual cell-level variable
    (control-score-based Monte Carlo tests using Pearson's correlation).

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    df_full_score : pd.DataFrame
        scDRS `.full_score` file for a given trait.
    var_cols : List[str]
        List of column names in `adata.obs` for continous cell-level variables.

    Returns
    -------    
    df_res : pd.DataFrame
        Correlation results (n_var, n_stats).
    """

    assert (
        len(set(var_cols) - set(adata.obs)) == 0
    ), "Missing `var_cols` variables from `adata.obs.columns`."

    cell_list = sorted(set(adata.obs_names) & set(df_full_score.index))
    control_list = [x for x in df_full_score.columns if x.startswith("ctrl_norm_score")]
    n_ctrl = len(control_list)
    df_reg = adata.obs.loc[cell_list, var_cols].copy()
    df_reg = df_reg.join(df_full_score.loc[cell_list, ["norm_score"] + control_list])

    # Variable-disease correlation
    col_list = ["n_ctrl", "corr_mcp", "corr_mcz"]
    df_res = pd.DataFrame(index=var_cols, columns=col_list, dtype=np.float32)
    for var_col in var_cols:
        corr_ = np.corrcoef(df_reg[var_col], df_reg["norm_score"])[0, 1]
        v_corr_ = np.array(
            [np.corrcoef(df_reg[var_col], df_reg[x])[0, 1] for x in control_list]
        )
        mc_p = ((v_corr_ >= corr_).sum() + 1) / (v_corr_.shape[0] + 1)
        mc_z = (corr_ - v_corr_.mean()) / v_corr_.std()
        df_res.loc[var_col] = [n_ctrl, mc_p, mc_z]

    return df_res


def downstream_gene_analysis(
    adata: anndata.AnnData, df_full_score: pd.DataFrame
) -> pd.DataFrame:
    """
    scDRS gene-level correlation analysis.

    Compute correlation between each gene and the scDRS disease score.

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    df_full_score : pd.DataFrame
        scDRS `.full_score` file for a given trait.

    Returns
    -------    
    df_res : pd.DataFrame
        Correlation results (n_gene, n_stats).
    """

    cell_list = sorted(set(adata.obs_names) & set(df_full_score.index))
    control_list = [x for x in df_full_score.columns if x.startswith("ctrl_norm_score")]
    df_reg = df_full_score.loc[cell_list, ["norm_score"]]

    mat_expr = adata[cell_list].X.copy()
    v_corr = _pearson_corr(mat_expr, df_reg["norm_score"].values)
    df_res = pd.DataFrame(
        index=adata.var_names, columns=["CORR", "RANK"], dtype=np.float32
    )
    df_res["CORR"] = v_corr
    df_res.sort_values("CORR", ascending=False, inplace=True)
    df_res["RANK"] = np.arange(df_res.shape[0])
    return df_res


##############################################################################
##################### Subroutines for downstream analysis ####################
##############################################################################
def test_gearysc(
    adata: anndata.AnnData,
    df_full_score: pd.DataFrame,
    groupby: str,
    opt="control_distribution_match",
) -> pd.DataFrame:
    """
    Compute significance level for Geary's C statistics.

    Parameters
    ----------
    adata : anndata.AnnData
        Must contain `connectivities` to compute the Geary's C statistic.
    df_full_score : DataFrame
        DataFrame with the scores of the cells, contains
        columns `zscore`, `norm_score`, `ctrl_norm_score_{i}`
    groupby : str
        Column name of the groupby variable.
    opt : str
        Options:
            - "control_distribution_match":
                The distribution of the scores of the control scores is similar to
                the distribution of the scores of the disease scores.

    Returns
    -------    
    df_rls : DataFrame
        DataFrame with the results of the test with `n_group` rows and 4 columns:
        
        - `pval`: significance level of Geary's C
        - `trait`: Geary's C test statistic of the trait scores
        - `ctrl_mean`: mean of the control scores
        - `ctrl_sd`: standard deviation of the control scores
    """
    assert np.all(
        df_full_score.index == adata.obs_names
    ), "adata.obs_names must match df_full_score.index"
    norm_score = df_full_score["norm_score"]
    ctrl_norm_score = df_full_score[
        [col for col in df_full_score.columns if col.startswith(f"ctrl_norm_score_")]
    ]
    n_null = ctrl_norm_score.shape[1]
    df_meta = adata.obs.copy()
    df_stats = pd.DataFrame(
        index=df_meta[groupby].unique(),
        columns=["trait"] + [f"null_{i_null}" for i_null in range(n_null)],
        data=np.nan,
    )

    for group, df_group in df_meta.groupby(groupby):
        group_index = df_group.index
        group_adata = adata[group_index]
        group_norm_score = norm_score[group_index]
        group_ctrl_norm_score = ctrl_norm_score.loc[group_index, :]

        if opt == "control_distribution_match":
            # control distribution match
            from scipy.stats import rankdata

            def distribution_match(v, ref):
                """
                Use order in `v` to match the distribution of `ref`
                """
                return np.sort(ref)[rankdata(v, method="ordinal") - 1]

            df_stats.loc[group, "trait"] = gearys_c(
                group_adata, group_norm_score.values
            )

            for i_null in range(n_null):
                df_stats.loc[group, f"null_{i_null}"] = gearys_c(
                    group_adata,
                    distribution_match(
                        group_ctrl_norm_score.iloc[:, i_null].values,
                        ref=group_norm_score,
                    ),
                )

        elif opt == "permutation":
            # permutation
            df_stats.loc[group, "trait"] = gearys_c(
                group_adata, group_norm_score.values
            )
            for i_null in range(n_null):
                df_stats.loc[group, f"null_{i_null}"] = gearys_c(
                    group_adata, np.random.permutation(group_norm_score.values)
                )
        elif opt == "control":
            # control
            df_stats.loc[group, "trait"] = gearys_c(
                group_adata, group_norm_score.values
            )
            for i_null in range(n_null):
                df_stats.loc[group, f"null_{i_null}"] = gearys_c(
                    group_adata, group_ctrl_norm_score.iloc[:, i_null].values
                )
        else:
            raise NotImplementedError

    # Summarize
    trait_col = "trait"
    ctrl_cols = [col for col in df_stats.columns if col.startswith("null_")]
    pval = (
        (df_stats[trait_col].values > df_stats[ctrl_cols].values.T).sum(axis=0) + 1
    ) / (len(ctrl_cols) + 1)
    pval[np.isnan(df_stats[trait_col])] = np.nan

    df_rls = pd.DataFrame(
        {
            "pval": pval,
            "trait": df_stats[trait_col].values,
            "ctrl_mean": df_stats[ctrl_cols].mean(axis=1).values,
            "ctrl_std": df_stats[ctrl_cols].std(axis=1).values,
        },
        index=df_stats.index,
    )

    df_rls["zsc"] = (
        -(df_rls[trait_col].values - df_rls["ctrl_mean"]) / df_rls["ctrl_std"]
    )
    return df_rls


def gearys_c(adata, vals):
    """
    Compute Geary's C statistics for an AnnData.
    
    Adopted from https://github.com/ivirshup/scanpy/blob/metrics/scanpy/metrics/_gearys_c.py
        
    :math:`C=\\frac{(N - 1)\\sum_{i,j} w_{i,j} (x_i - x_j)^2}{2W \\sum_i (x_i - \\bar{x})^2}`

    Parameters
    ----------
    adata : AnnData object
        adata.obsp["Connectivities] should contain the connectivity graph,
        with shape (n_obs, n_obs).
    vals : array-like
        Values to calculate Geary's C for. If one dimensional, should have
        shape (n_obs,).

    Returns
    -------    
    C : float
        Geary's C statistics.
    """
    graph = adata.obsp["connectivities"]
    assert graph.shape[0] == graph.shape[1]
    graph_data = graph.data.astype(np.float_, copy=False)
    assert graph.shape[0] == vals.shape[0]
    assert np.ndim(vals) == 1

    W = graph_data.sum()
    N = len(graph.indptr) - 1
    vals_bar = vals.mean()
    vals = vals.astype(np.float_)

    # numerators
    total = 0.0
    for i in range(N):
        s = slice(graph.indptr[i], graph.indptr[i + 1])
        # indices of corresponding neighbors
        i_indices = graph.indices[s]
        # corresponding connecting weights
        i_data = graph_data[s]
        total += np.sum(i_data * ((vals[i] - vals[i_indices]) ** 2))

    numer = (N - 1) * total
    denom = 2 * W * ((vals - vals_bar) ** 2).sum()
    C = numer / denom

    return C


def _pearson_corr(mat_X, mat_Y):
    """Pearson's correlation between every columns in mat_X and mat_Y.

    Parameters
    ----------
    mat_X : np.ndarray
        First matrix of shape (N,M1).
    mat_Y : np.ndarray
        Second matrix of shape (N,M2).

    Returns
    -------    
    mat_corr : np.ndarray
        Correlation matrix of shape (M1,M2).
    """
    # If sparse, use _pearson_corr_sparse
    if sp.sparse.issparse(mat_X) | sp.sparse.issparse(mat_Y):
        return _pearson_corr_sparse(mat_X, mat_Y)

    # Reshape
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])

    mat_X = (mat_X - mat_X.mean(axis=0)) / mat_X.std(axis=0).clip(min=1e-8)
    mat_Y = (mat_Y - mat_Y.mean(axis=0)) / mat_Y.std(axis=0).clip(min=1e-8)
    mat_corr = mat_X.T.dot(mat_Y) / mat_X.shape[0]
    mat_corr = np.array(mat_corr, dtype=np.float32)

    if (mat_X.shape[1] == 1) | (mat_Y.shape[1] == 1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr


def _pearson_corr_sparse(mat_X, mat_Y):
    """Pearson's correlation between every columns in mat_X and mat_Y (sparse matrix)

    Parameters
    ----------
    mat_X : np.ndarray
        First matrix of shape (N,M1).
    mat_Y : np.ndarray
        Second matrix of shape (N,M2).

    Returns
    -------    
    mat_corr : np.ndarray
        Correlation matrix of shape (M1,M2).
    """

    # Reshape
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])

    # Convert to sparse matrix if not already sparse
    if sp.sparse.issparse(mat_X) is False:
        mat_X = sp.sparse.csr_matrix(mat_X)
    if sp.sparse.issparse(mat_Y) is False:
        mat_Y = sp.sparse.csr_matrix(mat_Y)

    # Compute v_mean,v_var
    v_X_mean, v_X_var = scdrs.pp._get_mean_var(mat_X, axis=0)
    v_X_sd = np.sqrt(v_X_var).clip(min=1e-8)
    v_Y_mean, v_Y_var = scdrs.pp._get_mean_var(mat_Y, axis=0)
    v_Y_sd = np.sqrt(v_Y_var).clip(min=1e-8)

    mat_corr = mat_X.T.dot(mat_Y) / mat_X.shape[0]
    mat_corr = mat_corr - v_X_mean.reshape([-1, 1]).dot(v_Y_mean.reshape([1, -1]))
    mat_corr = mat_corr / v_X_sd.reshape([-1, 1]).dot(v_Y_sd.reshape([1, -1]))
    mat_corr = np.array(mat_corr, dtype=np.float32)

    if (mat_X.shape[1] == 1) | (mat_Y.shape[1] == 1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr


##############################################################################
############################## Archived functions ############################
##############################################################################
def correlate_gene(
    data, trs_name="trs_ez", suffix="", corr_opt="pearson", cov_list=None, copy=False
):

    """Compute the correlation between gene expressions and TRS

    Args
    ----
    data (n_cell, n_gene) : AnnData
        adata.X should contain size-normalized log1p transformed count data
    trs_name : str
        The variable to correlate gene expression with. Should be one column in data.obs.
    suffix : str
        The name of the added gene-wise correlation would be 'trs_corr'+suffix.
    corr_opt : str
        Option for computing the correlation
        'pearson': Pearson's correlation
        'spearman': Spearman's correlation
    cov_list : list of str
        Covariates to control for.
        The covariates are first centered and then regressed out from
            both trs_name and the gene expression before computing the correlation.
        Elements in cov_list should be present in data.obs.columns
    copy : bool
        If to make copy of the AnnData object

    Returns
    -------
    adata (AnnData):
        Add the columns 'scdrs_corr'+suffix to data.var
    """

    adata = data.copy() if copy else data

    # Check options
    corr_opt_list = ["pearson", "spearman"]
    if corr_opt not in corr_opt_list:
        raise ValueError(
            "# compute_scdrs_corr: corr_opt not in [%s]"
            % ", ".join([str(x) for x in corr_opt_list])
        )
    if trs_name not in adata.obs.columns:
        raise ValueError("# compute_scdrs_corr: %s not in data.obs.columns" % trs_name)
    if cov_list is not None:
        temp_list = list(set(cov_list) - set(adata.obs.columns))
        if len(temp_list) > 0:
            raise ValueError(
                "# compute_scdrs_corr: covariates %s not in data.obs.columns"
                % ",".join(temp_list)
            )

    # Get data
    mat_X = data.X.toarray()
    v_trs = data.obs[trs_name].values.copy()

    # Regress out covariates
    if cov_list is not None:
        mat_cov = adata.obs[cov_list].values.copy()
        mat_cov = mat_cov - mat_cov.mean(axis=0)
        v_trs = scdrs.pp.reg_out(v_trs, mat_cov)
        mat_X = scdrs.pp.reg_out(mat_X, mat_cov)

    # Compute correlation
    if corr_opt == "pearson":
        v_corr = _pearson_corr(mat_X, v_trs)

    if corr_opt == "spearman":
        v_corr = _spearman_corr(mat_X, v_trs)

    adata.var["scdrs_corr" + suffix] = v_corr

    return adata if copy else None


def _spearman_corr(mat_X, mat_Y):
    """Spearman's correlation between every columns in mat_X and mat_Y

    Args
    ----
    mat_X (N,M1): np.ndarray
    mat_Y (N,M2): np.ndarray

    Returns
    -------
    mat_corr (M1,M2): np.ndarray
        Correlation matrix
    """

    # Reshape
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])

    mat_X = _get_rank(mat_X, axis=0)
    mat_Y = _get_rank(mat_Y, axis=0)

    mat_X = (mat_X - mat_X.mean(axis=0)) / mat_X.std(axis=0).clip(min=1e-8)
    mat_Y = (mat_Y - mat_Y.mean(axis=0)) / mat_Y.std(axis=0).clip(min=1e-8)
    mat_corr = mat_X.T.dot(mat_Y) / mat_X.shape[0]

    if (mat_X.shape[1] == 1) | (mat_Y.shape[1] == 1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr


def _get_rank(mat_X, axis=0):
    """Get rank for each row/columns of the given matrix

    Args
    ----
    mat_X (N,M): np.ndarray
    axis: int
        axis=0: column-wise rank (across rows)
        axis=1: row-wise rank (across columns)
    Returns
    -------
    mat_rank  (N,M): np.ndarray
        Rank matrix
    """

    if axis == 0:
        mat_X = np.argsort(mat_X, axis=0)
        mat_rank = np.empty_like(mat_X)
        temp_v = np.arange(mat_X.shape[0])
        for i_col in range(mat_X.shape[1]):
            mat_rank[mat_X[:, i_col], i_col] = temp_v

    if axis == 1:
        mat_X = np.argsort(mat_X, axis=1)
        mat_rank = np.empty_like(mat_X)
        temp_v = np.arange(mat_X.shape[1])
        for i_row in range(mat_X.shape[0]):
            mat_rank[i_row, mat_X[i_row, :]] = temp_v

    return mat_rank


# @Kangcheng: this is the archived section.
# Remove it if you want but it can stay here for a while.
def group_stats(
    df_drs: pd.DataFrame,
    adata: anndata.AnnData,
    group_col: str,
    stats=["assoc", "hetero"],
) -> pd.DataFrame:
    """compute group-level statistics for scDRS results

    Parameters
    ----------
    df_drs : pd.DataFrame
        scDRS results dataframe
    adata : anndata.AnnData
        AnnData object
    group_col : str
        Column name of group column in adata.obs
    stats : list, optional
        statistics to compute, by default ["assoc", "hetero"]

    Returns
    -------
    pd.DataFrame
        Group-level statistics (n_group x n_stats)
    """

    assert group_col in adata.obs.columns, "group_col not in adata.obs"
    # stats must be contained in ["fdr_prop", "assoc", "hetero"]
    assert set(stats).issubset(
        set(["assoc", "hetero"])
    ), "stats must be contained in ['assoc', 'hetero']"
    group_list = np.unique(adata.obs[group_col])
    col_list = ["n_cell", "n_ctrl", "fdr_prop"]

    df = adata.obs[[group_col]].join(df_drs)
    for stat in stats:
        col_list.extend([stat + "_pval", stat + "_zsc"])
    df_stats = pd.DataFrame(index=group_list, columns=col_list, dtype=float)

    norm_score = df["norm_score"].values
    ctrl_norm_score = df[
        [col for col in df.columns if col.startswith(f"ctrl_norm_score_")]
    ].values
    n_ctrl = ctrl_norm_score.shape[1]

    v_fdr = multipletests(df["pval"].values, method="fdr_bh")[1]

    for group in group_list:
        # Basic info
        group_index = df[group_col].values == group
        df_stats.loc[group, "n_cell"] = sum(group_index)
        df_stats.loc[group, "n_ctrl"] = n_ctrl

        # FDR proportion
        df_stats.loc[group, "fdr_prop"] = (v_fdr[group_index] < 0.1).mean()

        # association
        if "assoc" in stats:
            score_q95 = np.quantile(norm_score[group_index], 0.95)
            v_ctrl_score_q95 = np.quantile(
                ctrl_norm_score[group_index, :], 0.95, axis=0
            )
            mc_p = ((v_ctrl_score_q95 >= score_q95).sum() + 1) / (
                v_ctrl_score_q95.shape[0] + 1
            )
            mc_z = (score_q95 - v_ctrl_score_q95.mean()) / v_ctrl_score_q95.std()
            df_stats.loc[group, ["assoc_pval", "assoc_zsc"]] = [mc_p, mc_z]

    # Heterogeneity
    if "hetero" in stats:
        df_rls = test_gearysc(adata, df, groupby=group_col)
        for group in group_list:
            mc_p, mc_z = df_rls.loc[group, ["pval", "zsc"]]
            df_stats.loc[group, ["hetero_pval", "hetero_zsc"]] = [mc_p, mc_z]

    return df_stats
