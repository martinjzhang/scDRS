import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
from skmisc.loess import loess
from tqdm import tqdm
import scdrs.pp as pp


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

    It operates in implicit-covariate-correction mode if both `FLAG_SPARSE`z
    and `FLAG_COV` are `True`, where computations are based on the implicit
    covariate-corrected data `CORRECTED_X = data.X + COV_MAT * COV_BETA + COV_GENE_MEAN`.

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
        Number of control gene sets
    n_genebin : int, default=200
        Number of bins for dividing genes by ctrl_match_key if
        `data.uns["SCDRS_PARAM"]["GENE_STATS"][ctrl_match_key]` is a continuous variable
    weight_opt : str, default="vs"
        Option for computing the raw score
        - 'uniform': average over the genes in the gene_list.
        - 'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct).
        - 'inv_std': weighted average with weights equal to 1/std.
        - 'od': overdispersion score.
    copy : bool, default=False
        If to make copy of the AnnData object to avoid writing on the orignal data
    return_raw_ctrl_score : bool, default=False
        If to return control scores.
    return_norm_ctrl_score : bool, default=False
        If to return control scores.
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
    if weight_opt in ["uniform", "vs", "inv_std"]:
        # For raw scores compuated as weighted average. estimate variance ratio assuming independence
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

    """Subroutine for score_cell. Select control gene sets.

    Args
    ----
    input_df_gene : pd.DataFrame
        Gene-level statistics of shape (n_gene, ):.
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
        cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])

        cov_mat = adata.uns["SCDRS_PARAM"]["COV_MAT"].values
        cov_beta = adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values
        gene_mean = adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list, :].values

        # Compute v_raw_score = transformed_X @ v_score_weight
        # where transformed_X = adata.X + cov_mat @ cov_beta + gene_mean
        v_raw_score = (
            adata[:, gene_list].X.dot(v_score_weight)
            + cov_mat @ (cov_beta.T @ v_score_weight)
            + gene_mean.T @ v_score_weight
        ).flatten()
    else:
        # Normal mode
        v_raw_score = adata[:, gene_list].X.dot(v_score_weight).reshape([-1])

    return v_raw_score, v_score_weight


# fixit: review
def _compute_overdispersion_score(adata, gene_list, gene_weight):
    """Compute overdispersion score

        Let w_g_raw = gene_weight / \sigma_{tech,g}^2
        Let w_g = w_g_raw / \sum_g w_g_raw
        s_c = \sum_g w_g * [(X_cg - \mu_g)^2 - \sigma_{tech,g}^2]

    Args
    ----
    adata (n_cell, n_gene) : AnnData
        adata.X should contain size-normalized log1p transformed count data
    gene_list (n_disease_gene) : list
        Trait gene list
    gene_weight (n_disease_gene) : list/np.ndarray
        Gene weights for genes in the gene_list

    Returns
    -------
    v_raw_score (n_cell,) : np.ndarray
        Raw score
    v_score_weight (n_disease_gene,) : np.ndarray
        Gene weights score
    """

    v_mean = adata.var.loc[gene_list, "mean"].values
    v_var_tech = adata.var.loc[gene_list, "var_tech"].values

    v_w = 1 / (v_var_tech + 1e-2)
    if gene_weight is not None:
        v_w = v_w * np.array(gene_weight)
    v_w = v_w / v_w.sum()

    # Compute overdispersion score
    if sp.sparse.issparse(adata.X):
        v_raw_score = (
            adata[:, gene_list].X.power(2).dot(v_w).reshape([-1])
        )  # Quadratic term
    else:
        v_raw_score = (
            (adata[:, gene_list].X ** 2).dot(v_w).reshape([-1])
        )  # Quadratic term
    v_raw_score = v_raw_score - adata[:, gene_list].X.dot(2 * v_w * v_mean).reshape(
        [-1]
    )  # Linear term
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

    # Calibrate gene sets (mean 0 and same independent variance)
    ind_zero_score = (v_raw_score == 0)
    ind_zero_ctrl_score = (mat_ctrl_raw_score == 0)

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

    # Cell-wise correction
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

    # Gene-set-wise correction
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
            save_intermediate
            + ".ctrl_raw_score.2nd_gs_alignment.tsv.gz",
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
    v_mean, v_var = pp._get_mean_var(adata.X, axis=1)

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
        v_trs = pp.reg_out(v_trs, mat_cov)
        mat_X = pp.reg_out(mat_X, mat_cov)

    # Compute correlation
    if corr_opt == "pearson":
        v_corr = _pearson_corr(mat_X, v_trs)

    if corr_opt == "spearman":
        v_corr = _spearman_corr(mat_X, v_trs)

    adata.var["scdrs_corr" + suffix] = v_corr

    return adata if copy else None


def _pearson_corr(mat_X, mat_Y):
    """Pearson's correlation between every columns in mat_X and mat_Y

    Args
    ----
    mat_X (N,M1): np.ndarray
    mat_Y (N,M2): np.ndarray

    Returns
    -------
    mat_corr: (M1,M2): np.ndarray
        Correlation matrix
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

    Args
    ----
    mat_X (N,M1): sp.sparse
    mat_Y (N,M2): sp.sparse

    Returns
    -------
    mat_corr: (M1,M2): np.ndarray
        Correlation matrix
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
    v_X_mean, v_X_var = pp._get_mean_var(mat_X, axis=0)
    v_X_sd = np.sqrt(v_X_var).clip(min=1e-8)
    v_Y_mean, v_Y_var = pp._get_mean_var(mat_Y, axis=0)
    v_Y_sd = np.sqrt(v_Y_var).clip(min=1e-8)

    mat_corr = mat_X.T.dot(mat_Y) / mat_X.shape[0]
    mat_corr = mat_corr - v_X_mean.reshape([-1, 1]).dot(v_Y_mean.reshape([1, -1]))
    mat_corr = mat_corr / v_X_sd.reshape([-1, 1]).dot(v_Y_sd.reshape([1, -1]))
    mat_corr = np.array(mat_corr, dtype=np.float32)

    if (mat_X.shape[1] == 1) | (mat_Y.shape[1] == 1):
        return mat_corr.reshape([-1])
    else:
        return mat_corr


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
    """Spearman's correlation between every columns in mat_X and mat_Y

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


# fixit: instead call scdrs.pp.compute_stats
# deprecated: use pp.preprocess instead
def compute_gene_contrib(data, gene_list, random_seed=0, copy=False, verbose=False):

    """Find the contribution of each gene to the scDRS score

    Args
    ----
    data (n_cell, n_gene) : AnnData
        adata.X should contain size-normalized log1p transformed count data
    gene_list (n_disease_gene) : list
        Trait gene list
    copy : bool
        If to make copy of the AnnData object

    Returns
    -------
    adata (AnnData):
        Add the columns 'trs_corr'+suffix to data.var
    """

    np.random.seed(random_seed)
    adata = data.copy() if copy else data

    # Pre-compute statistics
    if "var_tech" not in adata.var.columns:
        if verbose:
            print("# score_cell: recompute statistics using method.compute_stats")
        compute_stats(adata)

    # Compute contribution for each gene
    gene_list = sorted(set(gene_list) & set(adata.var_names))
    v_score_weight = 1 / np.sqrt(adata.var.loc[gene_list, "var_tech"].values + 1e-2)
    df_contrib = pd.DataFrame(
        index=adata.obs_names, columns=gene_list, data=adata[:, gene_list].X.toarray()
    )
    df_contrib = df_contrib * v_score_weight

    return df_contrib
