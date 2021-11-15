import numpy as np
import pandas as pd
from skmisc.loess import loess
import anndata


def sparse_compute_stats(
    data: anndata.AnnData, n_mean_bin: int = 20, n_var_bin: int = 20, copy: bool = False
):
    """
    Compute the mean, variance, and technical variance of the log1p transformed count data

    Parameters
    ----------
    data : anndata.AnnData
        Annotated data matrix (n_obs, n_vars)
    n_mean_bin : int
        Number of bins for the mean, default 20
    n_var_bin : int
        Number of bins for the variance, default 20
    copy : bool
        Whether to operate a copy of the data or inplace, default False

    Returns
    -------
    data : anndata.AnnData
        Annotated data matrix with mean, var calculated.
    """

    # Gene-wise statistics
    adata = data.copy() if copy else data
    adata.var["mean"], adata.var["var"] = sparse_get_mean_var(adata, axis=0)

    # Get the mean and var for the size-factor-normalized counts
    # It is highly correlated to the non-size-factor-normalized counts
    # exp(X)-1 to get ct matrix from logct

    adata.var["ct_mean"], adata.var["ct_var"] = sparse_get_mean_var(
        adata, transform_func=np.expm1, axis=0
    )

    # Borrowed from scanpy _highly_variable_genes_seurat_v3
    not_const = adata.var["ct_var"].values > 0
    estimat_var = np.zeros(adata.shape[1], dtype=np.float64)
    y = np.log10(adata.var["ct_var"].values[not_const])
    x = np.log10(adata.var["ct_mean"].values[not_const])
    model = loess(x, y, span=0.3, degree=2)
    model.fit()
    estimat_var[not_const] = model.outputs.fitted_values
    adata.var["ct_var_tech"] = 10 ** estimat_var
    # Recipe from Frost Nucleic Acids Research 2020
    adata.var["var_tech"] = (
        adata.var["var"] * adata.var["ct_var_tech"] / adata.var["ct_var"]
    )
    adata.var.loc[adata.var["var_tech"].isna(), "var_tech"] = 0

    # Add n_mean_bin*n_var_bin mean_var bins
    n_bin_max = np.floor(np.sqrt(adata.shape[1] / 10)).astype(int)
    if (n_mean_bin > n_bin_max) | (n_var_bin > n_bin_max):
        n_mean_bin, n_var_bin = n_bin_max, n_bin_max
        print(
            "Too few genes for 20*20 bins, setting n_mean_bin=n_var_bin=%d"
            % (n_bin_max)
        )
    v_mean_bin = pd.qcut(adata.var["mean"], n_mean_bin, labels=False, duplicates="drop")
    adata.var["mean_var"] = ""
    for bin_ in set(v_mean_bin):
        ind_select = v_mean_bin == bin_
        v_var_bin = pd.qcut(
            adata.var.loc[ind_select, "var"], n_var_bin, labels=False, duplicates="drop"
        )
        adata.var.loc[ind_select, "mean_var"] = [
            "%s.%s" % (x, y) for x, y in zip(v_mean_bin[ind_select], v_var_bin)
        ]

    # Cell-wise statistics
    adata.obs["mean"], adata.obs["var"] = sparse_get_mean_var(adata, axis=1)

    return adata if copy else None


def sparse_compute_raw_score(adata, gene_list, gene_weight, weight_opt):
    """Compute raw score
        v_score_weight = gene_weight * {uniform/vs/inv_std}
        `_SCDRS_SPARSE` is assumed to have been computed using `sparse_reg_out`

    Parameters
    ----------
    adata (n_cell, n_gene) : AnnData
        adata.X should contain raw count data
    gene_list (n_trait_gene) : list
        Trait gene list
    gene_weight (n_trait_gene) : list/np.ndarray
        Gene weights for genes in the gene_list
    weight_opt : str
        Option for computing the raw score
        - 'uniform': average over the genes in the gene_list
        - 'vs': weighted average with weights equal to 1/sqrt(technical_variance_of_logct)
        - 'inv_std': weighted average with weights equal to 1/std
        - 'od': overdispersion score

    Returns
    -------
    v_raw_score (n_cell,) : np.ndarray
        Raw score
    v_score_weight (n_trait_gene,) : np.ndarray
        Gene weights score
    """
    assert "_SCDRS_SPARSE" in adata.uns

    # Compute raw score (weighted average)
    if weight_opt == "uniform":
        v_score_weight = np.ones(len(gene_list))
    if weight_opt == "vs":
        v_score_weight = 1 / np.sqrt(adata.var.loc[gene_list, "var_tech"].values + 1e-2)
    if weight_opt == "inv_std":
        v_score_weight = 1 / np.sqrt(adata.var.loc[gene_list, "var"].values + 1e-2)

    if gene_weight is not None:
        v_score_weight = v_score_weight * np.array(gene_weight)
    v_score_weight = v_score_weight / v_score_weight.sum()

    cov_mat = adata.uns["_SCDRS_SPARSE"]["COV_MAT"]
    cov_beta = adata.uns["_SCDRS_SPARSE"]["COV_BETA"]
    gene_mean = adata.uns["_SCDRS_SPARSE"]["GENE_MEAN"]

    # calclulate v_raw_score = transformed_X @ v_score_weight
    # where: transformed_X = adata.X + cov_mat @ cov_beta + gene_mean
    v_raw_score = (
        adata[:, gene_list].X.dot(v_score_weight)
        + cov_mat @ (cov_beta.loc[gene_list, :].values.T @ v_score_weight)
        + gene_mean.loc[gene_list, :].values.T @ v_score_weight
    ).flatten()

    return v_raw_score, v_score_weight


def sparse_reg_out(adata, cov):
    """
    Compute the sparse regression output

    Parameters
    ----------
    adata (n_cell, n_gene) : AnnData
        adata.X should contain size-normalized log1p transformed count data
    cov (n_cell, n_gene) : np.ndarray
        Covariance matrix

    Returns
    -------
    adata.uns["_SCDRS_SPARSE"] : dict
        Dictionary containing the following keys:
        - "COV_MAT": Covariance matrix
        - "COV_BETA": Covariance beta
        - "GENE_MEAN": Gene mean
    """
    from scipy import sparse

    assert np.all(cov[:, 0] == 1), "first column of covariate matrix must be 1"
    n_obs, n_gene = adata.X.shape
    assert n_obs == cov.shape[0]
    n_cov = cov.shape[1]

    # calculate the gene mean
    mean = adata.X.mean(axis=0)

    # regress out the covariates
    # beta = (cov.T * cov)^-1 (cov.T * X)
    # X = (X - cov * beta) + m = X + cov * beta + m
    b = np.linalg.solve(
        np.dot(cov.T, cov) / n_obs, sparse.csr_matrix.dot(cov.T, adata.X) / n_obs
    )
    # transformed_X = X + cov * beta + m
    # or transformed_X = adata.X + COV_MAT * COV_BETA + GENE_MEAN
    # because `cov` is assumed to contain intercept, the centering of `adata.X` is
    # done by adding `COV_MAT * COV_BETA`

    adata.uns["_SCDRS_SPARSE"] = {
        "COV_BETA": pd.DataFrame(-b.T, index=adata.var_names),
        "COV_MAT": cov,
        "GENE_MEAN": pd.DataFrame(mean.A1, index=adata.var_names),
    }


def sparse_get_mean_var(
    adata: anndata.AnnData, axis: int, transform_func=None, n_chunk: int = 20
):
    """
    Compute mean and variance of sparse matrix iteratively over chunks of sparse matrix
    by converting to dense matrix and computing mean and variance of dense matrix

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix (n_obs, n_vars)
    axis : int
        Axis along which to compute mean and variance
    transform_func : function
        Function to transform the data before computing mean and variance
    n_chunk : int
        Number of chunks to split the data into when computing mean and variance
        this will determine the memory usage
    """

    assert axis in [0, 1], "axis must be one of [0, 1]"
    assert (
        "_SCDRS_SPARSE" in adata.uns
    ), "adata.uns['_SCDRS_SPARSE'] is not found, run `sparse_reg_out` before calling this function"

    n_obs, n_gene = adata.shape
    # cov_mat: (n, d) cov_beta (d, p)
    cov_mat = adata.uns["_SCDRS_SPARSE"]["COV_MAT"]
    cov_beta = adata.uns["_SCDRS_SPARSE"]["COV_BETA"].values.T
    gene_mean = adata.uns["_SCDRS_SPARSE"]["GENE_MEAN"].values.T

    if transform_func is None:
        transform_func = lambda x: x
    # mean / var of (X + cov_mat @ cov_beta + mean)

    if axis == 0:
        # output would be in the shape of (n_gene, )
        v_mean = np.zeros(n_gene)
        v_var = np.zeros(n_gene)
        start = 0
        chunk_size = n_gene // n_chunk
        while start < n_gene:
            stop = min(start + chunk_size, n_gene)
            chunk_X = (
                adata.X[:, start:stop]
                + cov_mat @ cov_beta[:, start:stop]
                + gene_mean[:, start:stop]
            )
            chunk_X = transform_func(chunk_X)
            v_mean[start:stop] = np.mean(chunk_X, axis).A1
            v_var[start:stop] = np.var(chunk_X, axis).A1
            start = stop

    elif axis == 1:
        v_mean = np.zeros(n_obs)
        v_var = np.zeros(n_obs)
        start = 0
        chunk_size = n_obs // n_chunk
        while start < n_obs:
            stop = min(start + chunk_size, n_obs)
            chunk_X = (
                adata.X[start:stop, :] + cov_mat[start:stop] @ cov_beta + gene_mean
            )
            chunk_X = transform_func(chunk_X)

            v_mean[start:stop] = np.mean(chunk_X, axis=axis).A1
            v_var[start:stop] = np.var(chunk_X, axis=axis).A1
            start = stop

    return v_mean, v_var
