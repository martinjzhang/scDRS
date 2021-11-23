import scanpy as sc
import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
from skmisc.loess import loess
from tqdm import tqdm
import scdrs


def preprocess(data, cov=None, n_mean_bin=20, n_var_bin=20, copy=False):
    """
    Preprocess single-cell data for scDRS analysis:
        1. Regress out covariates and add back the original mean for each gene.
        2. Compute gene-level and cell-level statistics.

    When data.X is sparse and cov is present, scDRS operates in the implicity
    covariate correction mode. In this mode, data.X is always sparse.
    Covariate correction information is stored in data.uns["SCDRS_PARAM"] but
    is not explicity applied to data.X.

    In other cases, scDRS operates in the normal mode, where covariate correction
    is applied to data.X.


    The preproecssing information is stored in data.uns["SCDRS_PARAM"].

    Parameters
    ----------
    data : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    cov : pandas.DataFrame, default=None
        Covariates of shape (n_cell, n_cov). Should contain
        a constant term and have values for at least 75% cells.
    n_mean_bin : int, default=20
        Number of mean-expression bins for matching control genes.
    n_var_bin : int, default=20
        Number of expression-variance bins for matching control genes.
    copy : bool, default=False
        Return a copy instead of writing to data.


    Returns
    -------
    FLAG_SPARSE : bool
        If data.X is sparse.  Stored as a key in data.uns["SCDRS_PARAM"].
    FLAG_COV : bool
        If covariate correction is performed. Stored as a key in
        data.uns["SCDRS_PARAM"].
    COV_MAT : pandas.DataFrame
        Covariance matrix of shape (n_cell, n_cov). Stored as a key
        in data.uns["SCDRS_PARAM"].
    COV_BETA: pandas.DataFrame
        Covariance effect sizes of shape (n_gene, n_cov). Stored as a key
        in data.uns["SCDRS_PARAM"].
    COV_GENE_MEAN: pandas.Series
        Mean
    GENE_STATS : pandas.DataFrame
        Gene-level statistics of shape (n_gene, n_stats). Stored as
        a key in data.uns["SCDRS_PARAM"].
    CELL_STATS : pandas.DataFrame
        Cell-level statistics of shape (n_cell, n_stats). Stored as
        a key in data.uns["SCDRS_PARAM"].


    Notes
    -----
    Covariate regression:
        adata.X =  cov * beta + resid_X.
    scDRS saves:
        COV_MAT = cov, COV_BETA = (-beta), COV_GENE_MEAN = adata.X.mean(axis=0)
    The scDRS covariate-corrected data:
        transformed_X = resid_X + GENE_MEAN = adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN.


    """

    adata = data.copy() if copy else data
    n_cell, n_gene = adata.shape

    # Parameters and flags
    flag_sparse = sp.sparse.issparse(adata.X)
    flag_cov = cov is not None
    adata.uns["SCDRS_PARAM"] = {"FLAG_SPARSE": flag_sparse, "FLAG_COV": flag_cov}

    # Update adata.X
    if flag_sparse:
        # Force sparse.csr_matrix for the sparse mode
        if not isinstance(adata.X, sparse.csr_matrix):
            adata.X = sparse.csr_matrix(adata.X)
    else:
        # Force np.ndarray for the dense mode
        if not isinstance(adata.X, np.ndarray):
            adata.X = np.array(adata.X)

    # Covariate correction
    if flag_cov:
        # Check if cells in data and cov are consistent
        assert (
            len(set(cov.index) & set(adata.obs_names)) > 0.75 * n_cell
        ), "cov does not match the cells in data"

        df_cov = pd.DataFrame(index=adata.obs_names)
        df_cov = df_cov.join(cov)
        df_cov.fillna(df_cov.mean(), inplace=True)

        # Add const term if df_cov does not already have it (or a linear combination of it)
        v_resid = reg_out(np.ones(n_cell), df_cov.values)
        if (v_resid ** 2).mean() > 0.01:
            df_cov["SCDRS_CONST"] = 1

        v_gene_mean = np.array(adata.X.mean(axis=0)).flatten() # numpy.ndarray of shape (n_gene,)
        if flag_sparse:
            # Sparse mode: save correction information
            mat_beta = np.linalg.solve(
                np.dot(df_cov.values.T, df_cov.values) / n_cell,
                sp.sparse.csr_matrix.dot(df_cov.values.T, adata.X) / n_cell,
            )

            adata.uns["SCDRS_PARAM"]["COV_MAT"] = df_cov
            adata.uns["SCDRS_PARAM"]["COV_BETA"] = pd.DataFrame(
                -mat_beta.T, index=adata.var_names, columns=df_cov.columns
            )
            adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"] = pd.DataFrame(
                v_gene_mean, index=adata.var_names
            )
        else:
            # Dense mode: regress out covariate and add back mean
            adata.X = reg_out(adata.X, df_cov.values)
            adata.X += v_gene_mean
            
#             # Note: this version (dense+cov) should produce the exact toydata results 
#             adata.var["mean"] = adata.X.mean(axis=0).T
#             adata.X -= adata.var["mean"].values
#             adata.X = reg_out(adata.X, df_cov[['SCDRS_CONST', 'covariate']].values)
#             adata.X += adata.var["mean"]


    # Precompute for each gene and mean&var for each cell
    if flag_sparse and flag_cov:
        implicit_cov_corr = True
    else:
        implicit_cov_corr = False

    df_gene, df_cell = compute_stats(
        adata,
        n_mean_bin=n_mean_bin,
        n_var_bin=n_var_bin,
        implicit_cov_corr=implicit_cov_corr,
    )

    adata.uns["SCDRS_PARAM"]["GENE_STATS"] = df_gene
    adata.uns["SCDRS_PARAM"]["CELL_STATS"] = df_cell

    return adata if copy else None


def compute_stats(adata, implicit_cov_corr=False, n_mean_bin=20, n_var_bin=20):
    """
    Compute gene-level and cell-level statstics used for scDRS analysis. `adata`
    should be log-scale. It has two modes. In the normal mode, it computes
    statistics for `adata.X`. In the implicit covariate correction mode, the
    covariate correction has not been performed on `adata.X` but the corresponding
    information is stored in `adata.uns["SCDRS_PARAM"]`. In this case, it computes
    statistics for the covariate-corrected data
    `transformed_X = adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`.

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed to be log-scale.
    implicit_cov_corr : bool, default=False
        If True, compute statistics for the implicit corrected data
        `adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`. Otherwise, compute
        for the original data `adata.X`.
    n_mean_bin : int, default=20
        Number of mean-expression bins for matching control genes.
    n_var_bin : int, default=20
        Number of expression-variance bins for matching control genes.

    Returns
    -------
    df_gene : pandas.DataFrame
        Gene-level statistics of shape (n_gene, 7). Contains
        ["mean", "var", "ct_mean", "ct_var", "ct_var_tech", "var_tech", "mean_var"].
    df_cell : pandas.DataFrame
        Cell-level statistics of shape (n_cell, 2). Contains ["mean", "var"].
    """

    if implicit_cov_corr:
        assert (
            "SCDRS_PARAM" in adata.uns
        ), "adata.uns['SCDRS_PARAM'] is not found, run `scdrs.pp.preprocess` before calling this function"

    df_gene = pd.DataFrame(
        index=adata.var_names,
        columns=[
            "mean",
            "var",
            "ct_mean",
            "ct_var",
            "ct_var_tech",
            "var_tech",
            "mean_var",
        ],
    )
    df_cell = pd.DataFrame(index=adata.obs_names, columns=["mean", "var"])
    # Gene-level statistics
    if not implicit_cov_corr:
        # Normal mode
        df_gene["mean"], df_gene["var"] = _get_mean_var(adata.X, axis=0)
        # Get the mean and var for the non-log-scale size-factor-normalized counts
        # It is highly correlated to the non-size-factor-normalized counts
        if sp.sparse.issparse(adata.X):  # sp sparse matrix
            temp_X = adata.X.copy().expm1()  # exp(X)-1 to get ct matrix from logct
        else:
            temp_X = np.expm1(adata.X)  # numpy ndarray
        df_gene["ct_mean"], df_gene["ct_var"] = _get_mean_var(temp_X, axis=0)
        del temp_X
    else:
        # Implicit covariate correction mode
        df_gene["mean"], df_gene["var"] = _get_mean_var_implicit_cov_corr(adata, axis=0)
        df_gene["ct_mean"], df_gene["ct_var"] = _get_mean_var_implicit_cov_corr(
            adata, transform_func=np.expm1, axis=0
        )

    # Borrowed from scanpy _highly_variable_genes_seurat_v3
    not_const = df_gene["ct_var"].values > 0
    estimat_var = np.zeros(adata.shape[1], dtype=np.float64)
    y = np.log10(df_gene["ct_var"].values[not_const])
    x = np.log10(df_gene["ct_mean"].values[not_const])
    model = loess(x, y, span=0.3, degree=2)
    model.fit()
    estimat_var[not_const] = model.outputs.fitted_values
    df_gene["ct_var_tech"] = 10 ** estimat_var
    # Recipe from Frost Nucleic Acids Research 2020
    df_gene["var_tech"] = df_gene["var"] * df_gene["ct_var_tech"] / df_gene["ct_var"]
    df_gene.loc[df_gene["var_tech"].isna(), "var_tech"] = 0

    # Add n_mean_bin*n_var_bin mean_var bins
    n_bin_max = np.floor(np.sqrt(adata.shape[1] / 10)).astype(int)
    if (n_mean_bin > n_bin_max) | (n_var_bin > n_bin_max):
        n_mean_bin, n_var_bin = n_bin_max, n_bin_max
        print(
            "Too few genes for 20*20 bins, setting n_mean_bin=n_var_bin=%d"
            % (n_bin_max)
        )
    v_mean_bin = pd.qcut(df_gene["mean"], n_mean_bin, labels=False, duplicates="drop")
    df_gene["mean_var"] = ""
    for bin_ in set(v_mean_bin):
        ind_select = v_mean_bin == bin_
        v_var_bin = pd.qcut(
            df_gene.loc[ind_select, "var"], n_var_bin, labels=False, duplicates="drop"
        )
        df_gene.loc[ind_select, "mean_var"] = [
            "%s.%s" % (x, y) for x, y in zip(v_mean_bin[ind_select], v_var_bin)
        ]

    # Cell-level statistics
    if not implicit_cov_corr:
        # Normal mode
        df_cell["mean"], df_cell["var"] = _get_mean_var(adata.X, axis=1)
    else:
        # Implicit covariate correction mode
        df_cell["mean"], df_cell["var"] = _get_mean_var_implicit_cov_corr(adata, axis=1)

    return df_gene, df_cell


##############################################################################
############################### Subroutines ##################################
##############################################################################
def reg_out(mat_Y, mat_X):
    """Regress mat_X out of mat_Y

    Args
    ----
    mat_Y (n_sample, n_response) : np.ndarray
        Response variable
    mat_X (n_sample, n_covariates) : np.ndarray
        Covariates

    Returns
    -------
    mat_Y_resid (n_sample, n_response) : np.ndarray
        Response variable residual
    """

    mat_X = np.array(mat_X)
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])
    mat_Y = np.array(mat_Y)
    if len(mat_Y.shape) == 1:
        mat_Y = mat_Y.reshape([-1, 1])

    n_sample = mat_Y.shape[0]
    mat_xtx = np.dot(mat_X.T, mat_X) / n_sample
    mat_xty = np.dot(mat_X.T, mat_Y) / n_sample
    mat_coef = np.linalg.solve(mat_xtx, mat_xty)
    mat_Y_resid = mat_Y - mat_X.dot(mat_coef)

    if mat_Y_resid.shape[1] == 1:
        mat_Y_resid = mat_Y_resid.reshape([-1])

    return mat_Y_resid


def _get_mean_var(sparse_X, axis=0):
    """
    Compute mean and var of a sparse / non-sparse matrix.
    """

    if sp.sparse.issparse(sparse_X):
        v_mean = sparse_X.mean(axis=axis)
        v_mean = np.array(v_mean).reshape([-1])
        v_var = sparse_X.power(2).mean(axis=axis)
        v_var = np.array(v_var).reshape([-1])
        v_var = v_var - v_mean ** 2
    else:
        v_mean = np.mean(sparse_X, axis=axis)
        v_var = np.var(sparse_X, axis=axis)

    return v_mean, v_var


def _get_mean_var_implicit_cov_corr(adata, axis=0, transform_func=None, n_chunk=20):
    """
    Compute mean and variance of sparse matrix of the form
    `adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`. Computed
    iteratively over chunks of sparse matrix by converting to
    dense matrix and computing mean and variance of dense matrix.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix (n_obs, n_vars)
    axis : {0, 1}, default=0
        Axis along which to compute mean and variance
    transform_func : function, default=None
        Function to transform the data before computing mean and variance
    n_chunk : int, default=20
        Number of chunks to split the data into when computing mean and variance
        this will determine the memory usage
    """

    assert axis in [0, 1], "axis must be one of [0, 1]"
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `preprocess` before calling this function"

    n_obs, n_gene = adata.shape
    # cov_mat: (n, d) cov_beta (d, p)
    cov_mat = adata.uns["SCDRS_PARAM"]["COV_MAT"].values
    cov_beta = adata.uns["SCDRS_PARAM"]["COV_BETA"].values.T
    gene_mean = adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].values.T

    if transform_func is None:
        transform_func = lambda x: x

    # mean / var of transform_func(X + cov_mat @ cov_beta + mean)
    if axis == 0:
        # output would have the shape of (n_gene, )
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
        # output would have the shape of (n_obs, )
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
