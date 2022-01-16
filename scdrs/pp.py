import numpy as np
import scipy as sp
from scipy import sparse
import pandas as pd
from skmisc.loess import loess


def preprocess(
    data, cov=None, adj_prop=None, n_mean_bin=20, n_var_bin=20, n_chunk=None, copy=False
):
    """
    Preprocess single-cell data for scDRS analysis.

        1. Correct covariates by regressing out the covariates (including
        a constant term) and adding back the original mean for each gene.

        2. Compute gene-level and cell-level statistics for the
        covariate-corrected data.

    Information is stored in `data.uns["SCDRS_PARAM"]`. It operates in
    implicit-covariate-correction mode when `data.X` is sparse and `cov`
    not `None` to improve memory efficiency; it operates in normal mode
    otherwise.

    In normal mode, `data.X` is replaced by the covariate-corrected data.

    In implicit-covariate-correction mode, the covariate correction information
    is stored in `data.uns["SCDRS_PARAM"]` but is not explicitly applied to
    `data.X`, so that `data.X` is always sparse. Subsequent computations on
    the covariate-corrected data are based on the original data `data.X` and
    the covariate correction information. Specifically,

        CORRECTED_X = data.X + COV_MAT * COV_BETA + COV_GENE_MEAN

    The `adj_prop` option is used for adjusting for cell group proportions,
    where each cell is inversely weighted proportional to its corresponding
    cell group size for computing expression mean and variance for genes.
    For stability, the smallest group size is set to be at least 1% of the
    largest group size.


    Parameters
    ----------
    data : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed
        to be size-factor-normalized and log1p-transformed.
    cov : pandas.DataFrame, default=None
        Covariates of shape (n_cell, n_cov). Should contain
        a constant term and have values for at least 75% cells.
    adj_prop : str, default=None
        Cell group annotation (e.g., cell type) used for adjusting for cell group proportions.
        `adj_prop` should be present in `data.obs.columns`.
    n_mean_bin : int, default=20
        Number of mean-expression bins for matching control genes.
    n_var_bin : int, default=20
        Number of expression-variance bins for matching control genes.
    n_chunk : int, default=None
        Number of chunks to split the data into when computing mean and variance
        using _get_mean_var_implicit_cov_corr. If n_chunk is None, set to 5/sparsity.
    copy : bool, default=False
        Return a copy instead of writing to data.


    Returns
    -------
    Overview:
        `data.X` will be updated as the covariate-corrected data in normal mode
        and will stay untouched in the implicit covariate correctoin mode.
        Preprocessing information is stored in `data.uns["SCDRS_PARAM"]`.
    FLAG_SPARSE : bool
        If data.X is sparse.
    FLAG_COV : bool
        If covariate correction is performed.
    COV_MAT : pandas.DataFrame
        Covariate matrix of shape (n_cell, n_cov).
    COV_BETA: pandas.DataFrame
        Covariate effect sizes of shape (n_gene, n_cov).
    COV_GENE_MEAN: pandas.Series
        Gene-level mean expression.
    GENE_STATS : pandas.DataFrame
        Gene-level statistics of shape (n_gene, 7):

        - "mean" : mean expression in log scale.
        - "var" : expression variance in log scale.
        - "var_tech" : technical variance in log scale.
        - "ct_mean" : mean expression in original non-log scale.
        - "ct_var" : expression variance in original non-log scale.
        - "ct_var_tech" : technical variance in original non-log scale.
        - "mean_var" : n_mean_bin * n_var_bin mean-variance bins

    CELL_STATS : pandas.DataFrame
        Cell-level statistics of shape (n_cell, 2):

        - "mean" : mean expression in log scale.
        - "var" : variance expression in log scale.


    Notes
    -----
    Covariate regression:
        adata.X =  cov * beta + resid_X.
    scDRS saves:
        COV_MAT = cov, COV_BETA = (-beta), COV_GENE_MEAN = adata.X.mean(axis=0)
    The scDRS covariate-corrected data:
        CORRECTED_X = resid_X + GENE_MEAN = adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN.


    """

    adata = data.copy() if copy else data
    n_cell, n_gene = adata.shape

    # Parameters and flags
    flag_sparse = sparse.issparse(adata.X)
    flag_cov = cov is not None
    flag_adj_prop = adj_prop is not None
    adata.uns["SCDRS_PARAM"] = {
        "FLAG_SPARSE": flag_sparse,
        "FLAG_COV": flag_cov,
        "FLAG_ADJ_PROP": flag_adj_prop,
    }

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

        # Gene mean: numpy.ndarray of shape (n_gene,)
        v_gene_mean = np.array(adata.X.mean(axis=0)).flatten()
        if flag_sparse:
            # Sparse mode: save correction information
            mat_beta = np.linalg.solve(
                np.dot(df_cov.values.T, df_cov.values) / n_cell,
                sparse.csr_matrix.dot(df_cov.values.T, adata.X) / n_cell,
            )

            adata.uns["SCDRS_PARAM"]["COV_MAT"] = df_cov
            adata.uns["SCDRS_PARAM"]["COV_BETA"] = pd.DataFrame(
                -mat_beta.T, index=adata.var_names, columns=df_cov.columns
            )
            adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"] = pd.Series(
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
        if n_chunk is None:
            n_chunk = 5 * adata.shape[0] * adata.shape[1] // adata.X.data.shape[0] + 1
    else:
        implicit_cov_corr = False
        if n_chunk is None:
            n_chunk = 20

    if flag_adj_prop:
        err_msg = "'adj_prop'=%s not in 'adata.obs.columns'" % adj_prop
        assert adj_prop in adata.obs, err_msg
        err_msg = (
            "On average <10 cells per group, maybe `%s` is not categorical?" % adj_prop
        )
        assert adata.obs[adj_prop].unique().shape[0] < 0.1 * n_cell, err_msg

        temp_df = adata.obs[[adj_prop]].copy()
        temp_df["cell"] = 1
        temp_df = temp_df.groupby(adj_prop).agg({"cell": len})
        temp_df["cell"].clip(lower=int(0.01 * temp_df["cell"].max()), inplace=True)
        temp_dic = {x: n_cell / temp_df.loc[x, "cell"] for x in temp_df.index}

        cell_weight = np.array([temp_dic[x] for x in adata.obs[adj_prop]])
        cell_weight = cell_weight / cell_weight.mean()
    else:
        cell_weight = None

    df_gene, df_cell = compute_stats(
        adata,
        implicit_cov_corr=implicit_cov_corr,
        cell_weight=cell_weight,
        n_mean_bin=n_mean_bin,
        n_var_bin=n_var_bin,
        n_chunk=n_chunk,
    )

    adata.uns["SCDRS_PARAM"]["GENE_STATS"] = df_gene
    adata.uns["SCDRS_PARAM"]["CELL_STATS"] = df_cell

    return adata if copy else None


def compute_stats(
    adata,
    implicit_cov_corr=False,
    cell_weight=None,
    n_mean_bin=20,
    n_var_bin=20,
    n_chunk=20,
):
    """
    Compute gene-level and cell-level statstics used for scDRS analysis. `adata`
    should be log-scale. It has two modes. In the normal mode, it computes
    statistics for `adata.X`. In the implicit covariate correction mode, the
    covariate correction has not been performed on `adata.X` but the corresponding
    information is stored in `adata.uns["SCDRS_PARAM"]`. In this case, it computes
    statistics for the covariate-corrected data

        `transformed_X = adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`

    Parameters
    ----------
    adata : anndata.AnnData
        Single-cell data of shape (n_cell, n_gene). Assumed to be log-scale.
    implicit_cov_corr : bool, default=False
        If True, compute statistics for the implicit corrected data
        `adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`. Otherwise, compute
        for the original data `adata.X`.
    cell_weight : array_like, default=None
        Cell weights of length `adata.shape[0]` for cells in `adata`,
        used for computing weighted gene-level statistics.
    n_mean_bin : int, default=20
        Number of mean-expression bins for matching control genes.
    n_var_bin : int, default=20
        Number of expression-variance bins for matching control genes.
    n_chunk : int, default=20
        Number of chunks to split the data into when computing mean and variance
        using _get_mean_var_implicit_cov_corr.

    Returns
    -------
    df_gene : pandas.DataFrame
        Gene-level statistics of shape (n_gene, 7):

        - "mean" : mean expression in log scale.
        - "var" : variance expression in log scale.
        - "var_tech" : technical variance in log scale.
        - "ct_mean" : mean expression in original non-log scale.
        - "ct_var" : variance expression in original non-log scale.
        - "ct_var_tech" : technical variance in original non-log scale.
        - "mean_var" : n_mean_bin * n_var_bin mean-variance bins
    df_cell : pandas.DataFrame
        Cell-level statistics of shape (n_cell, 2):

        - "mean" : mean expression in log scale.
        - "var" : variance expression in log scale.
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
            "var_tech",
            "ct_mean",
            "ct_var",
            "ct_var_tech",
            "mean_var",
        ],
    )
    df_cell = pd.DataFrame(index=adata.obs_names, columns=["mean", "var"])
    # Gene-level statistics
    if not implicit_cov_corr:
        # Normal mode
        df_gene["mean"], df_gene["var"] = _get_mean_var(
            adata.X, axis=0, weights=cell_weight
        )
        # Get the mean and var for the non-log-scale size-factor-normalized counts
        # It is highly correlated to the non-size-factor-normalized counts
        if sparse.issparse(adata.X):  # sp sparse matrix
            temp_X = adata.X.copy().expm1()  # exp(X)-1 to get ct matrix from logct
        else:
            temp_X = np.expm1(adata.X)  # numpy ndarray
        df_gene["ct_mean"], df_gene["ct_var"] = _get_mean_var(
            temp_X, axis=0, weights=cell_weight
        )
        del temp_X
    else:
        # Implicit covariate correction mode
        df_gene["mean"], df_gene["var"] = _get_mean_var_implicit_cov_corr(
            adata, axis=0, n_chunk=n_chunk, weights=cell_weight
        )
        df_gene["ct_mean"], df_gene["ct_var"] = _get_mean_var_implicit_cov_corr(
            adata, transform_func=np.expm1, axis=0, n_chunk=n_chunk, weights=cell_weight
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
    for bin_ in sorted(set(v_mean_bin)):
        ind_select = v_mean_bin == bin_
        v_var_bin = pd.qcut(
            df_gene.loc[ind_select, "var"], n_var_bin, labels=False, duplicates="drop"
        )
        df_gene.loc[ind_select, "mean_var"] = ["%d.%d" % (bin_, x) for x in v_var_bin]

    # Cell-level statistics
    if not implicit_cov_corr:
        # Normal mode
        df_cell["mean"], df_cell["var"] = _get_mean_var(adata.X, axis=1)
    else:
        # Implicit covariate correction mode
        df_cell["mean"], df_cell["var"] = _get_mean_var_implicit_cov_corr(
            adata, axis=1, n_chunk=n_chunk
        )

    return df_gene, df_cell


##############################################################################
######################## Preprocessing Subroutines ###########################
##############################################################################
def reg_out(mat_Y, mat_X):
    """Regress mat_X out of mat_Y.

    Parameters
    ----------
    mat_Y : np.ndarray
        Response variable of shape (n_sample, n_response).
    mat_X : np.ndarray
        Covariates of shape (n_sample, n_covariates).

    Returns
    -------
    mat_Y_resid : np.ndarray
        Response variable residual of shape (n_sample, n_response).
    """

    if sparse.issparse(mat_X):
        mat_X = mat_X.toarray()
    else:
        mat_X = np.array(mat_X)
    if len(mat_X.shape) == 1:
        mat_X = mat_X.reshape([-1, 1])

    if sparse.issparse(mat_Y):
        mat_Y = mat_Y.toarray()
    else:
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


def _get_mean_var(sparse_X, axis=0, weights=None):
    """
    Compute mean and var of a sparse / non-sparse matrix.

    Parameters
    ----------
    sparse_X : array_like
        Data matrix (can be dense/sparse).
    axis : {0, 1}, default=0
        Axis along which to compute mean and variance.
    weights : array_like, default=None
        Weights of length `sparse_X.shape[axis]`.

    Returns
    -------
    v_mean : np.ndarray
        Mean vector.
    v_var : np.ndarray
        Variance vector.

    """

    if weights is None:
        if sparse.issparse(sparse_X):
            # Sparse + unweighted
            v_mean = sparse_X.mean(axis=axis)
            v_mean = np.array(v_mean).reshape([-1])
            v_var = sparse_X.power(2).mean(axis=axis)
            v_var = np.array(v_var).reshape([-1])
            v_var = v_var - v_mean ** 2
        else:
            # Dense + unweighted
            sparse_X = np.array(sparse_X)
            v_mean = np.mean(sparse_X, axis=axis)
            v_var = np.var(sparse_X, axis=axis)

    else:
        weights = np.array(weights)
        if sparse.issparse(sparse_X):
            # Sparse + weighted
            v_mean = _weighted_sparse_average(sparse_X, weights, axis=axis)
            v_var = _weighted_sparse_average(sparse_X.power(2), weights, axis=axis)
            v_var = v_var - v_mean ** 2
        else:
            # Dense + weighted
            sparse_X = np.array(sparse_X)
            v_mean = np.average(sparse_X, axis=axis, weights=weights)
            v_var = np.average(np.square(sparse_X), axis=axis, weights=weights)
            v_var = v_var - v_mean ** 2

    return v_mean, v_var


def _weighted_sparse_average(sparse_X, weights, axis=0):
    """
    Compute weighted mean for a sparse matrix.

    Parameters
    ----------
    sparse_X : array_like
        Sparse data matrix.
    weights : array_like
        Weights of length `sparse_X.shape[axis]`.
    axis : {0, 1}, default=0
        Axis along which to compute mean.

    Returns
    -------
    v_mean : np.ndarray
        Mean vector.
    """

    if axis == 0:
        v_mean = sparse_X.T.multiply(weights).mean(axis=1)
        v_mean = np.array(v_mean).reshape([-1]) / np.mean(weights)
        return v_mean

    if axis == 1:
        v_mean = sparse_X.multiply(weights).mean(axis=1)
        v_mean = np.array(v_mean).reshape([-1]) / np.mean(weights)
        return v_mean


def _get_mean_var_implicit_cov_corr(
    adata, axis=0, weights=None, transform_func=None, n_chunk=20
):
    """
    Compute mean and variance of sparse matrix of the form

        `adata.X + COV_MAT * COV_BETA + COV_GENE_MEAN`.

    Computed iteratively over chunks of sparse matrix by converting to
    dense matrix and computing mean and variance of dense matrix.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix (n_obs, n_vars)
    axis : {0, 1}, default=0
        Axis along which to compute mean and variance
    weights : array_like, default=None
        Weights of length `adata.shape[axis]`.
    transform_func : function, default=None
        Function to transform the data before computing mean and variance
    n_chunk : int, default=20
        Number of chunks to split the data into when computing mean and variance
        this will determine the memory usage

    Returns
    -------
    v_mean : np.ndarray
        Mean vector.
    v_var : np.ndarray
        Variance vector.
    """

    assert axis in [0, 1], "axis must be one of [0, 1]"
    assert (
        "SCDRS_PARAM" in adata.uns
    ), "adata.uns['SCDRS_PARAM'] not found, run `preprocess` before calling this function"

    n_obs, n_gene = adata.shape

    # COV INFO: cov_mat: (n_cell, n_cov) cov_beta (n_cov, n_gene)
    cell_list = list(adata.obs_names)
    cov_list = list(adata.uns["SCDRS_PARAM"]["COV_MAT"])
    gene_list = list(adata.var_names)
    cov_mat = adata.uns["SCDRS_PARAM"]["COV_MAT"].loc[cell_list, cov_list].values
    cov_beta = adata.uns["SCDRS_PARAM"]["COV_BETA"].loc[gene_list, cov_list].values.T
    gene_mean = adata.uns["SCDRS_PARAM"]["COV_GENE_MEAN"].loc[gene_list].values

    if transform_func is None:
        transform_func = lambda x: x

    if weights is not None:
        weights = np.array(weights)

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
                + gene_mean[start:stop]
            )
            chunk_X = transform_func(chunk_X)
            v_mean[start:stop], v_var[start:stop] = _get_mean_var(
                chunk_X, axis=axis, weights=weights
            )
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
            v_mean[start:stop], v_var[start:stop] = _get_mean_var(
                chunk_X, axis=axis, weights=weights
            )
            start = stop

    return v_mean, v_var
