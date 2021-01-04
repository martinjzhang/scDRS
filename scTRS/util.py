import numpy as np
import scipy as sp
import pandas as pd
import numbers
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests

def test_sctrs():
    print('# test_sctrs')
    return 

def get_sparse_var(sparse_X, axis=0):
    
    v_mean = sparse_X.mean(axis=axis)
    v_mean = np.array(v_mean).reshape([-1])
    v_var = sparse_X.power(2).mean(axis=axis)
    v_var = np.array(v_var).reshape([-1])
    v_var = v_var - v_mean**2
    
    return v_mean,v_var

# https://stats.stackexchange.com/questions/403652/two-sample-quantile-quantile-plot-in-python
def qqplot(x, y, quantiles=None, interpolation='nearest', ax=None, **kwargs):
    """Draw a quantile-quantile plot for `x` versus `y`.

    Parameters
    ----------
    x, y : array-like
        One-dimensional numeric arrays.

    ax : matplotlib.axes.Axes, optional
        Axes on which to plot. If not provided, the current axes will be used.

    quantiles : int or array-like, optional
        Quantiles to include in the plot. This can be an array of quantiles, in
        which case only the specified quantiles of `x` and `y` will be plotted.
        If this is an int `n`, then the quantiles will be `n` evenly spaced
        points between 0 and 1. If this is None, then `min(len(x), len(y))`
        evenly spaced quantiles between 0 and 1 will be computed.

    interpolation : {‘linear’, ‘lower’, ‘higher’, ‘midpoint’, ‘nearest’}
        Specify the interpolation method used to find quantiles when `quantiles`
        is an int or None. See the documentation for numpy.quantile().

    kwargs : dict of keyword arguments
        Keyword arguments to pass to matplotlib.axes.Axes.scatter() when drawing
        the q-q plot.
    """
    # Get current axes if none are provided
    if ax is None:
        ax = plt.gca()

    if quantiles is None:
        quantiles = min(len(x), len(y))

    # Compute quantiles of the two samples
    if isinstance(quantiles, numbers.Integral):
        quantiles = np.linspace(start=0, stop=1, num=int(quantiles))
    else:
        quantiles = np.atleast_1d(np.sort(quantiles))
    x_quantiles = np.quantile(x, quantiles, interpolation=interpolation)
    y_quantiles = np.quantile(y, quantiles, interpolation=interpolation)

    # Draw the q-q plot
    ax.scatter(x_quantiles, y_quantiles, **kwargs)
    
def empirical_zsc(score, null):
    """
        Calculate the empirical p-value for two arrays. 
        For each element in `score`, count how many elements in the null are below the score, 
            and derive a p-value.
        This is for a somewhat better implementation.
        
        score: Array-like, e.g., the TRS score
        null: Array-like, e.g., null distribution
        
        Returns empirical p-value, same-length as null
    """
    df = pd.DataFrame({
        'id': np.concatenate([
            [f'score_{i}' for i in np.arange(len(score))],
            [f'null_{i}' for i in np.arange(len(null))]]),
        'val': np.concatenate([score, null]),
        'is_null': np.concatenate([np.zeros(len(score)), np.ones(len(null))])
    })
    df = df.sort_values('val')
    # TODO: does these pseudo-count makes sense?
    df['cum_prop'] = (np.cumsum(df['is_null']) + 1) / (len(null) + 2)
    df = df[df['id'].str.startswith('score_')].sort_index()
    return sp.stats.norm.ppf(df['cum_prop'].values)


# https://github.com/scikit-learn/scikit-learn/blob/0fb307bf3/sklearn/preprocessing/_data.py#L1092
def sparse_robust_scale(X, quantile_range=(25.0, 75.0), copy=True):
    """
        sparse_robust_scale:
        Compute the robust scale for matri X
        
        X: sparse matrix, rows are individuals, columns are genes,
        
        Return: Compute a scale for each gene
    """
    assert sp.sparse.issparse(X), "X must be a sparse matrix"
    # at fit, convert sparse matrices to csc for optimized computation of the quantiles
    X = X.tocsc(copy=copy)
    q_min, q_max = quantile_range
    if not 0 <= q_min <= q_max <= 100:
        raise ValueError("Invalid quantile range: %s" %
                         str(quantile_range))
    
    quantiles = []
    for feature_idx in range(X.shape[1]):
        column_nnz_data = X.data[X.indptr[feature_idx]:
                                 X.indptr[feature_idx + 1]]
        column_data = np.zeros(shape=X.shape[0], dtype=X.dtype)
        column_data[:len(column_nnz_data)] = column_nnz_data

        quantiles.append(np.nanpercentile(column_data, quantile_range))

    quantiles = np.transpose(quantiles)

    return quantiles[1] - quantiles[0]


def gearys_c(adata, vals):
    """Compute Geary's C statistics for an AnnData
    Adopted from https://github.com/ivirshup/scanpy/blob/metrics/scanpy/metrics/_gearys_c.py
    
    C =
        \frac{
            (N - 1)\sum_{i,j} w_{i,j} (x_i - x_j)^2
        }{
            2W \sum_i (x_i - \bar{x})^2
        }
    
    Args:
        adata (AnnData): AnnData object
            adata.obsp["Connectivities] should contain the connectivity graph, 
            with shape `(n_obs, n_obs)`
        vals (Array-like):
            Values to calculate Geary's C for. If one dimensional, should have
            shape `(n_obs,)`. 
    Returns:
        C: the Geary's C statistics
    """
    graph = adata.obsp["connectivities"]
    assert graph.shape[0] == graph.shape[1]
    graph_data = graph.data.astype(np.float_, copy=False)
    assert graph.shape[0] == vals.shape[0]
    assert(np.ndim(vals) == 1)
    
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

def calculate_trs_stats(zsc_dict, zsc_index, stats_dict, adata, stratify_by):
    """
    Calculate statistics of TRS stratified by covariate, e.g., celltype, 
    this function only cope with simple statistics such as mean, sd that
    takes an array of unordered TRS and outputs a single scalar.
    # Arguments
    zsc_dict: z-score dict, trait -> z-score array, sorted according to zsc_index
    zsc_index: index for the z-score in `zsc_dict`
    stats_dict: statistics to be computed, name -> function,
        the function takes an array and output a scalar
    adata: AnnData for the meta information
    stratify_by: which column in `meta_df` to stratify
    # TODO
    1. Add permutation procedure
    """
    assert stratify_by in ["tissue_celltype", "cell_ontology_class"]
    
    meta_df = adata.obs.copy()
    trait_list = list(zsc_dict.keys())
    zsc_df = meta_df.join(pd.DataFrame(zsc_dict, index=zsc_index), how="inner")
    
    stratify_list = sorted(zsc_df[stratify_by].unique())
    rls_df_dict = {stats: pd.DataFrame(index=stratify_list, columns=trait_list, data=0) for stats in stats_dict}
    
    for strata, strata_df in zsc_df.groupby(stratify_by):
        strata_n_cell = strata_df.shape[0]
        for stats_name in stats_dict:
            stats_fun = stats_dict[stats_name]
            if stats_fun is not None:
                for trait in trait_list:
                    rls_df_dict[stats_name].loc[strata, trait] = stats_fun(strata_df[trait].values)
            else:
                if stats_name == "gearysc":
                    for trait in trait_list:
                        rls_df_dict[stats_name].loc[strata, trait] = gearys_c(adata[strata_df.index], strata_df[trait].values)
    return rls_df_dict

#         gc_pval = get_p_from_empi_null(-gc, -np.array(gc_null_dist))
#         rls_dict.setdefault('gearys_c', []).append(gc)
#         rls_dict.setdefault('gearys_c_pval', []).append(gc_pval)
    
#             for perm_i in range(num_perms):
#                 perm_cells = np.random.choice(trs_df.index, size=ct_num_cells, replace=False)
#                 fun_null_dist.append(fun(trs_df.loc[perm_cells, 'trs_ez']))
                
#             fun_pval = get_p_from_empi_null(fun_rls, fun_null_dist)
            
#             rls_dict.setdefault(f'{fun_name}', []).append(fun_rls)
#             rls_dict.setdefault(f'{fun_name}_pval', []).append(fun_pval)

# ================================================================================
# Plotting utilities
# ================================================================================
def plot_assoc_matrix(pval_dict, pval_index, meta_df, stratify_by, fdr_level=0.2):
    """
    plot (tissue / tissue-celltype) x traits association matrix
    # Arguments:
    - pval_dict: dictionary of the p-value: trait -> array, array has been ordered
    - pval_index: index of the p-values
    - meta_df: providing metadata for `stratify_by` use
    - stratify_by: the column in `meta_df` to stratify the association.
    """
    
    def num2str(x):
        if x>1000:
            return '%0.1fk'%(x/1000)
        elif x>0:
            return '%d'%x
        else:
            return ''
    
    trait_list = list(pval_dict.keys())
    pval_df = pd.DataFrame(pval_dict, index=pval_index).join(meta_df, how="inner")
    
    assert stratify_by in ['tissue', 'tissue_celltype', 'celltype']
    stratify_list = sorted(pval_df[stratify_by].unique())
    
    # Dataframe for plotting
    df_plot = pd.DataFrame(index=stratify_list, columns=trait_list, data=0)

    for trait in trait_list:
        pval = pval_df[trait].values
        fdr = multipletests(pval, method='fdr_bh')[1]
        temp_df = pval_df.loc[fdr<fdr_level].groupby([stratify_by]).agg({stratify_by:len})
        temp_df = temp_df.loc[~temp_df[stratify_by].isna()]
        df_plot.loc[temp_df.index, trait] = temp_df[stratify_by].values

    df_plot = df_plot.loc[df_plot.max(axis=1)>10]
    df_plot = df_plot.T
    df_plot[df_plot<10] = 0
    if df_plot.size == 0:
        print('No association')
        return
    mat_annot = np.zeros(df_plot.shape, dtype=object)
    for i_col,col in enumerate(df_plot.columns):
        mat_annot[:,i_col] = [num2str(x) for x in df_plot[col].values]
    df_plot = np.log10(df_plot+1)

    plt.figure(figsize=[0.4*df_plot.shape[1]+3, 0.25*df_plot.shape[0]+3])
    sns.heatmap(df_plot, annot=mat_annot, fmt='s', cbar=False)
    plt.xticks(np.arange(df_plot.shape[1])+0.5, df_plot.columns, rotation=45, ha='right')
    plt.show()
    
def plot_qq(pval_dict, num_cols=6):
    """
    plot quantile-quantile figures in batches
    # Arguments:
    - pval_dict: dictionary of the p-value: trait -> array
    - num_cols: number of columns per row
    """
    zsc_dict = {trait: -sp.stats.norm.ppf(pval_dict[trait]).clip(min=-10,max=10) for trait in pval_dict}
    fdr_dict = {trait: multipletests(pval_dict[trait], method='fdr_bh')[1] for trait in pval_dict}

    # qqplot
    plot_trait_list = list(pval_dict.keys())

    normal_x = np.random.rand(50000)
    normal_x = -np.log10(normal_x)

    plt.figure(figsize=[20, 2+3*len(plot_trait_list) / num_cols])
    for trait_i, trait in enumerate(plot_trait_list):

        trait_logpval = -np.log10(pval_dict[trait])
        plt.subplot(int(np.ceil(len(plot_trait_list) / num_cols)), num_cols, trait_i + 1)
        qqplot(x=normal_x, y=trait_logpval, quantiles=2000, s=3, alpha=0.5)
        plt.axline((1, 1), slope=1, linestyle='--', color='k', alpha=0.5)

        plt.title('%s\nFDR<0.2: %d cells'%(trait,(fdr_dict[trait]<0.2).sum()))
        plt.ylabel('%s\nObserved[-log10(P)]'%trait)
        plt.xlabel('Expected[-log10(P)]')
    plt.tight_layout()
    plt.show()