import numpy as np
import scipy as sp
import pandas as pd
import numbers
import matplotlib.pyplot as plt

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
    