import numpy as np
import scipy as sp
import pandas as pd
import numbers
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import scanpy as sc
import os
import anndata
import matplotlib.transforms as mtrans
import matplotlib.pyplot as plt
from typing import List, Dict
from anndata import read_h5ad
import fnmatch
import matplotlib.patches as patches


def convert_species_name(species):
    if species in ["Mouse", "mouse", "Mus_musculus", "mus_musculus", "mmusculus"]:
        return "mmusculus"
    if species in ["Human", "human", "Homo_sapiens", "homo_sapiens", "hsapiens"]:
        return "hsapiens"
    raise ValueError("species name '%s' is not supported" % species)


def str_or_list_like(x):
    """Determine if x is list-list (list, tuple) using duck-typing.

    Here is a set of Attributes for different classes
    |         x          |      type(x)       |    x.strip    | x.__getitem__ |  x.__iter__   |
    |         aa         |   <class 'str'>    |     True      |     True      |     True      |
    |     ['a', 'b']     |   <class 'list'>   |     False     |     True      |     True      |
    |     ('a', 'b')     |  <class 'tuple'>   |     False     |     True      |     True      |
    |     {'b', 'a'}     |   <class 'set'>    |     False     |     False     |     True      |
    |  {'a': 1, 'b': 2}  |   <class 'dict'>   |     False     |     True      |     True      |
    """

    if hasattr(x, "strip"):
        return "str"
    elif hasattr(x, "__getitem__") or hasattr(x, "__iter__"):
        return "list_like"
    else:
        return "others"


def load_h5ad(
    h5ad_file: str, flag_filter_data: bool = False, flag_raw_count: bool = True
) -> anndata.AnnData:
    """Load h5ad file and optionally filter out cells and perform normalization.

    Parameters
    ----------
    h5ad_file : str
        Path to h5ad file
    flag_filter_data : bool
        If True, filter out cells with
        
        - sc.pp.filter_cells(adata, min_genes=250)
        - sc.pp.filter_genes(adata, min_cells=50)
    flag_raw_count : bool
        If True, perform size-factor normalization and log1p transformation.

    Returns
    -------    
    adata : anndata.AnnData
        Single-cell data.
    """
    adata = read_h5ad(h5ad_file)
    if flag_filter_data:
        sc.pp.filter_cells(adata, min_genes=250)
        sc.pp.filter_genes(adata, min_cells=50)
    if flag_raw_count:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
    return adata


def load_scdrs_score(
    score_file: str, obs_names: List[str] = None
) -> Dict[str, pd.DataFrame]:
    """Load scDRS scores. 
    
    Use "@" to specify multiple files, e.g., `score_folder/@.full_score.gz`

    Parameters
    ----------
    score_file : str
        Path to scDRS `.full_score.gz` file. Use '@' to specify multiple file names,
        e.g., `score_folder/@.full_score.gz`. However, `score_folder` should
        not contain '@'.
    obs_names : List[str]
        Expected list of cells. Score files with less than 10% overlap with this list
        will be skipped.

    Returns
    -------    
    dict_score : Dict[str, pd.DataFrame]
        Dictionary of scDRS full score DataFrames, keyed by trait name.
    """
    assert score_file.endswith(
        "full_score.gz"
    ), "Expect scDRS .full_score.gz files for score_file"

    # Get score_dir and score_file_list for potentially multiple score files
    score_file_pattern = score_file.split(os.path.sep)[-1]
    score_dir = score_file.replace(os.path.sep + score_file_pattern, "")
    score_file_list = [
        x
        for x in os.listdir(score_dir)
        if fnmatch.fnmatch(x, score_file_pattern.replace("@", "*"))
    ]
    print("Score file folder: %s" % score_dir)
    print("Find %d score files: %s" % (len(score_file_list), ",".join(score_file_list)))

    dict_score = {}
    for score_file in score_file_list:
        temp_df = pd.read_csv(
            score_dir + os.path.sep + score_file, sep="\t", index_col=0
        )
        if obs_names is not None:
            # Check overlap of cells between score_file and obs_names
            n_cell_overlap = len(set(obs_names) & set(temp_df.index))
            if n_cell_overlap < 0.1 * len(obs_names):
                print(
                    "WARNING: %s skipped, containing %d/%d target cells"
                    % (score_file, n_cell_overlap, len(obs_names))
                )
                continue
        dict_score[score_file.replace(".full_score.gz", "")] = temp_df.copy()

    return dict_score


def load_homolog_mapping(src_species: str, dst_species: str) -> dict:
    """Load gene homologs between mouse and human.

    Parameters
    ----------
    src_species : str
        Source species. One of 'mmusculus', 'mouse', 'hsapiens', or 'human'.
    dst_species : str
        Destination species. One of 'mmusculus', 'mouse', 'hsapiens', or 'human'.
        Cannot be the same as `src_species`.

    Returns
    -------    
    dic_map : dict
        Dictionary of gene homologs (gene symbol).
    """

    src_species = convert_species_name(src_species)
    dst_species = convert_species_name(dst_species)

    assert src_species != dst_species, "src and dst cannot be the same"

    df_hom = pd.read_csv(
        os.path.join(os.path.dirname(__file__), "data/mouse_human_homologs.txt"),
        sep="\t",
    )
    if (src_species == "hsapiens") & (dst_species == "mmusculus"):
        dic_map = {
            x: y for x, y in zip(df_hom["HUMAN_GENE_SYM"], df_hom["MOUSE_GENE_SYM"])
        }
    elif (src_species == "mmusculus") & (dst_species == "hsapiens"):
        dic_map = {
            x: y for x, y in zip(df_hom["MOUSE_GENE_SYM"], df_hom["HUMAN_GENE_SYM"])
        }
    else:
        raise NotImplementedError(
            f"gene conversion from {src_species} to {dst_species} is not supported"
        )

    return dic_map


def load_gs(
    gs_path: str,
    src_species: str = None,
    dst_species: str = None,
    to_intersect: List[str] = None,
) -> dict:
    """Load the gene set file (.gs file).

    Parameters
    ----------
    gs_path : str
        Path to the gene set file with the following two columns, separated by tab:
        
        - 'TRAIT'
        - 'GENESET':
          (1) <gene1>,<gene2>,... each gene will be weighted uniformly or
          (2) <gene1>:<weight1>,<gene2>:<weight2>,... each gene will be weighted by its weight.
    src_species : str, default=None
        Source species, must be either 'mmusculus' or 'hsapiens' if not None
    dst_species : str, default=None
        Destination species, must be either 'mmusculus' or 'hsapiens' if not None
    to_intersect : List[str], default None.
        Gene list to intersect with the input .gs file.


    Returns
    -------    
    dict_gs : dict
        Dictionary of gene sets: {
            trait1: (gene_list, gene_weight_list),
            trait2: (gene_list, gene_weight_list),
            ...
        }
    """

    assert (src_species is None) == (
        dst_species is None
    ), "src_species and dst_species must be both None or not None"

    # Load homolog map dict_map; only needed when src_species and dst_species
    # are not None and different.
    if ((src_species is not None) & (dst_species is not None)) and (
        src_species != dst_species
    ):
        dict_map = load_homolog_mapping(src_species, dst_species)  # type: ignore
    else:
        dict_map = None  # type: ignore

    # Load gene set file
    dict_gs = {}
    df_gs = pd.read_csv(gs_path, sep="\t")
    for i, (trait, gs) in df_gs.iterrows():
        gs_info = [g.split(":") for g in gs.split(",")]
        if np.all([len(g) == 1 for g in gs_info]):
            # if all genes are weighted uniformly
            dict_weights = {g[0]: 1.0 for g in gs_info}
        elif np.all([len(g) == 2 for g in gs_info]):
            # if all genes are weighted by their weights
            dict_weights = {g[0]: float(g[1]) for g in gs_info}
        else:
            raise ValueError(f"gene set {trait} contains genes with invalid format")

        # Convert species if needed
        # convert gene list to homologs, if gene can not be mapped, remove it
        # in both gene list and gene weight
        if dict_map is not None:
            dict_weights = {
                dict_map[g]: w for g, w in dict_weights.items() if g in dict_map
            }

        # Intersect with other gene sets
        if to_intersect is not None:
            to_intersect = set(to_intersect)
            dict_weights = {g: w for g, w in dict_weights.items() if g in to_intersect}

        gene_list = list(dict_weights.keys())
        dict_gs[trait] = (
            gene_list,
            [dict_weights[g] for g in gene_list],
        )

    return dict_gs


def test_overlap(list1, list2, list_background):
    """
    Test overlap of two gene sets using Fisher's exact test
    """

    set1 = set(list1)
    set2 = set(list2)
    set_background = set(list_background)

    n1 = len(set1)
    n2 = len(set2)
    n_overlap = len(set1 & set2)
    n_other = len(set_background - set1 - set2)

    oddsratio, pvalue = sp.stats.fisher_exact(
        [[n_other, n1 - n_overlap], [n2 - n_overlap, n_overlap]]
    )

    if (
        (n_overlap == 0)
        | (n_other == 0)
        | ((n2 - n_overlap) == 0)
        | ((n1 - n_overlap) == 0)
    ):
        return pvalue, oddsratio, 0, 0
    else:
        se_log_or = np.sqrt(
            1 / (n1 - n_overlap) + 1 / (n2 - n_overlap) + 1 / n_overlap + 1 / n_other
        )
        or_ub = np.exp(np.log(oddsratio) + 1.96 * se_log_or)
        or_lb = np.exp(np.log(oddsratio) - 1.96 * se_log_or)
        return pvalue, oddsratio, or_ub, or_lb
    

def meta_analysis(effects, se, method='random', weights=None):
    """
    Random effect meta-analysis. From Omer Weissbrod
    """

    assert method in ['fixed', 'random']
    d = effects
    variances = se**2
    
    #compute random-effects variance tau2
    vwts = 1.0 / variances
    fixedsumm = vwts.dot(d) / vwts.sum()    
    Q = np.sum(((d - fixedsumm)**2) / variances)
    df = len(d)-1
    tau2 = np.maximum(0, (Q-df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))
    
    #defing weights
    if weights is None:
        if method == 'fixed':
            wt = 1.0 / variances
        else:
            wt = 1.0 / (variances + tau2)
    else:
        wt = weights
    
    #compute summtest
    summ = wt.dot(d) / wt.sum()
    if method == 'fixed':
        varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
    else:
        varsum = np.sum(wt*wt*(variances+tau2)) / (np.sum(wt)**2)
    ###summtest = summ / np.sqrt(varsum)
    
    summary=summ
    se_summary=np.sqrt(varsum)
    
    return summary, se_summary


# https://stats.stackexchange.com/questions/403652/two-sample-quantile-quantile-plot-in-python
def qqplot(x, y, quantiles=None, interpolation="nearest", ax=None, **kwargs):
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


def zsc2pval(zsc):
    """
    Convert z-score to one-sided p-value. Accurate up to `zsc=36` and `pval=4.2e-284`.
    """
    #     return 1 - sp.stats.norm.cdf(zsc)
    return sp.stats.norm.cdf(-zsc)  # This is more accurate


def pval2zsc(pval):
    """
    Convert one-sided p-value to z-score. Accurate up to `zsc=36` and `pval=4.2e-284`.
    """
    return -sp.stats.norm.ppf(pval)


# ================================================================================
# Plotting utilities
# ================================================================================
def small_squares(ax, pos, size=1, linewidth=0.8):
    """
    Draw many small squares on ax, given the positions of
    these squares.

    """
    for xy in pos:
        x, y = xy
        margin = (1 - size) / 2
        rect = patches.Rectangle(
            (x + margin, y + margin),
            size,
            size,
            linewidth=linewidth,
            edgecolor="k",
            facecolor="none",
            zorder=20,
        )
        ax.add_patch(rect)


def plot_group_stats(
    dict_df_stats=None,
    df_fdr_prop=None,
    df_assoc_fdr=None,
    df_hetero_fdr=None,
):
    """plot group-level statistics for scDRS results

    Parameters
    ----------
    dict_df_stats: dict trait -> pd.DataFrame
        group-level statistics from `group_stats` function for each trait
    df_fdr_prop : pd.DataFrame
        dataframe of proportion of cells with FDR < 0.1
    df_assoc : pd.DataFrame
        dataframe of group-level association statistics
    df_hetero : pd.DataFrame
        dataframe of group-level heterogeneity statistics
    """
    if dict_df_stats is not None:
        trait_list = list(dict_df_stats.keys())
        # compile df_fdr_prop, df_assoc_fdr, df_hetero_fdr from dict_df_stats
        df_fdr_prop = pd.concat(
            [dict_df_stats[trait]["fdr_prop"] for trait in trait_list], axis=1
        ).T
        df_assoc_fdr = pd.concat(
            [dict_df_stats[trait]["assoc_pval"] for trait in trait_list], axis=1
        ).T
        df_assoc_fdr = pd.DataFrame(
            multipletests(df_assoc_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_assoc_fdr.shape
            ),
            index=df_assoc_fdr.index,
            columns=df_assoc_fdr.columns,
        )
        df_hetero_fdr = pd.concat(
            [dict_df_stats[trait]["hetero_pval"] for trait in trait_list], axis=1
        ).T
        df_hetero_fdr = pd.DataFrame(
            multipletests(df_hetero_fdr.values.flatten(), method="fdr_bh")[1].reshape(
                df_hetero_fdr.shape
            ),
            index=df_hetero_fdr.index,
            columns=df_hetero_fdr.columns,
        )
        df_fdr_prop.index = trait_list
        df_assoc_fdr.index = trait_list
        df_hetero_fdr.index = trait_list

    df_hetero_fdr = df_hetero_fdr.applymap(lambda x: "×" if x < 0.05 else "")
    df_hetero_fdr[df_assoc_fdr > 0.05] = ""

    fig, ax = plot_heatmap(
        df_fdr_prop,
        squaresize=30,
        heatmap_annot=df_hetero_fdr,
        heatmap_annot_kws={"color": "black", "size": 8},
        heatmap_cbar_kws=dict(
            use_gridspec=False, location="top", fraction=0.1, pad=0.05, drawedges=True
        ),
        heatmap_vmin=0,
        heatmap_vmax=0.2,
        colormap_n_bin=5,
    )

    small_squares(
        ax,
        pos=[(y, x) for x, y in zip(*np.where(df_assoc_fdr < 0.05))],
        size=0.6,
        linewidth=0.5,
    )

    cb = ax.collections[0].colorbar
    cb.ax.tick_params(labelsize=8)

    cb.ax.set_title("Prop. of sig. cells (FDR < 0.1)", fontsize=8)
    cb.outline.set_edgecolor("black")
    cb.outline.set_linewidth(1)

    plt.tight_layout()


def discrete_cmap(N, base_cmap=None, start_white=True):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    if start_white:
        color_list[0, :] = 1.0
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def plot_heatmap(
    df,
    dpi=150,
    squaresize=20,
    heatmap_annot=None,
    heatmap_annot_kws={"color": "black", "size": 4},
    heatmap_linewidths=0.5,
    heatmap_linecolor="gray",
    heatmap_xticklabels=True,
    heatmap_yticklabels=True,
    heatmap_cbar=True,
    heatmap_cbar_kws=dict(use_gridspec=False, location="top", fraction=0.03, pad=0.01),
    heatmap_vmin=0.0,
    heatmap_vmax=1.0,
    xticklabels_rotation=45,
    colormap_n_bin=10,
):
    figwidth = df.shape[1] * squaresize / float(dpi)
    figheight = df.shape[0] * squaresize / float(dpi)
    fig, ax = plt.subplots(1, figsize=(figwidth, figheight), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_facecolor("silver")
    sns.heatmap(
        df,
        annot=heatmap_annot,
        annot_kws=heatmap_annot_kws,
        fmt="",
        cmap=discrete_cmap(colormap_n_bin, "Reds"),
        linewidths=heatmap_linewidths,
        linecolor=heatmap_linecolor,
        square=True,
        ax=ax,
        xticklabels=heatmap_xticklabels,
        yticklabels=heatmap_yticklabels,
        cbar=heatmap_cbar,
        cbar_kws=heatmap_cbar_kws,
        vmin=heatmap_vmin,
        vmax=heatmap_vmax,
    )

    plt.yticks(fontsize=8)
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=xticklabels_rotation,
        va="top",
        ha="right",
        fontsize=8,
    )
    ax.tick_params(left=False, bottom=False, pad=-2)
    trans = mtrans.Affine2D().translate(5, 0)
    for t in ax.get_xticklabels():
        t.set_transform(t.get_transform() + trans)
    return fig, ax


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
        if x > 1000:
            return "%0.1fk" % (x / 1000)
        elif x > 0:
            return "%d" % x
        else:
            return ""

    trait_list = list(pval_dict.keys())
    pval_df = pd.DataFrame(pval_dict, index=pval_index).join(meta_df, how="inner")

    assert stratify_by in ["tissue", "tissue_celltype", "celltype"]
    stratify_list = sorted(pval_df[stratify_by].unique())

    # Dataframe for plotting
    df_plot = pd.DataFrame(index=stratify_list, columns=trait_list, data=0)

    for trait in trait_list:
        pval = pval_df[trait].values
        fdr = multipletests(pval, method="fdr_bh")[1]
        temp_df = (
            pval_df.loc[fdr < fdr_level].groupby([stratify_by]).agg({stratify_by: len})
        )
        temp_df = temp_df.loc[~temp_df[stratify_by].isna()]
        df_plot.loc[temp_df.index, trait] = temp_df[stratify_by].values

    df_plot = df_plot.loc[df_plot.max(axis=1) > 10]
    df_plot = df_plot.T
    df_plot[df_plot < 10] = 0
    if df_plot.size == 0:
        print("No association")
        return
    mat_annot = np.zeros(df_plot.shape, dtype=object)
    for i_col, col in enumerate(df_plot.columns):
        mat_annot[:, i_col] = [num2str(x) for x in df_plot[col].values]
    df_plot = np.log10(df_plot + 1)

    plt.figure(figsize=[0.4 * df_plot.shape[1] + 3, 0.25 * df_plot.shape[0] + 3])
    sns.heatmap(df_plot, annot=mat_annot, fmt="s", cbar=False)
    plt.xticks(
        np.arange(df_plot.shape[1]) + 0.5, df_plot.columns, rotation=45, ha="right"
    )
    plt.show()


def plot_qq(pval_dict, num_cols=6):
    """
    plot quantile-quantile figures in batches
    Arguments
    ---------
    - pval_dict: dictionary of the p-value: trait -> array
    - num_cols: number of columns per row
    """
    zsc_dict = {
        trait: -sp.stats.norm.ppf(pval_dict[trait]).clip(min=-10, max=10)
        for trait in pval_dict
    }
    fdr_dict = {
        trait: multipletests(pval_dict[trait], method="fdr_bh")[1]
        for trait in pval_dict
    }

    # qqplot
    plot_trait_list = list(pval_dict.keys())

    normal_x = np.random.rand(50000)
    normal_x = -np.log10(normal_x)

    plt.figure(figsize=[20, 2 + 3 * len(plot_trait_list) / num_cols])
    for trait_i, trait in enumerate(plot_trait_list):

        trait_logpval = -np.log10(pval_dict[trait])
        plt.subplot(
            int(np.ceil(len(plot_trait_list) / num_cols)), num_cols, trait_i + 1
        )
        qqplot(x=normal_x, y=trait_logpval, quantiles=2000, s=3, alpha=0.5)
        plt.axline((1, 1), slope=1, linestyle="--", color="k", alpha=0.5)

        plt.title("%s\nFDR<0.2: %d cells" % (trait, (fdr_dict[trait] < 0.2).sum()))
        plt.ylabel("%s\nObserved[-log10(P)]" % trait)
        plt.xlabel("Expected[-log10(P)]")
    plt.tight_layout()
    plt.show()


def plot_score_umap(
    score_dict, score_index, umap_adata, umap_color=["cell_ontology_class"], n_col=5
):
    """
    Overlay score on UMAP
    ---
    Args
    score_dict: label -> list of scores
    score_index: index of the scores in `score_dict`
    umap_adata: Anndata containing the umap
    umap_color: which attributes to plot before hand
    n_col: number of columns per row
    """
    umap_adata = umap_adata.copy()

    sc.pl.umap(umap_adata, color=umap_color, size=20, ncols=1)
    df_plot = pd.DataFrame(index=umap_adata.obs.index)
    df_plot["UMAP1"] = umap_adata.obsm["X_umap"][:, 0]
    df_plot["UMAP2"] = umap_adata.obsm["X_umap"][:, 1]
    df_plot = df_plot.join(pd.DataFrame(score_dict, index=score_index))

    # Trait TRS plot
    plt.figure(figsize=[15, 2 + 3 * len(score_dict) / n_col])

    for trait_i, trait in enumerate(score_dict.keys()):
        plt.subplot(int(np.ceil(len(score_dict) / n_col)), n_col, trait_i + 1)
        # max_ = np.quantile(np.absolute(df_plot[trait].values), 0.99)
        # min_ = np.quantile(np.absolute(df_plot[trait].values), 0.01)
        plt.scatter(
            df_plot["UMAP1"], df_plot["UMAP2"], c=df_plot[trait], cmap="RdBu_r", s=4
        )
        plt.colorbar()
        plt.clim(-4, 4)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title("%s" % trait)
    plt.tight_layout()
    plt.show()


def p_2_str(p_):
    if p_ > 0.05:
        return ""
    elif p_ > 0.005:
        return "*"
    else:
        return "**"


def p_2_str_num(p_, n_ctrl):
    if p_ > 1 / (n_ctrl + 0.5):
        return "P=%0.3f" % p_
    else:
        return "P<%0.3f" % (1 / n_ctrl)
    
    
def reorder_col(df_sim):
    mat_condensed_dist = sp.spatial.distance.squareform(1-df_sim)
    mat_linkage = sp.cluster.hierarchy.linkage(mat_condensed_dist, method='average')
    reordered_col_list = list(df_sim.columns[sp.cluster.hierarchy.leaves_list(mat_linkage)])
    return reordered_col_list
