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


def check_import():
    print("# test: scdrs")
    return


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


def convert_gs_species(
    df_gs: pd.DataFrame, src: str = "hsapiens", dst: str = "mmusculus"
) -> pd.DataFrame:
    """Convert species of genes in a gene set dataframe

    Parameters
    ----------
    gs_species : pd.DataFrame
        Gene set dataframe with columns:
        - 'TRAIT'
        - 'GENESET'
    src : str, default 'hsapiens'
        Source species, must be either 'mmusculus' or 'hsapiens'
    dst : str, default 'mmusculus'
        Destination species, must be either 'mmusculus' or 'hsapiens'


    Returns
    -------
    gs_species_converted : pd.DataFrame
    """

    assert src in [
        "mmusculus",
        "hsapiens",
    ], "src must be either 'mmusculus' or 'hsapiens'"
    assert dst in [
        "mmusculus",
        "hsapiens",
    ], "dst must be either 'mmusculus' or 'hsapiens'"
    assert src != dst, "src and dst cannot be the same"

    dirname = os.path.dirname(__file__)
    df_hom = pd.read_csv(
        os.path.join(dirname, "data/mouse_human_homologs.txt"), sep="\t"
    )
    if (src == "hsapiens") & (dst == "mmusculus"):
        dic_map = {
            x: y for x, y in zip(df_hom["HUMAN_GENE_SYM"], df_hom["MOUSE_GENE_SYM"])
        }
    elif (src == "mmusculus") & (dst == "hsapiens"):
        dic_map = {
            x: y for x, y in zip(df_hom["MOUSE_GENE_SYM"], df_hom["HUMAN_GENE_SYM"])
        }
    else:
        raise ValueError(
            "# gene conversion from %s to %s is not supported" % (src, dst)
        )
    df_gs = df_gs.copy()
    for trait in df_gs.index:
        src_gene_list = df_gs.loc[trait, "GENESET"].split(",")
        dst_gene_list = [dic_map[x] for x in set(src_gene_list) & set(dic_map.keys())]
        df_gs.loc[trait, "GENESET"] = ",".join(dst_gene_list)
    return df_gs


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


def test_gearysc(
    adata: anndata.AnnData,
    df_score_full: pd.DataFrame,
    groupby: str,
    opt="control_distribution_match",
) -> pd.DataFrame:
    """
    Compute significance level for Geary's C statistics

    Args
    ----
    adata: anndata.AnnData
        must contain `connectivities` to compute the Geary's C statistic
    df_score_full: DataFrame
        DataFrame with the scores of the cells, contains
        columns `zscore`, `norm_score`, `ctrl_norm_score_{i}`
    groupby: str
        Column name of the groupby variable.
    opt: str
        Options:
            - "control_distribution_match":
                The distribution of the scores of the control scores is similar to
                the distribution of the scores of the disease scores.

    Returns
    -------
    df_rls: DataFrame
        DataFrame with the results of the test with `n_group` rows and 4 columns:
            - `pval`: significance level of Geary's C
            - `trait`: Geary's C test statistic of the trait scores
            - `ctrl_mean`: mean of the control scores
            - `ctrl_sd`: standard deviation of the control scores
    """

    df_score_full = df_score_full.reindex(adata.obs.index).dropna()
    norm_score = df_score_full["norm_score"]
    ctrl_norm_score = df_score_full[
        [col for col in df_score_full.columns if col.startswith(f"ctrl_norm_score_")]
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


def zsc2pval(zsc):
    return 1 - sp.stats.norm.cdf(zsc)


def pval2zsc(pval):
    return -sp.stats.norm.ppf(pval)


# ================================================================================
# Plotting utilities
# ================================================================================
import matplotlib.patches as patches


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
    # Arguments:
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
