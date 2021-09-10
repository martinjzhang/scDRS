import pandas as pd
import matplotlib
import matplotlib.transforms as mtrans
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from os.path import join
import submitit
import os
from scipy.stats import norm


def pval2zsc(pval):
    import scipy

    return -scipy.stats.norm.ppf(pval)

def zsc2pval(zsc):
    import scipy
    return 1 - scipy.stats.norm.cdf(zsc)


def agg_trs_results(df_obs, df_score, trait_list, stats, groupby="tissue_celltype"):
    """
    Extract TRS results
    df_obs: from AnnData.obs
    df_score: cell-specific score
    trait_list: list of traits
    groupby: aggregate by either tissue_celltype or cell_ontology_class
    """
    assert stats in ["pval", "mean", "sd", "fdr", "group_fdr"]
    assert groupby in ["tissue_celltype", "cell_ontology_class"]
    # extract pval or fdr
    if stats == "fdr":
        df_obs = df_obs.join(
            df_score[[x for x in df_score.columns if x.split(".")[-1] == "fdr"]]
        )
    else:
        df_obs = df_obs.join(
            df_score[[x for x in df_score.columns if x.split(".")[-1] == "pval"]]
        )

    # filter groups with small number of cells
    temp_df = df_obs.groupby(groupby).agg({"cell": len})
    temp_key_list = temp_df.index[temp_df["cell"] > 100]

    if stats == "pval":
        temp_df = df_obs.groupby(groupby).agg(
            {
                "%s.pval" % x: lambda a: np.maximum(-np.log10(len(a) * np.min(a)), 0.0)
                for x in trait_list
            }
        )
    elif stats == "mean":
        temp_df = df_obs.groupby(groupby).agg(
            {"%s.pval" % x: lambda a: np.mean(-np.log10(a)) for x in trait_list}
        )
    elif stats == "sd":
        temp_df = df_obs.groupby(groupby).agg(
            {"%s.pval" % x: lambda a: np.std(-np.log10(a)) for x in trait_list}
        )
    elif stats == "fdr":
        temp_df = df_obs.groupby(groupby).agg(
            {"%s.fdr" % x: lambda a: np.mean(a < 0.2) for x in trait_list}
        )
    else:
        raise NotImplementedError

    temp_df = temp_df.loc[temp_key_list]
    if stats == "fdr":
        temp_df.columns = [x.replace(f".fdr", "") for x in temp_df.columns]
    else:
        temp_df.columns = [x.replace(f".pval", "") for x in temp_df.columns]
    trs_df = temp_df.T.copy()
    trs_df.columns = [c.replace(" ", "_").replace(",", "") for c in trs_df.columns]
    return trs_df


from matplotlib import colors

def discrete_cmap(N, base_cmap=None, start_white=True):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    if start_white:
        color_list[0, :] = 1.
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def my_plot_heatmap(
    df,
    dpi=150,
    squaresize=20,
    heatmap_linewidths=0.5,
    heatmap_linecolor="gray",
    heatmap_xticklabels=True,
    heatmap_yticklabels=True,
    heatmap_cbar=True,
    heatmap_cbar_kws=dict(use_gridspec=False, location="top", fraction=0.03, pad=0.01),
    heatmap_vmin=0.0,
    heatmap_vmax=1.0,
    xticklabels_rotation=45,
    colormap_n_bin=10
):
    figwidth = df.shape[1] * squaresize / float(dpi)
    figheight = df.shape[0] * squaresize / float(dpi)
    fig, ax = plt.subplots(1, figsize=(figwidth, figheight), dpi=dpi)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    # discrete_cmap(10, "Reds")
    # colors.LinearSegmentedColormap.from_list("my_reds", [(1, 1, 1), (1, 0, 0)], N=10),
    sns.heatmap(
        df,
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
        ax.get_xticklabels(), rotation=xticklabels_rotation, va="top", ha="right", fontsize=8
    )
    ax.tick_params(left=False, bottom=False, pad=-2)
    trans = mtrans.Affine2D().translate(5, 0)
    for t in ax.get_xticklabels():
        t.set_transform(t.get_transform() + trans)
    return fig, ax


def cluster_viz(
    df1,
    df2=None,
    name1=None,
    name2=None,
    n_cluster=3,
    vmin=0,
    vmax=7,
    figsize_scale=(0.3, 0.3),
):
    """
    Use df1 to determine the cluster, and visualize both df1 and df2
    """
    # cluster with df1
    model = SpectralCoclustering(n_clusters=n_cluster, random_state=0)
    model.fit(df1.values)
    plot_ct_list = df1.columns[np.argsort(model.column_labels_)]

    # visualize df1 results
    plot_df1 = df1.loc[:, plot_ct_list]
    plt.figure(
        figsize=[
            figsize_scale[0] * plot_df1.shape[1],
            figsize_scale[1] * plot_df1.shape[0],
        ]
    )
    sns.heatmap(
        plot_df1,
        cmap="Reds",
        cbar=True,
        cbar_kws={"fraction": 0.01},
        vmin=vmin,
        vmax=vmax,
        linewidths=0.1,
        linecolor="gray",
    )
    plt.title(name1)
    plt.tight_layout()
    plt.show()

    if df2 is not None:
        # visualize df2 results
        plot_df2 = df2.loc[:, plot_ct_list]
        plt.figure(
            figsize=[
                figsize_scale[0] * plot_df2.shape[1],
                figsize_scale[1] * plot_df2.shape[0],
            ]
        )
        sns.heatmap(
            plot_df2,
            cmap="Reds",
            cbar=True,
            cbar_kws={"fraction": 0.01},
            vmin=vmin,
            vmax=vmax,
            linewidths=0.1,
            linecolor="gray",
        )
        plt.title(name2)
        plt.tight_layout()
        plt.show()


def viz(
    df1,
    df2,
    name1=None,
    name2=None,
    n_cluster=3,
    vmin=0,
    vmax=7,
    figsize_scale=(0.3, 0.3),
):
    """
    Use df1 to determine the cluster, and visualize both df1 and df2
    """

    fig, axes = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=[figsize_scale[0] * df1.shape[1], figsize_scale[1] * df1.shape[0]],
        sharex=True,
        dpi=125,
    )
    h = sns.heatmap(
        df1,
        cmap="Reds",
        cbar=True,
        cbar_kws={"fraction": 0.01},
        xticklabels=True,
        yticklabels=True,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.1,
        linecolor="gray",
        ax=axes[0],
    )
    h.set_yticklabels(h.get_yticklabels(), size=8)
    axes[0].set_title(name1)
    h = sns.heatmap(
        df2,
        cmap="Reds",
        cbar=True,
        cbar_kws={"fraction": 0.01},
        xticklabels=True,
        yticklabels=True,
        vmin=vmin,
        vmax=vmax,
        linewidths=0.1,
        linecolor="gray",
        ax=axes[1],
    )
    h.set_yticklabels(h.get_yticklabels(), size=8)
    axes[1].set_title(name2)
    fig.tight_layout()
    plt.show()

    
def plot_diagonal_block(xs, ys, ax, linewidth=1, linestyle="--", color="black"):
    assert len(xs) == len(ys)
    start_x = 0
    start_y = 0
    for x, y in zip(xs, ys):
        ax.hlines(
            y=start_y,
            xmin=start_x,
            xmax=start_x + x,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
        )
        ax.hlines(
            y=start_y + y,
            xmin=start_x,
            xmax=start_x + x,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
        )
        ax.vlines(
            x=start_x,
            ymin=start_y,
            ymax=start_y + y,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
        )
        ax.vlines(
            x=start_x + x,
            ymin=start_y,
            ymax=start_y + y,
            linewidth=linewidth,
            linestyle=linestyle,
            color=color,
        )
        start_x += x
        start_y += y

from statsmodels.stats.multitest import multipletests
def agg_trs_pval(df_obs, df_pval, stats, groupby, fdr_prop_threshold=None, pval_group_quantile=None):
    """
    Extract TRS results
    df_obs: from AnnData.obs
    df_score: cell-specific score
    trait_list: list of traits
    groupby: aggregate by either tissue_celltype or cell_ontology_class
    """
    assert stats in ["zsc_mean", "zsc_sd", "fdr_prop", "pval_group"]

    list_trait = df_pval.columns
    # extract pval or fdr
    if stats.startswith("zsc"):
        df_zsc = df_pval.transform(lambda x: pval2zsc(x))
        df_obs = df_obs.join(df_zsc)
        if stats == "zsc_mean":
            df_obs = df_obs.groupby(groupby).agg(
                {trait: lambda x: np.mean(x) for trait in list_trait}
            )
        elif stats == "zsc_sd":
            df_obs = df_obs.groupby(groupby).agg(
                {trait: lambda x: np.std(x) for trait in list_trait}
            )
        else:
            raise NotImplementedError

    elif stats == "fdr_prop":
        df_fdr = df_pval.transform(lambda x: multipletests(x, method="fdr_bh")[1])
        assert fdr_prop_threshold is not None
        df_obs = df_obs.join(df_fdr)
        df_obs = df_obs.groupby(groupby).agg(
            {trait: lambda x: np.mean(x < fdr_prop_threshold) for trait in list_trait}
        )
    
    elif stats == "pval_group":
        assert pval_group_quantile is not None
        df_obs = df_obs.join(df_pval)
        df_obs = df_obs.groupby(groupby).agg(
            {trait: lambda x : np.maximum(-np.log10(1 / pval_group_quantile * np.quantile(x, pval_group_quantile)), 0.0) for trait in list_trait}
        )

    else:
        raise NotImplementedError
    df_obs = df_obs.T.copy()
    df_obs.columns = [c.replace(" ", "_").replace(",", "") for c in df_obs.columns]
    return df_obs