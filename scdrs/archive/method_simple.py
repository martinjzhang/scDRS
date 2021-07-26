# data
import scanpy as sc
import pandas as pd
import numpy as np
import scipy as sp

import os
import time
from os.path import join

# stats
from statsmodels.stats.multitest import multipletests

# plotting
import matplotlib.pyplot as plt
import seaborn as sns

import itertools
from gprofiler import GProfiler
import submitit
import itertools

# scTRS tools
import scTRS.util as util
import scTRS.data_loader as dl
import scTRS.method as md

from scTRS.method import get_p_from_empi_null
from IPython.display import Markdown, display


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def load_score_dataset(name, filter_homolog=True):
    DATA_PATH = '/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data'
    assert name in ['tms_facs', 'tms_droplet', 'aizarani']
    if name in ['tms_facs', 'tms_droplet']:
        if name == 'tms_facs':
            data = dl.load_tms_ct(DATA_PATH, data_name='facs')
        elif name == 'tms_droplet':
            data = dl.load_tms_ct(DATA_PATH, data_name='droplet')
        else:
            raise NotImplementedError
            
        if filter_homolog:
            GENE_SCORE_PATH = join(DATA_PATH, 'trs_gene_scores')
            hsapiens_mmusculus_mapping = pd.read_csv(join(GENE_SCORE_PATH, 'meta_data', 'hsapiens_mmusculus_mapping.csv'))
            hsapiens_genes = sorted(list(set(hsapiens_mmusculus_mapping['mmusculus'].values) & set(data.var_names)))
            data = data[:, hsapiens_genes].copy()
        return data
    elif name == 'aizarani':
        data = dl.load_aizarani_raw_data(opt='processed')
        # TODO: more filtering?
        return data
    else:
        raise NotImplementedError

        
def _compute_trs(adata, gene_list, gene_weight):

    gene_list = list(gene_list)
    # the gene_weight includes both vst weight and user-specified weights.
    gene_weight /= gene_weight.sum()
    temp_v = adata[:, gene_list].X.dot(gene_weight) - adata[:, gene_list].var['mean'].dot(gene_weight)
    v_trs = np.array(temp_v, dtype=np.float64).reshape([-1])
        
    return v_trs

def _select_ctrl_genes(trait_genes, all_gene_df, bin_name, num_ctrl, num_bins):
    # calculate the bin statistics
    df = all_gene_df[[bin_name]].reset_index().rename(columns={'index': 'gene', bin_name: 'value'}).copy()
    df['bin'] = pd.qcut(df['value'], q=num_bins, duplicates='drop')
    df = df.groupby('bin').agg({'gene': set})
    
    # calculate number of overlap
    overlap_num = df['gene'].apply(lambda x: len(x & set(trait_genes)))

    ctrl_genes = [[] for i in range(num_ctrl)]
    for index, bin_row in df.iterrows():
        bin_genes = sorted(list(bin_row['gene']))
        if overlap_num.loc[index]==0: continue
        for ctrl_i in range(num_ctrl):
            ctrl_genes[ctrl_i].extend(np.random.choice(bin_genes, size=overlap_num.loc[index], replace=False))
    return ctrl_genes
    
    
def calculate_norm_factor(v, robust=True, axis=None):
    if robust:
        med = np.median(v, axis=axis)
        q1, q3 = np.percentile(v, (25, 75), axis=axis)
        return med, q3 - q1
    else:
        mean = np.mean(v, axis=axis)
        std = np.std(v, axis=axis)
        return mean, std
    

def score_cell(dataset, gene_list, num_ctrl = 500, num_bins=200, random_state=1234):
    """
    dataset: for what dataset to derive TRS ['tms_facs', random, liver] Either a string or anndata
    gene_weight: use what kind of gene weights [all, subset]
    permutation_gene: how to select the genes for permutation [random, whole_mean_match, subset_mean_match, whole_weighted_mean_match, subset_weighted_mean_match]
    permutation_weight: how the weight for permutation genes is set [trait, own]
    
    """
    np.random.seed(random_state)
    if isinstance(dataset, str):
        score_dataset = load_score_dataset(name=dataset)
    else:
        score_dataset = dataset.copy()
    md.compute_stats(score_dataset)
    
    all_gene_df = score_dataset.var[['mean', 'var', 'var_tech']].copy()
    all_gene_df['weight'] = 1. / np.sqrt(all_gene_df['var_tech'].values.clip(min=1e-2))

    # construct trait_gene_df and control_gene_df
    trait_gene_df = all_gene_df.reindex(gene_list).dropna()
    
    ctrl_genes_list = _select_ctrl_genes(gene_list, all_gene_df, bin_name='mean', num_ctrl=num_ctrl, num_bins=num_bins)
    ctrl_gene_df_list = [all_gene_df[['weight']].loc[ctrl_genes, :].copy() for ctrl_genes in ctrl_genes_list]

    # compute TRS for both trait gene set and control gene sets
    trait_trs = _compute_trs(score_dataset, trait_gene_df.index.values, trait_gene_df.weight.values)

    ctrl_trs_list = []
    for ctrl_gene_df in ctrl_gene_df_list: 
        ctrl_trs_list.append(_compute_trs(score_dataset, ctrl_gene_df.index.values, ctrl_gene_df.weight.values))
    ctrl_trs = np.vstack(ctrl_trs_list)
    
    def _normalize_across_cell(trait_trs, ctrl_trs):
        """
        trait_trs: (#cells, ) array
        ctrl_trs: (#controls, #cells) array
        """
        # TODO: add a normalization within tissue option
        ctrl_trs_mean, ctrl_trs_sd = np.mean(ctrl_trs, axis=0), np.std(ctrl_trs, axis=0)
        trait_zsc = (trait_trs - ctrl_trs_mean) / ctrl_trs_sd
        ctrl_zsc = (ctrl_trs - ctrl_trs_mean) / ctrl_trs_sd
        return trait_zsc, ctrl_zsc

    
    trait_zsc, ctrl_zsc = _normalize_across_cell(trait_trs, ctrl_trs)
    trait_ep = get_p_from_empi_null(trait_zsc, ctrl_zsc.flatten())
    
    

    return_dict = {'trait_trs': trait_trs,
                   'pval': (np.sum(trait_trs < ctrl_trs, axis=0) + 1) / (num_ctrl + 1),
                   'trait_ep': trait_ep,
                   'index': score_dataset.obs.index,
                   'gene_list': gene_list}
    
    if isinstance(dataset, str) and dataset.startswith('tms'):
        # Within tissue normalization
        within_tissue_trait_ep = np.zeros_like(trait_ep)
        for tissue in score_dataset.obs['tissue'].unique():
            tissue_index = np.where(score_dataset.obs['tissue'] == tissue)[0]
            within_tissue_trait_ep[tissue_index] = get_p_from_empi_null(trait_zsc[tissue_index], ctrl_zsc[:, tissue_index].flatten())
        return_dict['within_tissue_trait_ep'] = within_tissue_trait_ep
    
    def _calculate_norm_factor(v, robust_opt, axis=None):
        if robust_opt:
            med = np.median(v, axis=axis)
            q1, q3 = np.percentile(v, (25, 75), axis=axis)
            return med, q3 - q1
        else:
            mean = np.mean(v, axis=axis)
            std = np.std(v, axis=axis)
            return mean, std
    
    def _normalize_across_control(trait_zsc, ctrl_zsc, robust_opt):
        """
        robust_opt in [False, True]
        """
        assert robust_opt in [False, True]
        trait_mu, trait_sigma = _calculate_norm_factor(trait_zsc, robust_opt=robust_opt)
        norm_trait_zsc = (trait_zsc - trait_mu) / trait_sigma
        ctrl_mu, ctrl_sigma = _calculate_norm_factor(ctrl_zsc, robust_opt=robust_opt, axis=1)
        norm_ctrl_zsc = (ctrl_zsc - ctrl_mu.reshape(-1, 1)) / ctrl_sigma.reshape(-1, 1)
        return norm_trait_zsc, norm_ctrl_zsc
    
    for robust_opt, prefix in [[True, 'robust'],
                               [False, 'normal']]:
        norm_trait_zsc, norm_ctrl_zsc = _normalize_across_control(trait_zsc, ctrl_zsc, robust_opt=robust_opt)
        norm_trait_ep = get_p_from_empi_null(norm_trait_zsc, norm_ctrl_zsc.flatten())
        return_dict[f'{prefix}_trait_ep'] = norm_trait_ep
      
    return return_dict

def load_gene_score(path, gene_id_col, score_col, ascending=True, num_genes=None, hsapiens_mmusculus_mapping=None):
    """load gene score from file and sort 
    Args:
        path: the file path to the gene score, in csv format
        gene_id_col: the column corresponding the gene identifier
        score_col: the column corresponding to the gene score
        ascending: whether to rank genes by score in ascending order
        convert_mmusculus: whether to convert to mmusculus gene symbols, input a dataframe if needed
        
    Returns:
        gene_score (pd.DataFrame): a dataframe with sorted scores
    """
    df = pd.read_csv(path, usecols=[gene_id_col, score_col])
    df.columns =['GENE', 'SCORE']
    if hsapiens_mmusculus_mapping is not None:
        df = pd.merge(df, hsapiens_mmusculus_mapping, left_on='GENE', right_on='hsapiens')[['mmusculus', 'SCORE']].rename(columns={'mmusculus': 'GENE'})
    df = df.sort_values('SCORE', ascending=ascending)
    if num_genes is not None:
        df = df.iloc[0 : num_genes, :].reset_index(drop=True)
    return df

def zsc2pval(zsc):
    return 1 - sp.stats.norm.cdf(zsc)

def pval2zsc(pval):
    return -sp.stats.norm.ppf(pval).clip(min=-10,max=10)

def plot_simulation_score_hist(score_list, score_index, adata, stratify_by_tissue=True, num_bins=20, num_cols=6):
    
    num_sim = len(score_list)
    
    # plot the overall distribution
    hists = [np.histogram(score, bins=num_bins, range=(0, 1), density=True)[0] for score in score_list]
    xpos = np.arange(0, 1, 1 / num_bins) + 0.5 / num_bins
    for hist in hists:
        plt.plot(xpos, hist, 'k--', alpha=0.2)
    plt.errorbar(xpos,
             np.mean(hists, axis=0), 
             np.std(hists, axis=0) * 2 / np.sqrt(num_sim), 
             linestyle='None', marker='^', capsize=5)
    plt.ylim(0, 2)
    plt.axhline(y=1., linestyle='--', color='r')
    plt.xlabel('p-value bin')
    plt.ylabel('Density')
    plt.title('All')
    if stratify_by_tissue:
        tissue_df = adata[score_index].obs.tissue
        tissue_list = tissue_df.unique()
        plt.figure(figsize=[20, 2+3*len(tissue_list)/num_cols])
        for tissue_i, tissue in enumerate(tissue_list):
            tissue_index = (tissue_df == tissue)
            num_bins = 20
            plt.subplot(int(np.ceil(len(tissue_list) / num_cols)), num_cols, tissue_i + 1)
            hists = [np.histogram(score[tissue_index], bins=num_bins, range=(0, 1), density=True)[0] for score in score_list]

            xpos = np.arange(0, 1, 1 / num_bins) + 0.5 / num_bins
            for hist in hists:
                plt.plot(xpos, hist, 'k--', alpha=0.1)

            plt.errorbar(xpos, 
                         np.mean(hists, axis=0), 
                         np.std(hists, axis=0) * 2 / np.sqrt(num_sim), 
                         linestyle='None', marker='^', capsize=5)
            plt.ylim(0, 2)
            plt.axhline(y=1., linestyle='--', color='r')
            plt.xlabel('p-value bin')
            plt.ylabel('Density')
            plt.title(tissue)
        plt.tight_layout()
    plt.show()
    
def plot_simulation_pval_calibration(pval_dict):
    """
    score_dict: method -> #sims-length list of scores
    """
    plt.figure(figsize=(10, 5))
    cutoff_list = [2, 3, 4]
    
    for cutoff_i, cutoff in enumerate(cutoff_list):
        plt.subplot(1, len(cutoff_list), cutoff_i + 1)
        thres = 5 * 10 ** (-cutoff)
        
        sig_prop_dict = dict()
        for name in pval_dict:
            sig_prop_dict[name] = [np.mean(pval < thres) for pval in pval_dict[name]]
        xticks = np.arange(len(sig_prop_dict))
        plt.errorbar(xticks, [np.mean(sig_prop_dict[name]) for name in sig_prop_dict],
             yerr=[2 * np.std(sig_prop_dict[name]) / np.sqrt(len(sig_prop_dict[name])) for name in sig_prop_dict],
             color='k', fmt='.', capsize=6)
        plt.xticks(xticks, [name for name in pval_dict], rotation=90)
        plt.ylim(0, thres * 2)
        plt.axhline(thres, linestyle='--', color='r', )
        plt.title(f'Cutoff p<1e-{cutoff}')
        plt.ylabel(f'Prop. of cells\nwith p<1e-{cutoff}')
    plt.tight_layout()
    plt.show()
    
def plot_assoc_matrix(score_dict, score_index, adata, stratify_by):
    """
    plot tissue x traits association matrix 
    """


    def num2str(x):
        if x>1000:
            return '%0.1fk'%(x/1000)
        elif x>0:
            return '%d'%x
        else:
            return ''
    
    assert stratify_by in ['tissue', 'tissue_celltype', 'celltype']
    stratify_list = sorted(list(set(adata.obs[stratify_by])))
    
    # make df for plotting
    df_plot = pd.DataFrame(index=stratify_list, columns=list(score_dict.keys()), data=0)
    df_obs = adata.obs.copy()

    for trait in score_dict:

        pval = score_dict[trait]
        fdr = multipletests(pval, method='fdr_bh')[1]
        # TODO: check the concordance between score_index and adata
        # TODO: Deal with the len properly??
        temp_df = df_obs.loc[fdr<0.2].copy()
        temp_df = temp_df.groupby([stratify_by]).agg({stratify_by:len})
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
    if stratify_by == 'tissue-celltype':
        plt.figure(figsize=[0.4*df_plot.shape[1]+3, 0.25*df_plot.shape[0]+2+1])
    else:
        plt.figure(figsize=[0.4*df_plot.shape[1]+3, 0.25*df_plot.shape[0]+2])
    sns.heatmap(df_plot, annot=mat_annot, fmt='s', cbar=False)
    plt.xticks(np.arange(df_plot.shape[1])+0.5, df_plot.columns, rotation=45, ha='right')
    plt.show()
    
    
def plot_qqplot(pval_dict, num_cols=6):
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
        util.qqplot(x=normal_x, y=trait_logpval, quantiles=2000, s=3, alpha=0.5)
        plt.axline((1, 1), slope=1, linestyle='--', color='k', alpha=0.5)

        plt.title('%s\nFDR<0.2: %d cells'%(trait,(fdr_dict[trait]<0.2).sum()))
        plt.ylabel('%s\nObserved[-log10(P)]'%trait)
        plt.xlabel('Expected[-log10(P)]')
    plt.tight_layout()
    plt.show()

def plot_pval_umap(pval_dict, score_index, umap_adata_dict, umap_color=['cell_ontology_class'], num_cols=5):
    """
    Overlay p-value on UMAP
    """
    for tissue in umap_adata_dict:
        
        umap_adata = umap_adata_dict[tissue].copy() 
        fig, ax = plt.subplots(figsize=(6, 6))
        sc.pl.umap(umap_adata, color=umap_color, size=20, ax=ax)
        df_plot = pd.DataFrame(index=umap_adata.obs.index)
        df_plot['UMAP1'] = umap_adata.obsm['X_umap'][:,0]
        df_plot['UMAP2'] = umap_adata.obsm['X_umap'][:,1]

        # Trait TRS plot
        plt.figure(figsize=[15, 2+3*len(pval_dict)/num_cols])
        
        for trait_i, trait in enumerate(pval_dict.keys()):
            plt.subplot(int(np.ceil(len(pval_dict) / num_cols)), num_cols, trait_i + 1)
            temp_df = pd.DataFrame(index=score_index)
            temp_df[trait] = -sp.stats.norm.ppf(pval_dict[trait]).clip(min=-10,max=10)
            df_plot = df_plot.join(temp_df[trait])

            max_ = np.quantile(np.absolute(df_plot[trait].values), 0.99)
            min_ = np.quantile(np.absolute(df_plot[trait].values), 0.01)
            plt.scatter(df_plot['UMAP1'], df_plot['UMAP2'], c=df_plot[trait],
                        cmap='RdBu_r', vmax=max_, vmin=-max_, s=4)
            plt.colorbar()
            plt.clim(-4,4)
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP2')
            plt.title('%s'%trait)

        plt.tight_layout()
        plt.show()
        

def plot_score_correlation(score_dict, score_index, adata, stratify_by):
    """
    Plot the correlation matrix between pairs of traits
    stratify_by: a list of configurations: e.g. stratify_by = {'all': True, tissue': ['Liver'], 'tissue_celltype', 'Heart.B cell']]
    """
    df = pd.DataFrame()
    for trait in score_dict:
        df[trait] = score_dict[trait]
    
    figsize = len(score_dict) / 3 * 2
    
    for group_name in stratify_by:
        if group_name == 'all':
            
            df_corr = df.corr()
            plt.figure(figsize=[figsize, figsize])
            sns.heatmap(df_corr, annot=df_corr, fmt='0.2f', xticklabels=False, center=0.0, vmax=1.0, vmin=-1.0)
            plt.title('TRS z-score correlation across all cells ')
            plt.show()
        else:
            for group in stratify_by[group_name]:
                index = np.where(adata.obs[group_name] == group)[0]
                group_df = df.iloc[index, :]
                df_corr = group_df.corr()

                plt.figure(figsize=[figsize, figsize])
                sns.heatmap(df_corr, annot=df_corr, fmt='0.2f', xticklabels=False, center=0.0, vmax=1.0, vmin=-1.0)
                plt.title(f'TRS z-score correlation across {group} cells ')
                plt.show()
    
#     if 'tissue' in stratify_by:
#         for tissue in stratify_by['tissue']:
#             index = np.where(adata.obs['tissue'] == tissue)[0]
#             subset_df = df.iloc[index, :]
#             df_corr = subset_df.corr()
            
#             plt.figure(figsize=[figsize, figsize])
#             sns.heatmap(df_corr, annot=df_corr, fmt='0.2f', xticklabels=False, center=0.0, vmax=1.0, vmin=-1.0)
#             plt.title(f'TRS z-score correlation across {tissue} cells ')
#             plt.show()
            
#     if 'tissue_celltype' in stratify_by:
#         for tc in stratify_by['tissue_celltype']:
#             index = np.where(adata.obs['tissue_celltype'] == tc)[0]
#             subset_df = df.iloc[index, :]
#             df_corr = subset_df.corr()
            
#             plt.figure(figsize=[figsize, figsize])
#             sns.heatmap(df_corr, annot=df_corr, fmt='0.2f', xticklabels=False, center=0.0, vmax=1.0, vmin=-1.0)
#             plt.title(f'TRS z-score correlation across {tc} cells ')
#             plt.show()
            