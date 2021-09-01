import scanpy as sc
from anndata import read_h5ad
import pandas as pd
import numpy as np
import scipy as sp
import os
from os.path import join
import time
import argparse
from statsmodels.stats.multitest import multipletests

# inhouse tools
import scdrs.util as util
import scdrs.data_loader as dl
import scdrs.method as md


"""
# Fixit
- Warning for compute_score: Trying to set attribute `.X` of view, copying.


# Todo
- Support other data formats in addition to h5ad
- Implement a memory efficient version
- "gene_weight" argument needs to be tested 

# Finished
- Add --n_ctrl (default value 500) 
- Add --cov_file option to regress out covariates stored in COV_FILE before feeding into the score function 
- Add --ctrl_match_opt='mean_var': use mean- and var- matched control genes 
- Change name from scTRS to scdrs (072721)

"""

VERSION='0.0.1'

def convert_species_name(species):
    if species in ['Mouse', 'mouse', 'Mus_musculus', 'mus_musculus', 'mmusculus']:
        return 'mmusculus'
    if species in ['Human', 'human', 'Homo_sapiens', 'homo_sapiens', 'hsapiens']:
        return 'hsapiens'
    raise ValueError('# compute_score: species name %s not supported'%species)


def main(args):
    sys_start_time = time.time()
    
    MASTHEAD = "******************************************************************************\n"
    MASTHEAD += "* Single-cell disease relevance score (scDRS)\n"
    MASTHEAD += "* Version %s\n"%VERSION
    MASTHEAD += "* Martin Jinye Zhang and Kangcheng Hou\n"
    MASTHEAD += "* HSPH / Broad Institute / UCLA\n"
    # MASTHEAD += "* GNU General Public License v3\n"
    MASTHEAD += "******************************************************************************\n"
        
    ###########################################################################################    
    ######                                    Parse Options                              ######
    ###########################################################################################    
    H5AD_FILE=args.h5ad_file
    H5AD_SPECIES=args.h5ad_species
    COV_FILE=args.cov_file
    GS_FILE=args.gs_file
    GS_SPECIES=args.gs_species
    CTRL_MATCH_OPT=args.ctrl_match_opt
    WEIGHT_OPT=args.weight_opt
    FLAG_FILTER=args.flag_filter=='True'
    FLAG_RAW_COUNT=args.flag_raw_count=='True'
    N_CTRL=int(args.n_ctrl)
    FLAG_RETURN_CTRL_RAW_SCORE=args.flag_return_ctrl_raw_score=='True'
    FLAG_RETURN_CTRL_NORM_SCORE=args.flag_return_ctrl_norm_score=='True'
    OUT_FOLDER=args.out_folder
    
    if H5AD_SPECIES!=GS_SPECIES:
        H5AD_SPECIES=convert_species_name(H5AD_SPECIES)
        GS_SPECIES=convert_species_name(GS_SPECIES)
        
    header = MASTHEAD
    header += 'Call: ./compute_score.py \\\n'
    header += '--h5ad_file %s\\\n'%H5AD_FILE
    header += '--h5ad_species %s\\\n'%H5AD_SPECIES
    header += '--cov_file %s\\\n'%COV_FILE
    header += '--gs_file %s\\\n'%GS_FILE
    header += '--gs_species %s\\\n'%GS_SPECIES
    header += '--ctrl_match_opt %s\\\n'%CTRL_MATCH_OPT
    header += '--weight_opt %s\\\n'%WEIGHT_OPT
    header += '--flag_filter %s\\\n'%FLAG_FILTER
    header += '--flag_raw_count %s\\\n'%FLAG_RAW_COUNT
    header += '--n_ctrl %d\\\n'%N_CTRL
    header += '--flag_return_ctrl_raw_score %s\\\n'%FLAG_RETURN_CTRL_RAW_SCORE
    header += '--flag_return_ctrl_norm_score %s\\\n'%FLAG_RETURN_CTRL_NORM_SCORE
    header += '--out_folder %s\n'%OUT_FOLDER
    print(header)
    
    # Check options 
    if H5AD_SPECIES!=GS_SPECIES:
        if H5AD_SPECIES not in ['mmusculus', 'hsapiens']:
            raise ValueError('H5AD_SPECIES needs to be one of [mmusculus, hsapiens] '
                             'unless H5AD_SPECIES==GS_SPECIES')
        if GS_SPECIES not in ['mmusculus', 'hsapiens']:
            raise ValueError('GS_SPECIES needs to be one of [mmusculus, hsapiens] '
                             'unless H5AD_SPECIES==GS_SPECIES')
    if CTRL_MATCH_OPT not in ['mean', 'mean_var']:
        raise ValueError('CTRL_MATCH_OPT needs to be one of [mean, mean_var]')
    if WEIGHT_OPT not in ['uniform', 'vs', 'inv_std', 'adapt', 'od']:
        raise ValueError('WEIGHT_OPT needs to be one of [uniform, vs, inv_std, adapt, od]')
    
    ###########################################################################################    
    ######                                     Load data                                 ######
    ###########################################################################################
    print('Load data:')
    
    # Load .h5ad file 
    adata = read_h5ad(H5AD_FILE)
    if FLAG_FILTER:
        sc.pp.filter_cells(adata, min_genes=250)
        sc.pp.filter_genes(adata, min_cells=50)
    if FLAG_RAW_COUNT:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
    print('--h5ad_file loaded: n_cell=%d, n_gene=%d (sys_time=%0.1fs)'
          %(adata.shape[0], adata.shape[1], time.time()-sys_start_time))
    
    # adata = adata[0:500,:].copy()
    
    # Load .cov file and regress out covariates 
    if COV_FILE is not None:
        df_cov = pd.read_csv(COV_FILE, sep='\t', index_col=0)
        cov_list = list(df_cov.columns)
        adata.obs.drop([x for x in cov_list if x in adata.obs.columns], axis=1, inplace=True)
        adata.obs = adata.obs.join(df_cov)
        adata.obs.fillna(adata.obs[cov_list].mean(), inplace=True)
        print('--cov_file loaded: covariates=[%s] (sys_time=%0.1fs)'
              %(', '.join(cov_list), time.time()-sys_start_time))
    
        adata.var['mean'] = adata.X.mean(axis=0).T
        if sp.sparse.issparse(adata.X):
            adata.X = adata.X.toarray()
        adata.X -= adata.var['mean'].values
        adata.X = md.reg_out(adata.X, adata.obs[cov_list].values)
        adata.X += adata.var['mean']
        print('Regress out covariates from --h5ad_file (sys_time=%0.1fs)'
              %(time.time()-sys_start_time))
    
    # Load .gs file 
    df_gs = pd.read_csv(GS_FILE, sep='\t')
    df_gs.index = df_gs['TRAIT']
    print('--gs_file loaded: n_geneset=%d (sys_time=%0.1fs)'
          %(df_gs.shape[0], time.time()-sys_start_time))
    
    # Convert df_gs genes to H5AD_SPECIES genes
    if H5AD_SPECIES!=GS_SPECIES:
        # Load homolog file 
        df_hom = pd.read_csv('/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/'
                             'gene_annotation/mouse_human_homologs.txt', sep='\t')
        if (GS_SPECIES=='hsapiens') & (H5AD_SPECIES=='mmusculus'):
            dic_map = {x:y for x,y in zip(df_hom['HUMAN_GENE_SYM'], df_hom['MOUSE_GENE_SYM'])}
        elif (GS_SPECIES=='mmusculus') & (H5AD_SPECIES=='hsapiens'):
            dic_map = {x:y for x,y in zip(df_hom['MOUSE_GENE_SYM'], df_hom['HUMAN_GENE_SYM'])}
        else:
            raise ValueError('# compute_score: gene conversion from %s to %s is not supported'
                             %(GS_SPECIES, H5AD_SPECIES))
    
        for trait in df_gs.index:
            gs_gene_list = df_gs.loc[trait, 'GENESET'].split(',')
            h5ad_gene_list = [dic_map[x] for x in set(gs_gene_list)&set(dic_map.keys())]
            df_gs.loc[trait, 'GENESET'] = ','.join(h5ad_gene_list)
        print('--gs_file converted from %s to %s genes (sys_time=%0.1fs)'
              %(GS_SPECIES, H5AD_SPECIES, time.time()-sys_start_time))
    print('')
    ###########################################################################################    
    ######                                  Computation                                  ######
    ###########################################################################################
    
    # Compute statistics, including the 20*20 mean-var bins
    md.compute_stats(adata)
    print('Compute basic statistics (sys_time=%0.1fs)'%(time.time()-sys_start_time))
    
    # Compute score 
    for trait in df_gs.index:
        gene_list = df_gs.loc[trait,'GENESET'].split(',')
        gene_list = sorted(set(gene_list) & set(adata.var_names))
        if len(gene_list)<10:
            print('Gene set %s: skipped due to small size (n_gene=%d, sys_time=%0.1fs)'
                  %(trait, len(gene_list), time.time()-sys_start_time))
            continue
            
        df_res = md.score_cell(adata, gene_list, ctrl_match_key=CTRL_MATCH_OPT, n_ctrl=N_CTRL,
                               weight_opt=WEIGHT_OPT,
                               return_ctrl_raw_score=FLAG_RETURN_CTRL_RAW_SCORE, 
                               return_ctrl_norm_score=FLAG_RETURN_CTRL_NORM_SCORE, verbose=False)
        df_res.iloc[:,0:6].to_csv(join(OUT_FOLDER, '%s.score.gz'%trait), sep='\t',
                                  index=True, compression='gzip')
        if FLAG_RETURN_CTRL_RAW_SCORE|FLAG_RETURN_CTRL_NORM_SCORE:
            df_res.to_csv(join(OUT_FOLDER, '%s.full_score.gz'%trait), sep='\t',
                          index=True, compression='gzip')
        v_fdr = multipletests(df_res['pval'].values, method='fdr_bh')[1]
        n_rej_01 = (v_fdr<0.1).sum()
        n_rej_02 = (v_fdr<0.2).sum()
        print('Gene set %s (n_gene=%d): %d/%d FDR<0.1 cells, %d/%d FDR<0.2 cells (sys_time=%0.1fs)'
              %(trait, len(gene_list), n_rej_01, df_res.shape[0], n_rej_02, df_res.shape[0],
                time.time()-sys_start_time))
            
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='compute score')
    
    parser.add_argument('--h5ad_file', type=str, required=True)
    parser.add_argument('--h5ad_species', type=str, required=True, help='one of [hsapiens, mmusculus]')
    parser.add_argument('--cov_file', type=str, required=False, default=None)    
    parser.add_argument('--gs_file', type=str, required=True)
    parser.add_argument('--gs_species', type=str, required=True, help='one of [hsapiens, mmusculus]')
    parser.add_argument('--ctrl_match_opt', type=str, required=False, default='mean_var')
    parser.add_argument('--weight_opt', type=str, required=False, default='vs')
    parser.add_argument('--flag_filter', type=str, required=False, default='True')
    parser.add_argument('--flag_raw_count', type=str, required=False, default='True')
    parser.add_argument('--n_ctrl', type=int, required=False, default=1000)
    parser.add_argument('--flag_return_ctrl_raw_score', type=str, required=False, default='False')
    parser.add_argument('--flag_return_ctrl_norm_score', type=str, required=False, default='False')
    parser.add_argument('--out_folder', type=str, required=True)
    
    args = parser.parse_args()
    
    main(args)