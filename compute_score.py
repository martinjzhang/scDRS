import scanpy as sc
from anndata import read_h5ad
import pandas as pd
import numpy as np
import scipy as sp
import os
from os.path import join
import time
import argparse

# inhouse tools
import scTRS.util as util
import scTRS.data_loader as dl
import scTRS.method as md


"""
# Fixit
- 

# Todo
- Support other data formats in addition to h5ad

# Finished
- Add --n_ctrl (default value 500) 
- Add --cov_file option to regress out covariates stored in COV_FILE before feeding into the score function 
- Add --ctrl_match_opt='mean_var': use mean- and var- matched control genes 

"""

def convert_species_name(species):
    if species in ['Mouse', 'mouse', 'Mus_musculus', 'mus_musculus', 'mmusculus']:
        return 'mmusculus'
    if species in ['Human', 'human', 'Homo_sapiens', 'homo_sapiens', 'hsapiens']:
        return 'hsapiens'
    raise ValueError('# compute_score: species name %s not supported'%species)


def main(args):
    sys_start_time = time.time()
        
    ###########################################################################################    
    ######                                    Parse Options                              ######
    ###########################################################################################    
    H5AD_FILE=args.h5ad_file
    H5AD_SPECIES=args.h5ad_species
    COV_FILE=args.cov_file
    GS_FILE=args.gs_file
    GS_SPECIES=args.gs_species
    CTRL_MATCH_OPT=args.ctrl_match_opt
    FLAG_FILTER=args.flag_filter=='True'
    FLAG_RAW_COUNT=args.flag_raw_count=='True'
    N_CTRL=args.n_ctrl
    FLAG_RETURN_CTRL_RAW_SCORE=args.flag_return_ctrl_raw_score=='True'
    FLAG_RETURN_CTRL_NORM_SCORE=args.flag_return_ctrl_norm_score=='True'
    OUT_FOLDER=args.out_folder
    
    if H5AD_SPECIES!=GS_SPECIES:
        H5AD_SPECIES=convert_species_name(H5AD_SPECIES)
        GS_SPECIES=convert_species_name(GS_SPECIES)
    
    print('# compute_score: H5AD_FILE: ', H5AD_FILE)
    print('# compute_score: H5AD_SPECIES: ', H5AD_SPECIES)
    print('# compute_score: COV_FILE: ', COV_FILE)
    print('# compute_score: GS_FILE: ', GS_FILE)
    print('# compute_score: GS_SPECIES: ', GS_SPECIES)
    print('# compute_score: CTRL_MATCH_OPT: ', CTRL_MATCH_OPT)
    print('# compute_score: FLAG_FILTER: ', FLAG_FILTER)
    print('# compute_score: FLAG_RAW_COUNT: ', FLAG_RAW_COUNT)
    print('# compute_score: N_CTRL: ', N_CTRL)
    print('# compute_score: FLAG_RETURN_CTRL_RAW_SCORE: ', FLAG_RETURN_CTRL_RAW_SCORE)
    print('# compute_score: FLAG_RETURN_CTRL_NORM_SCORE: ', FLAG_RETURN_CTRL_NORM_SCORE)
    print('# compute_score: OUT_FOLDER: ', OUT_FOLDER)
    
    if CTRL_MATCH_OPT not in ['mean', 'mean_var']:
        raise ValueError('# compute_score: CTRL_MATCH_OPT needs to be one of [mean, mean_var]')
            
    ###########################################################################################    
    ######                                   Data Loading                                ######
    ###########################################################################################
    
    # Load .h5ad file 
    adata = read_h5ad(H5AD_FILE)
    if FLAG_FILTER:
        sc.pp.filter_cells(adata, min_genes=250)
        sc.pp.filter_genes(adata, min_cells=50)
    if FLAG_RAW_COUNT:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
    print('# compute_score: H5AD_FILE loaded: n_cell=%d, n_gene=%d'%adata.shape)
    
    # adata = adata[0:500,:].copy()
    
    # Load .cov file and regress out covariates 
    if COV_FILE is not None:
        df_cov = pd.read_csv(COV_FILE, sep='\t', index_col=0)
        cov_list = list(df_cov.columns)
        adata.obs.drop([x for x in cov_list if x in adata.obs.columns], axis=1, inplace=True)
        adata.obs = adata.obs.join(df_cov)
        adata.obs.fillna(adata.obs[cov_list].mean(), inplace=True)
        print('# compute_score: COV_FILE loaded, cov_list=%s'%cov_list)
    
        adata.var['mean'] = adata.X.mean(axis=0).T
        if sp.sparse.issparse(adata.X):
            adata.X = adata.X.toarray()
        adata.X -= adata.var['mean'].values
        adata.X = md.reg_out(adata.X, adata.obs[cov_list].values)
        adata.X += adata.var['mean']
        print('# compute_score: covariates regressed out from adata.X')
    
    # Load .gs file 
    df_gs = pd.read_csv(GS_FILE, sep='\t')
    df_gs.index = df_gs['TRAIT']
    print('# compute_score: GS_FILE loaded: ', df_gs.shape)
    
    # Convert df_gs genes to H5AD_SPECIES genes
    if H5AD_SPECIES!=GS_SPECIES:
        # Load homolog file 
        df_hom = pd.read_csv('/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/gene_annotation/'
                             'mouse_human_homologs.txt', sep='\t')
        if (GS_SPECIES=='hsapiens') & (H5AD_SPECIES=='mmusculus'):
            dic_map = {x:y for x,y in zip(df_hom['HUMAN_GENE_SYM'], df_hom['MOUSE_GENE_SYM'])}
        elif (GS_SPECIES=='mmusculus') & (H5AD_SPECIES=='hsapiens'):
            dic_map = {x:y for x,y in zip(df_hom['MOUSE_GENE_SYM'], df_hom['HUMAN_GENE_SYM'])}
        else:
            raise ValueError('# compute_score: gene conversion from %s to %s is not supported'%(GS_SPECIES, H5AD_SPECIES))
    
        for trait in df_gs.index:
            gs_gene_list = df_gs.loc[trait, 'GENESET'].split(',')
            h5ad_gene_list = [dic_map[x] for x in set(gs_gene_list)&set(dic_map.keys())]
            df_gs.loc[trait, 'GENESET'] = ','.join(h5ad_gene_list)
        print('# compute_score: GS_FILE converted from %s to %s genes'%(GS_SPECIES, H5AD_SPECIES))
    print('# compute_score: sys_time=%0.1fs'%(time.time()-sys_start_time))
    
    ###########################################################################################    
    ######                                  Computation                                  ######
    ###########################################################################################
    
    # Divide genes into mean-var bins
    md.compute_stats(adata)
    if CTRL_MATCH_OPT=='mean_var':
        v_mean_bin = pd.qcut(adata.var['mean'], 20, labels=False, duplicates='drop')
        adata.var['mean_var'] = ''
        for bin_ in set(v_mean_bin):
            ind_select = (v_mean_bin==bin_)
            v_var_bin = pd.qcut(adata.var.loc[ind_select, 'var'], 20, labels=False, duplicates='drop')
            adata.var.loc[ind_select, 'mean_var'] = ['%s.%s'%(x,y) for x,y in zip(v_mean_bin[ind_select],v_var_bin)]
    
    # Compute score 
    for trait in df_gs.index:
        gene_list = df_gs.loc[trait,'GENESET'].split(',')
        gene_list = sorted(set(gene_list) & set(adata.var_names))
        if len(gene_list)<10:
            print('# %s skipped due to small size (n_gene=%d)'%(trait, len(gene_list)))
            continue
            
        df_res = md.score_cell(adata, gene_list, ctrl_match_key=CTRL_MATCH_OPT, n_ctrl=N_CTRL, 
                               return_ctrl_raw_score=FLAG_RETURN_CTRL_RAW_SCORE, 
                               return_ctrl_norm_score=FLAG_RETURN_CTRL_NORM_SCORE, verbose=False)
        df_res.iloc[:,0:6].to_csv(join(OUT_FOLDER, '%s.score.gz'%trait), sep='\t', index=True, compression='gzip')
        if FLAG_RETURN_CTRL_RAW_SCORE|FLAG_RETURN_CTRL_NORM_SCORE:
            df_res.to_csv(join(OUT_FOLDER, '%s.full_score.gz'%trait), sep='\t', index=True, compression='gzip')
        print('# compute_score: score computed for %s (%d genes), sys_time=%0.1fs'
              %(trait, len(gene_list), time.time()-sys_start_time))
            
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='compute score')
    
    parser.add_argument('--h5ad_file', type=str, required=True)
    parser.add_argument('--h5ad_species', type=str, required=True, help='one of [hsapiens, mmusculus]')
    parser.add_argument('--cov_file', type=str, required=False, default=None)    
    parser.add_argument('--gs_file', type=str, required=True)
    parser.add_argument('--gs_species', type=str, required=True, help='one of [hsapiens, mmusculus]')
    parser.add_argument('--ctrl_match_opt', type=str, required=False, default='mean_var')
    parser.add_argument('--flag_filter', type=str, required=False, default='True')
    parser.add_argument('--flag_raw_count', type=str, required=False, default='True')
    parser.add_argument('--n_ctrl', type=int, required=False, default=500)
    parser.add_argument('--flag_return_ctrl_raw_score', type=str, required=False, default='False')
    parser.add_argument('--flag_return_ctrl_norm_score', type=str, required=False, default='False')
    parser.add_argument('--out_folder', type=str, required=True)
    
    args = parser.parse_args()
    main(args)