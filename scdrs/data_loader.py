import scanpy as sc
import numpy as np
import time
from anndata import read_h5ad
from anndata import AnnData
import scdrs.method as md

def load_tms_processed(file_path, data_name='facs', tissue='all'):
    
    """load processed tms data
    
    Args:
        file_path (str): scDRS_data path
        data_name (str): one of 'facs', 'droplet'
        tissue (str): 'all' or one of the facs or droplet tissues
        

    Returns:
        adata (AnnData): Combined data for FACS and droplet
    """
    
    dic_tissue_list = {'facs': ['Aorta', 'BAT', 'Bladder', 'Brain_Myeloid', 'Brain_Non-Myeloid',
                                'Diaphragm', 'GAT', 'Heart', 'Kidney', 'Large_Intestine',
                                'Limb_Muscle', 'Liver', 'Lung', 'MAT', 'Mammary_Gland', 'Marrow', 
                                'Pancreas', 'SCAT', 'Skin', 'Spleen', 'Thymus', 'Tongue', 'Trachea'],
                       'droplet': ['Bladder', 'Fat', 'Heart_and_Aorta', 'Kidney', 'Large_Intestine',
                                   'Limb_Muscle', 'Liver', 'Lung', 'Mammary_Gland', 'Marrow',
                                   'Pancreas', 'Skin', 'Spleen', 'Thymus', 'Tongue', 'Trachea']}
    
    if data_name not in ['facs', 'droplet']:
        raise ValueError("# load_tms_processed: data_name must be one of [facs, droplet]")
    
    if (tissue!='all') & (tissue not in dic_tissue_list[data_name]):
        raise ValueError("# load_tms_processed: tissue for %s must be 'all' or one of [%s]"
                         %(data_name, ', '.join(dic_tissue_list[data_name])))
    
    if tissue == 'all':
        tissue_list = dic_tissue_list[data_name]
    else:
        tissue_list = [tissue]
        
    print("# load_tms_processed: load %s data, tissue=[%s]"
          %(data_name, ', '.join(tissue_list)))
            
    # Load filtered data
    dic_data = {}
    for x in tissue_list:
        dic_data[x] = read_h5ad('%s/tabula_muris_senis/tabula-muris-senis-%s-processed-official-annotations-%s.h5ad'
                                %(file_path, data_name, x))
    return dic_data


def load_tms_ct(file_path, data_name='facs', flag_size_factor=True, 
                total_ct_per_cell=1e4, flag_log1p=True,
                flag_scale=False, flag_filter_genes=False):
    
    """load normalized tms ct data
    1. Load raw data ('facs' or 'droplet')
    2. Size factor normalization to counts total_ct_per_cell
    3. log(x+1) transform
    4. [Optional] scale each gene to have unit scale

    Args:
        file_path (str): file path. Should contain both 
            "tabula-muris-senis-droplet-official-raw-obj.h5ad" and 
            "tabula-muris-senis-facs-official-raw-obj.h5ad"

    Returns:
        adata (AnnData): Combined data for FACS and droplet
    """
    
    if data_name not in ['facs', 'droplet']:
        raise ValueError('# load_tms_ct: %s need to be one of [facs, droplet]'%(data_name))
        
    if data_name=='facs':
        file_name = 'tabula_muris_senis/tabula-muris-senis-facs-official-raw-obj.h5ad'
    elif data_name=='droplet':
        file_name = 'tabula_muris_senis/tabula-muris-senis-droplet-official-raw-obj.h5ad'
    else:
        return None
            
    # Load filtered data
    adata = read_h5ad(file_path+'/'+file_name)
    
    if flag_filter_genes == True:
        sc.pp.filter_genes(adata, min_cells=3)
        
    # Size factor normalization
    if flag_size_factor == True:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        
    # log(x+1) transform
    if flag_log1p == True:
        sc.pp.log1p(adata)
    
    # if flag_scale == True:
    #     sc.pp.scale(adata, max_value=10, zero_center=False)
    
    if 'facs' in data_name:
        ind_select = adata.obs['age'].isin(['3m', '18m', '24m'])
        adata = adata[ind_select,]
        
    # Add extra annotations
    adata.obs['age_num'] = [int(x.replace('m','')) for x in adata.obs['age']]
    adata.obs['tissue_celltype'] = ['%s.%s'%(x,y) for x,y in zip(adata.obs['tissue'],
                                                                 adata.obs['cell_ontology_class'])]
    return adata


def load_ts_facs(file_path, opt='raw',
                 flag_size_factor=True,
                 total_ct_per_cell=1e4,
                 flag_log1p=True):
    
    if opt=='raw':
        adata = read_h5ad(file_path+'/single_cell_data/tabula_sapiens/obj_smartseq2_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        pass
        
    return adata


def load_canogamez(file_path, opt='raw',
                   flag_size_factor=True,
                   total_ct_per_cell=1e4,
                   flag_log1p=True):
    
    if opt=='raw':
        adata = read_h5ad(file_path+'/single_cell_data/tcell_canogamez_nc_2020/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        adata = read_h5ad(file_path+'/single_cell_data/tcell_canogamez_nc_2020/obj_processed.h5ad')
        
    return adata


def load_nathan(file_path, opt='raw',
                flag_size_factor=True,
                total_ct_per_cell=1e4,
                flag_log1p=True):
    
    if opt=='raw':
        adata = read_h5ad(file_path+'/single_cell_data/nathan_ni_2021/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        adata = read_h5ad(file_path+'/single_cell_data/nathan_ni_2021/obj_processed.h5ad')
        
    return adata


def load_aizarani(file_path, opt='raw',
                  flag_size_factor=True,
                  total_ct_per_cell=1e4,
                  flag_log1p=True):
    
    if opt=='raw':
        adata = read_h5ad(file_path+'/single_cell_data/liver_atlas/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/liver_atlas'
        adata = read_h5ad(file_path+'/single_cell_data/liver_atlas/obj_processed.h5ad')
        
    return adata


def load_halpern(file_path, opt='raw',
                 flag_size_factor=True,
                 total_ct_per_cell=1e4,
                 flag_log1p=True):
    
    if opt=='raw':
        adata = read_h5ad(file_path+'/single_cell_data/mouse_liver_halpern_nature_2017/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        adata = read_h5ad(file_path+'/single_cell_data/mouse_liver_halpern_nature_2017/obj_processed.h5ad')
        
    return adata


def load_richter(file_path, opt='raw',
                 flag_size_factor=True,
                 total_ct_per_cell=1e4,
                 flag_log1p=True):
    
    if opt=='raw':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/richter_biorxiv_2020'
        adata = read_h5ad(file_path+'/single_cell_data/richter_biorxiv_2020/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        adata = read_h5ad(file_path+'/single_cell_data/richter_biorxiv_2020/obj_processed.h5ad')
        
    return adata


def load_rao(file_path, opt='raw',
             flag_size_factor=True,
             total_ct_per_cell=1e4,
             flag_log1p=True):
    
    if opt=='raw':
        adata = read_h5ad(file_path+'/single_cell_data/rao_wang_2021/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        adata = read_h5ad(file_path+'/single_cell_data/rao_wang_2021/obj_processed.h5ad')
        
    return adata


def load_asp_st_raw_data(opt='raw',
                         flag_size_factor=True,
                         total_ct_per_cell=1e4,
                         flag_log1p=True):
    
    if opt=='raw':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/heart_asp_cell_2019'
        adata = read_h5ad(DATA_PATH+'/obj_raw_st.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/heart_asp_cell_2019'
        adata = read_h5ad(DATA_PATH+'/obj_processed_st.h5ad')
        
    return adata


def load_asp_sc_raw_data(opt='raw',
                         flag_size_factor=True,
                         total_ct_per_cell=1e4,
                         flag_log1p=True):
    
    if opt=='raw':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/heart_asp_cell_2019'
        adata = read_h5ad(DATA_PATH+'/obj_raw_sc.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/heart_asp_cell_2019'
        adata = read_h5ad(DATA_PATH+'/obj_processed_sc.h5ad')
        
    return adata


def load_xin_raw_data(opt='raw',
                      flag_size_factor=True,
                      total_ct_per_cell=1e4,
                      flag_log1p=True):
    
    if opt=='raw':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/xin_diabetes_2018'
        adata = read_h5ad(DATA_PATH+'/obj_raw.h5ad')
        
        if flag_size_factor == True:
            sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
        if flag_log1p == True:
            sc.pp.log1p(adata)
        
    if opt=='processed':
        DATA_PATH='/n/holystore01/LABS/price_lab/Users/mjzhang/scDRS_data/single_cell_data/xin_diabetes_2018'
        adata = read_h5ad(DATA_PATH+'/obj_processed.h5ad')
        
    return adata
