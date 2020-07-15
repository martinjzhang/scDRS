import scanpy as sc
import numpy as np
import time
from anndata import read_h5ad
from anndata import AnnData

def load_tms_raw(file_path, data_name='facs', 
                 flag_size_factor=True, 
                 total_ct_per_cell=1e4, 
                 flag_log1p=True):
    
    """load normalized data
    1. Load filtered data for both FACS and droplet
    2. Size factor normalization to counts per 1 million (total_ct_per_cell)
    3. log(x+1) transform
    4. Combine the data 

    Args:
        file_path (str): file path. Should contain both 
            "tabula-muris-senis-droplet-official-raw-obj.h5ad" and 
            "tabula-muris-senis-facs-official-raw-obj.h5ad"

    Returns:
        adata (AnnData): Combined data for FACS and droplet
    """
    
    if data_name not in ['facs', 'droplet']:
        raise ValueError('# load_tms: %s need to be one of [facs, droplet]')
        
    if data_name=='facs':
        file_name = 'tabula-muris-senis-facs-official-raw-obj.h5ad'
    elif data_name=='droplet':
        file_name = 'tabula-muris-senis-droplet-official-raw-obj.h5ad'
    else:
        return None
            
    # Load filtered data
    adata = read_h5ad(f'{file_path}/{file_name}')
    # Size factor normalization
    if flag_size_factor == True:
        # sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
    # log(x+1) transform
    if flag_log1p == True:
        sc.pp.log1p(adata)
    # Combine the data 
    if 'facs' in data_name:
        ind_select = adata.obs['age'].isin(['3m', '18m', '24m'])
        adata = adata[ind_select,]
    # Add extra annotations
    adata.obs['age_num'] = [int(x.replace('m','')) for x in adata.obs['age']]
    return adata

def load_tms_processed(file_path, data_name='facs'):
    
    """load normalized data
    1. Load filtered data for both FACS and droplet
    2. Size factor normalization to counts per 1 million (total_ct_per_cell)
    3. log(x+1) transform
    4. Combine the data 

    Args:
        file_path (str): file path. Should contain both 
            "tabula-muris-senis-droplet-official-raw-obj.h5ad" and 
            "tabula-muris-senis-facs-official-raw-obj.h5ad"

    Returns:
        adata (AnnData): Combined data for FACS and droplet
    """
    
    # if data_name not in ['facs', 'droplet']:
    #     raise ValueError('# load_tms: %s need to be one of [facs, droplet]')
    #     
    # if data_name=='facs':
    #     file_name = 'tabula-muris-senis-facs-official-raw-obj.h5ad'
    # elif data_name=='droplet':
    #     file_name = 'tabula-muris-senis-droplet-official-raw-obj.h5ad'
    # else:
    #     return None
            
    # Load filtered data
    print('%s/tms_processed_data/%s.h5ad'%(file_path, data_name))
    adata = read_h5ad('%s/tms_processed_data/%s.h5ad'%(file_path, data_name))
    return adata