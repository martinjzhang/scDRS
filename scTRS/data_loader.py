import scanpy as sc
import numpy as np
import time
from anndata import read_h5ad
from anndata import AnnData

def load_tms_processed(file_path, data_name='facs', tissue='all'):
    
    """load processed tms data
    
    Args:
        file_path (str): scTRS_data path
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
          %(data_name, ', '.join(dic_tissue_list[data_name])))
            
    # Load filtered data
    dic_data = {}
    for x in tissue_list:
        dic_data[x] = read_h5ad('%s/tabula_muris_senis/tabula-muris-senis-%s-processed-official-annotations-%s.h5ad'
                                %(file_path, data_name, x))
    return dic_data



# def load_tms_raw(file_path, data_name='facs', 
#                  flag_size_factor=True, 
#                  total_ct_per_cell=1e4, 
#                  flag_log1p=True):
#     
#     """load normalized data
#     1. Load filtered data for both FACS and droplet
#     2. Size factor normalization to counts per 1 million (total_ct_per_cell)
#     3. log(x+1) transform
#     4. Combine the data 
# 
#     Args:
#         file_path (str): file path. Should contain both 
#             "tabula-muris-senis-droplet-official-raw-obj.h5ad" and 
#             "tabula-muris-senis-facs-official-raw-obj.h5ad"
# 
#     Returns:
#         adata (AnnData): Combined data for FACS and droplet
#     """
#     
#     if data_name not in ['facs', 'droplet']:
#         raise ValueError('# load_tms: %s need to be one of [facs, droplet]')
#         
#     if data_name=='facs':
#         file_name = 'tabula-muris-senis-facs-official-raw-obj.h5ad'
#     elif data_name=='droplet':
#         file_name = 'tabula-muris-senis-droplet-official-raw-obj.h5ad'
#     else:
#         return None
#             
#     # Load filtered data
#     adata = read_h5ad(f'{file_path}/{file_name}')
#     # Size factor normalization
#     if flag_size_factor == True:
#         # sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#         sc.pp.normalize_per_cell(adata, counts_per_cell_after=total_ct_per_cell)
#     # log(x+1) transform
#     if flag_log1p == True:
#         sc.pp.log1p(adata)
#     # Combine the data 
#     if 'facs' in data_name:
#         ind_select = adata.obs['age'].isin(['3m', '18m', '24m'])
#         adata = adata[ind_select,]
#     # Add extra annotations
#     adata.obs['age_num'] = [int(x.replace('m','')) for x in adata.obs['age']]
#     return adata
