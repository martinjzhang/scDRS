import scTRS.util as util
import scTRS.method as md
import scTRS.data_loader as dl
import submitit
from anndata import read_h5ad
import numpy as np
import pandas as pd
from os.path import join, exists
import itertools
import glob

def load_adata(dataset):
    """
    Given the dataset, return the AnnData
    
    Args
    ---
    dataset: dataset name, can be ["tms_droplet.{tissue}", "tms_facs.{tissue}", "mouse_liver_halpern", "liver_atlas"]
    """
    if dataset.startswith("tms_droplet."):
        tissue = dataset.split('.')[1]
        adata = read_h5ad(f"/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tabula-muris-senis-droplet-processed-official-annotations-{tissue}.h5ad")
    elif dataset.startswith("tms_facs."):
        tissue = dataset.split('.')[1]
        adata = read_h5ad(f"/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tabula-muris-senis-facs-processed-official-annotations-{tissue}.h5ad")
    elif dataset == "mouse_liver_halpern":
        adata = read_h5ad("/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/single_cell_data/mouse_liver_halpern_nature_2017/obj_raw.h5ad")
    elif dataset == "mouse_liver_halpern_full":
        adata = read_h5ad("/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/single_cell_data/mouse_liver_halpern_nature_2017/obj_raw_full.h5ad")
    elif dataset == "liver_atlas":
        adata = read_h5ad("/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/single_cell_data/liver_atlas/obj_processed.h5ad")
        adata.obs["cell_ontology_class"] = adata.obs["celltype"] 
    else:
        raise NotImplementedError
    return adata

def load_score(root_dir, dataset, trait):
    """
    Given the dataset, return the AnnData
    
    Args
    ---
    dataset: dataset name, can be ["tms_droplet.{tissue}", "tms_facs.{tissue}", "mouse_liver_halpern", "liver_atlas"]
    trait: trait name
    """
    
    if dataset.startswith("tms_droplet."):
        score_dir = "score.tms_droplet.gwas_max_abs_z.top500.weight_1en2"
    elif dataset.startswith("tms_facs."):
        score_dir = "score.tms_facs.gwas_max_abs_z.top500"
    elif dataset in ["mouse_liver_halpern", "mouse_liver_halpern_full", "liver_atlas"]:
        score_dir = f"score.{dataset}.gwas_max_abs_z.top500"
    else:
        raise NotImplementedError
    score = pd.read_csv(join(root_dir, score_dir, f"{trait}.full_score.gz"), sep='\t', index_col=0)
    return score

def wrapper(trait, dataset, out_prefix):
    DATA_PATH = '/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data'
    SCORE_FILE_DIR = join(DATA_PATH, "score_file")
    
    adata = load_adata(dataset)
    score = load_score(SCORE_FILE_DIR, dataset, trait)
    trait_col = "norm_score"
    ctrl_cols = [col for col in score.columns if col.startswith(f"ctrl_norm_score_")]
    
    # takes the relevant cells
    
    score = score.reindex(adata.obs.index).dropna()
    
    zsc_dict = dict()
    zsc_dict[trait_col] = score[trait_col]
    for col in ctrl_cols:
        zsc_dict[col] = score[col]

    trs_stats = util.calculate_trs_stats(zsc_dict=zsc_dict, zsc_index=score.index, 
                                    stats_dict = {"mean": np.mean, "sd": np.std, "gearysc": None},
                                    adata=adata, stratify_by="cell_ontology_class")
    for name in trs_stats:
        trs_stats[name].to_csv(out_prefix + f".{name}.gz", index=True, sep='\t', compression='gzip')
        

executor = submitit.AutoExecutor(folder="~/submitit_log/")
executor.update_parameters(timeout_min=50, mem_gb=16, slurm_partition="serial_requeue")
out_dir = "out/celltype_hetero/group_stats"

trait_list = ['PASS_Schizophrenia_Ruderfer2018',
              'PASS_BipolarDisorder_Ruderfer2018',
              'PASS_Alzheimers_Jansen2019', 
              'PASS_AdultOnsetAsthma_Ferreira2019',
              'PASS_Coronary_Artery_Disease', 
              'PASS_LargeArteryStroke_Malik2018', 
              'PASS_HDL', 'PASS_LDL',
              'PASS_Rheumatoid_Arthritis', 'PASS_Lupus', 
              'PASS_FastingGlucose_Manning',
              'PASS_IBD_deLange2017', 
              'PASS_Type_1_Diabetes', 
              'PASS_Type_2_Diabetes']

dataset_list = []
# tms_facs
file_list = glob.glob("/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tabula-muris-senis-facs-processed-official-annotations-*")
dataset_list.extend(["tms_facs." + f.split('annotations-')[-1].split('.')[0] for f in file_list])
# tms_droplet
file_list = glob.glob("/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data/tabula_muris_senis/tabula-muris-senis-droplet-processed-official-annotations-*")
dataset_list.extend(["tms_droplet." + f.split('annotations-')[-1].split('.')[0] for f in file_list])
# rest of single tissue dataset
# dataset_list.extend(["mouse_liver_halpern", "mouse_liver_halpern_full", "liver_atlas"])
dataset_list.extend(["liver_atlas"])
        
submit_func = lambda trait_dataset: wrapper(trait_dataset[0], trait_dataset[1], join(out_dir, '@'.join(trait_dataset)))
trait_dataset_list = list(itertools.product(trait_list, dataset_list))
todo_trait_dataset_list = [trait_dataset for trait_dataset in trait_dataset_list if not exists(join(out_dir, '@'.join(trait_dataset) + ".mean.gz"))]

jobs = executor.map_array(submit_func, todo_trait_dataset_list)