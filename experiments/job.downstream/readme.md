# Downstream analysis with TRS
- `celltype_association.ipynb`: TRS can identify trait - celltype interesting association pairs by aggregating the cell-wise p-values
- `celltype_hetero.ipynb`: TRS can identify heterogeneous trait signal patterns within a cell type.

## Cell type association

## Cell type heterogeneity
We first compute the z-score for all the control scores as well, using the following code
```python
from scipy.stats import rankdata
import scipy as sp
import pandas as pd
import scipy as sp
import numpy as np
from os.path import join
import submitit
import scTRS.method as md

def compute_control_zsc(prefix, out):
    score = pd.read_csv(f"{prefix}.score.gz", sep='\t')
    full_score = pd.read_csv(f"{prefix}.full_score.gz", sep='\t')

    control_cols = [col for col in full_score.columns if col.startswith("ctrl_norm_score_")]
    v_norm_score = full_score["norm_score"].values
    mat_ctrl_norm_score = full_score[control_cols].values
    n_ctrl = len(control_cols)

    cat_score = np.vstack([v_norm_score, mat_ctrl_norm_score.T])
    rank_score = rankdata(cat_score, axis=0)
    mat_mc_pval = 1 - (rank_score - 1) / rank_score.shape[0]
    mat_pooled_p = md.get_p_from_empi_null(cat_score.flatten(), mat_ctrl_norm_score.flatten()).reshape(cat_score.shape).T
    mat_pooled_z = -sp.stats.norm.ppf(mat_pooled_p).clip(min=-10,max=10)

    df = pd.DataFrame(index=score["index"], data=mat_pooled_z, dtype=np.float32)
    df.columns = ["trait", *[f"ctrl_{i}" for i in range(n_ctrl)]]

    df.to_csv(out, sep='\t', index=True, compression='gzip')
    
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

DATA_PATH = '/n/holystore01/LABS/price_lab/Users/mjzhang/scTRS_data'
score_file_dir = join(DATA_PATH, "score_file/score.facs.gwas_max_abs_z.top500")
out_dir = "out/celltype_hetero"
compute_func = lambda trait: compute_control_zsc(join(score_file_dir, trait), join(out_dir, f"{trait}.ctrl_zsc.gz"))

executor = submitit.AutoExecutor(folder="~/submitit_log/")
executor.update_parameters(timeout_min=30, mem_gb=12, slurm_partition="serial_requeue")
jobs = executor.map_array(compute_func, trait_list)
```
