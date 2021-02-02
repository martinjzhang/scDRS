
"""
This script download the AML dataset from van Galen et al. 2019.
It parses the dataset to be an h5ad file without any further processing.

"""
import scanpy as sc
import numpy as np
from os.path import join
import pandas as pd
import anndata
import glob
from scipy.sparse import csr_matrix

data_dir = "data/"

AML_files = glob.glob(join(data_dir, "*AML*D0.dem.txt.gz"))
healthy_files = glob.glob(join(data_dir, "*BM[1-4].dem.txt.gz"))

AML_donors = [f.split('_')[1].split('-')[0] for f in AML_files]
healthy_donors = [f.split('_')[1].split('.')[0] for f in healthy_files]


adata_list = []
for group, donors in zip(["AML", "healthy"], [AML_donors, healthy_donors]):
    for donor in donors:
        prefix = f"*{donor}"
        if group == "AML":
            prefix += "-D0"
        mat_file = glob.glob(join(data_dir, f"*{prefix}.dem.txt.gz"))
        anno_file = glob.glob(join(data_dir, f"*{prefix}.anno.txt.gz"))
        assert len(mat_file) == 1, len(anno_file) == 1
        mat_file, anno_file = mat_file[0], anno_file[0]
        print(donor, mat_file, anno_file)
        adata = sc.read_csv(mat_file, delimiter='\t').T
        anno = pd.read_csv(anno_file, sep='\t', index_col=0)
        assert np.all(adata.obs_names == anno.index)
        adata.obs = anno
        adata.obs['donor'] = donor
        adata_list.append(adata)

adata = anndata.concat(adata_list)
adata.X = csr_matrix(adata.X)
adata.write_h5ad(join(data_dir, "vangalen2019.h5ad"))