import scanpy as sc
import numpy as np
from os.path import join
import pandas as pd
import anndata

data_dir = "data/"

mat1 = sc.read_mtx(join(data_dir, "GSM2560248_2.1.mtx.gz")).T
mat2 = sc.read_mtx(join(data_dir, "GSM2560249_2.2.mtx.gz")).T

cell_md = pd.read_csv(join(data_dir, "GSE96583_batch2.total.tsne.df.tsv.gz"), sep='\t').rename(columns={"cell": "celltype"})
gene_md = pd.read_csv(join(data_dir, "GSE96583_batch2.genes.tsv.gz"), sep='\t', names=["ENSEMBL", "SYMBOL"])

cell_md["cluster"] = cell_md["cluster"].astype("category")
cell_md["ind"] = cell_md["ind"].astype("category")

gene_md.index = gene_md["SYMBOL"].values

adata = anndata.concat([mat1, mat2])
adata.obs = cell_md
adata.var = gene_md
adata.var_names = adata.var["SYMBOL"].values

sc.pp.filter_genes(adata, min_cells=10)
# remove duplicated genes
adata = adata[:, ~adata.var_names.duplicated()]
adata.write_h5ad(join(data_dir, "kang2018.h5ad"))
