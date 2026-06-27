import scdrs
import pandas as pd
import scanpy as sc
from anndata import AnnData
import os, subprocess, tempfile

CORTEX_URL = "https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt"

# save to files
if not os.path.exists("data/"):
    os.makedirs("data/")

# download gene set
with tempfile.TemporaryDirectory() as tmpdir:
    subprocess.check_call(
        f"wget https://figshare.com/ndownloader/files/34300898 -O {tmpdir}/geneset.zip "
        + f"&& unzip {tmpdir}/geneset.zip -d {tmpdir} "
        + f"&& mv {tmpdir}/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs data/gwas.gs",
        shell=True,
    )

with tempfile.TemporaryDirectory() as tmpdir:
    subprocess.check_call(
        f"wget https://figshare.com/ndownloader/files/34300925 -O {tmpdir}/data.zip "
        + f"&& unzip {tmpdir}/data.zip -d {tmpdir} "
        + f"&& mv {tmpdir}/single_cell_data/zeisel_2015/geneset.gs data/spatial.gs",
        shell=True,
    )

# meta information
df_meta = pd.read_csv(CORTEX_URL, nrows=10, sep="\t", header=None)
df_meta = df_meta.iloc[:, 1:].T
columns = df_meta.iloc[0, :]
df_meta = df_meta.iloc[1:, :]
df_meta.columns = columns
df_meta = df_meta.set_index("cell_id")
df_meta.columns.name = None
df_meta["total mRNA mol"] = df_meta["total mRNA mol"].astype(float)

# expression information
df_expr = pd.read_csv(CORTEX_URL, skiprows=11, sep="\t", header=None).set_index(0)
df_expr.index.name = "gene"
# 1st column in backspin cluster, we don't need it here
df_expr = df_expr.iloc[:, 1:]
df_expr.columns = df_meta.index
df_expr = df_expr.T
raw_adata = AnnData(df_expr, obs=df_meta)

# assemble AnnData
sc.pp.filter_cells(raw_adata, min_genes=0)
sc.pp.filter_genes(raw_adata, min_cells=30)

raw_adata.raw = raw_adata

sc.pp.normalize_total(raw_adata, target_sum=1e4)
sc.pp.log1p(raw_adata)
sc.pp.highly_variable_genes(raw_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

raw_adata = raw_adata[:, raw_adata.var.highly_variable]

sc.pp.scale(raw_adata, max_value=10)
sc.tl.pca(raw_adata, svd_solver="arpack")
sc.pp.neighbors(raw_adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(raw_adata)
sc.tl.leiden(raw_adata)

adata = raw_adata.raw.to_adata()
adata.obsp = raw_adata.obsp
adata.X = adata.X

# compile covariates
df_cov = pd.DataFrame(index=adata.obs.index)
df_cov["const"] = 1
df_cov["n_genes"] = adata.obs["n_genes"]

# format gene sets
dict_gwas_gs = scdrs.util.load_gs(
    "data/gwas.gs", src_species="human", dst_species="mouse"
)
dict_spatial_gs = scdrs.util.load_gs("data/spatial.gs")

dict_gs = {
    "SCZ": dict_gwas_gs["PASS_Schizophrenia_Pardinas2018"],
    "Height": dict_gwas_gs["UKB_460K.body_HEIGHTz"],
    "Dorsal": dict_spatial_gs["spatial_dorsal"],
}

adata.write_h5ad("data/expr.h5ad")
df_cov.to_csv("data/cov.tsv", sep="\t")
scdrs.util.save_gs("data/geneset.gs", dict_gs)