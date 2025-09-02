[![DOI](https://zenodo.org/badge/278546123.svg)](https://zenodo.org/badge/latestdoi/278546123)

<!--# scDRS <!-- omit in toc -->

scDRS (single-cell disease-relevance score) is a method for associating individual cells in single-cell RNA-seq data with disease GWASs, built on top of [AnnData](https://anndata.readthedocs.io/en/latest/) and [Scanpy](https://scanpy.readthedocs.io/en/stable/). 

Read the [documentation](https://martinjzhang.github.io/scDRS/): [installation](https://martinjzhang.github.io/scDRS/index.html#installation), [usage](https://martinjzhang.github.io/scDRS/index.html#usage), [command-line interface (CLI)](https://martinjzhang.github.io/scDRS/reference_cli.html#), [file formats](https://martinjzhang.github.io/scDRS/file_format.html), etc. 

Check out [instructions](https://github.com/martinjzhang/scDRS/issues/2) for making customized gene sets using MAGMA. 

### Reference
[Zhang*, Hou*, et al. "Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data"](https://www.nature.com/articles/s41588-022-01167-z), Nature Genetics, 2022.

### Versions
- [v1.0.4](https://github.com/martinjzhang/scDRS): development version. Convert `adata.X` to `csc` matrix for compute raw score speed up https://github.com/martinjzhang/scDRS/pull/118; 
- [v1.0.3b](https://github.com/martinjzhang/scDRS/releases/tag/v1.0.3b): latest stable version. Fixing a bug of negative values of `ct_mean` when `--adj-prop` and `--cov` are on and there are genes extremely low expression; print `--adj-prop` info in `scdrs compute-score`; check p-value and z-score files that the gene column should have header `GENE`; force index in df_cov and df_score to be str; add --min-genes and --min-cells in CLI for customized filtering; adjustable FDR threshold for plot_group_stats https://github.com/martinjzhang/scDRS/pull/75.
- [v1.0.2](https://github.com/martinjzhang/scDRS/releases/tag/v1.0.2): stable version in 09/22 right after publication. Bug fixes on `scdrs.util.plot_group_stats`; input checks in `scdrs munge-gs` and `scdrs.util.load_h5ad`.
#### Older versions
- [v1.0.1](https://github.com/martinjzhang/scDRS/releases/tag/v1.0.1): stable version used in publication. Identical to `v1.0.0` except documentation.
- [v1.0.0](https://github.com/martinjzhang/scDRS/releases/tag/v1.0.0): stable version used in revision 1. Results are identical to `v0.1` for binary gene sets. Changes with respect to `v0.1`:
    -  scDRS command-line interface (CLI) instead of `.py` scripts for calling scDRS in bash, including `scdrs munge-gs`, `scdrs compute-score`, and `scdrs perform-downstream`.
    -  More efficient in memory use due to the use of sparse matrix throughout the computation.
    -  Allow the use of quantitative weights.
    -  New feature `--adj-prop` for adjusting for cell type-proportions.
- [v0.1](https://github.com/martinjzhang/scDRS/tree/v0.1): stable version used in the initial submission.

### Code and data to reproduce results of the paper
See [scDRS_paper](https://github.com/martinjzhang/scDRS_paper) for more details ([experiments](./experiments) folder is deprecated). Data are at [figshare](https://figshare.com/projects/Single-cell_Disease_Relevance_Score_scDRS_/118902). 
- Download [GWAS gene sets](https://figshare.com/articles/dataset/scDRS_data_release_030122/19312583?file=34300898) (**.gs** files) for 74 diseases and complex traits.
- Download [scDRS results](https://figshare.com/articles/dataset/scDRS_data_release_030122_score_file_tmsfacs/19312607) (**.score.gz** and **.full_score.gz** files) for TMS FACS + 74 diseases/trait.

Older versions
- Initial submission: [GWAS gene sets](https://figshare.com/articles/dataset/scDRS_data_release_092121/16664080?file=30853708) and [scDRS results](https://figshare.com/articles/dataset/scDRS_data_release_092121_score_file_tmsfacs/16664077).


# Explore scDRS results via [CELLxGENE](https://cellxgene.cziscience.com/docs/01__CellxGene)
<!--- 
- Demo for [3 TMS FACS cell types and 3 diseases/traits](https://scdrs-demo.herokuapp.com/).
- Results for [110,096 TMS FACS cells and 74 diseases/traits](https://scdrs-tms-facs.herokuapp.com/).
-->
- [h5ad files](https://figshare.com/articles/dataset/scdrs_cellxgene/15065061) compatible with CELLxGENE
- [Instructions](https://cellxgene.cziscience.com/docs/01__CellxGene) on running CELLxGENE

| <img src="./docs/images/cellxgene.tmsfacs.png" alt="cellxgene" style="width:100%"> | <img src="./docs/images/cellxgene.tmsfacs_ibd.png" alt="cellxgene" style="width:100%"> |
| :--------------------------------------------------------------------------------: | :------------------------------------------------------------------------------------: |
|                   110,096 cells from 120 cell types in TMS FACS                    |                                  IBD-associated cells                                  |



# scDRS scripts (deprecated) 
---
**NOTE:** scDRS scripts are still maintained but deprecated. Consider using [scDRS command-line interface](https://martinjzhang.github.io/scDRS/reference_cli.html) instead.

---

### scDRS script for score calculation 
Input: scRNA-seq data (**.h5ad** file) and gene set file (**.gs** file)

Output: scDRS score file (**{trait}.score.gz** file) and full score file (**{trait}.full_score.gz** file) for each trait in the **.gs** file

```sh
h5ad_file=your_scrnaseq_data
cov_file=your_covariate_file
gs_file=your_gene_set_file
out_dir=your_output_folder

python compute_score.py \
    --h5ad_file ${h5ad_file}.h5ad\
    --h5ad_species mouse\
    --cov_file ${cov_file}.cov\
    --gs_file ${gs_file}.gs\
    --gs_species human\
    --flag_filter True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder ${out_dir}
```

- `--h5ad_file` (**.h5ad** file) : scRNA-seq data
- `--h5ad_species` ("hsapiens"/"human"/"mmusculus"/"mouse") : species of the scRNA-seq data samples
- `--cov_file` (**.cov** file) : covariate file (optional, .tsv file, see [file format](https://martinjzhang.github.io/scDRS/file_format.html))
- `--gs_file` (**.gs** file) : gene set file (see [file format](https://martinjzhang.github.io/scDRS/file_format.html))
- `--gs_species` ("hsapiens"/"human"/"mmusculus"/"mouse") : species for genes in the gene set file 
- `--flag_filter` ("True"/"False") : if to perform minimum filtering of cells and genes
- `--flag_raw_count` ("True"/"False") : if to perform normalization (size-factor + log1p)
- `--n_ctrl` (int) : number of control gene sets (default 1,000)
- `--flag_return_ctrl_raw_score` ("True"/"False") : if to return raw control scores
- `--flag_return_ctrl_norm_score` ("True"/"False") : if to return normalized control scores
- `--out_folder` : output folder. Score files will be saved as `{out_folder}/{trait}.score.gz` (see [file format](https://martinjzhang.github.io/scDRS/file_format.html))

### scDRS script for downsteam applications

Input: scRNA-seq data (**.h5ad** file), gene set file (**.gs** file), and scDRS full score files (**.full_score.gz** files)

Output: **{trait}.scdrs_ct.{cell_type}** file (same as the new **{trait}.scdrs_group.{cell_type}** file) for cell type-level analyses (association and heterogeneity); **{trait}.scdrs_var** file (same as the new **{trait}.scdrs_cell_corr** file) for cell variable-disease association; **{trait}.scdrs_gene** file for disease gene prioritization.

```sh
h5ad_file=your_scrnaseq_data
out_dir=your_output_folder
python compute_downstream.py \
    --h5ad_file ${h5ad_file}.h5ad \
    --score_file @.full_score.gz \
    --cell_type cell_type \
    --cell_variable causal_variable,non_causal_variable,covariate\
    --flag_gene True\
    --flag_filter False\
    --flag_raw_count False\ # flag_raw_count is set to `False` because the toy data is already log-normalized, set to `True` if your data is not log-normalized
    --out_folder ${out_dir}
```

- `--h5ad_file` (**.h5ad** file) : scRNA-seq data
- `--score_file` (**.full_score.gz** files) : scDRS full score files; supporting use of "@" to match strings
- `--cell_type` (str) : cell type column (supporting multiple columns separated by comma); must be present in `adata.obs.columns`; used for cell type-disease association analyses (5% quantile as test statistic) and detecting association heterogeneity within cell type (Geary's C as test statistic)
- `--cell_variable` (str) : cell-level variable columns (supporting multiple columns separated by comma); must be present in `adata.obs.columns`; used for cell variable-disease association analyses (Pearson's correlation as test statistic)
- `--flag_gene` ("True"/"False") : if to correlate scDRS disease scores with gene expression
- `--flag_filter` ("True"/"False") : if to perform minimum filtering of cells and genes
- `--flag_raw_count` ("True"/"False") : if to perform normalization (size-factor + log1p)
- `--out_folder` : output folder. Score files will be saved as `{out_folder}/{trait}.scdrs_ct.{cell_type}` for cell type-level analyses (association and heterogeneity); `{out_folder}/{trait}.scdrs_var` file for cell variable-disease association; `{out_folder}/{trait}.scdrs_var.{trait}.scdrs_gene` file for disease gene prioritization. (see [file format](https://martinjzhang.github.io/scDRS/file_format.html))
    

<!---
# File formats
## scDRS input files
**.h5ad** file (scRNA-seq data, compatible with [scanpy](https://scanpy.readthedocs.io/en/stable/index.html))
            
**.cov** file (.tsv file) : covariate file for the scRNA-seq data
- index: cell names, consistent with `adata.obs_names`
- other comlumns: covariate with numerical values

  |          index          | const | n_genes | sex_male |  age  |
  | :---------------------: | :---: | :-----: | :------: | :---: |
  | A10_B000497_B009023_S10 |   1   |  2706   |    1     |  18   |
  | A10_B000756_B007446_S10 |   1   |  3212   |    1     |  18   |
        
**.gs** file (.tsv file) : one can use [precomputed GWAS gene sets](https://figshare.com/articles/dataset/scDRS_data_release_092121/16664080?file=30853708) for 74 diseases and complex traits or [compute gene sets directly from GWAS](docs/compute_magma_gs.md).
- TRAIT: trait name
- GENESET: a comma-separated string of genes 

  |           TRAIT           |         GENESET          |
  | :-----------------------: | :----------------------: |
  |        PASS_HbA1C         |   FN3KRP,FN3K,HK1,GCK    |
  | PASS_MedicationUse_Wu2019 | FTO,SEC16B,ADCY3,DNAJC27 |

## scDRS output files
**{trait}.score.gz** file : scDRS score file
- index: cell names, should be the same as adata.obs_names
- raw_score: raw disease score
- norm_score: normalized disease score
- mc_pval: cell-level MC p-value
- pval: cell-level scDRS p-value
- nlog10_pval: -log10(pval)
- zscore: z-score converted from pval

  |          index          | raw_score | norm_score | mc_pval |  pval   | nlog10_pval | zscore |
  | :---------------------: | :-------: | :--------: | :-----: | :-----: | :---------: | :----: |
  | A10_B000497_B009023_S10 |   0.730   |    7.04    | 0.0476  | 0.00166 |    2.78     |  2.94  |
  | A10_B000756_B007446_S10 |   0.725   |    7.30    | 0.0476  | 0.00166 |    2.78     |  2.94  |
        
**{trait}.full_score.gz** file : scDRS full score file
- All contents of **{trait}.score.gz** file
- ctrl_norm_score_{i_ctrl} : raw control scores (specified by `--flag_return_ctrl_raw_score True`) 
- ctrl_norm_score_{i_ctrl} : normalized control scores (specified by `--flag_return_ctrl_norm_score True`) .

**{trait}.scdrs_ct.{cell_type}** file : cell type-level analysis (association and heterogeneity)
- {trait} : trait name consistent with **{trait}.full_score.gz** file
- {cell_type} : the cell type column in `adata.obs.columns`
- first column: cell types in the `adata.obs` `cell_type` column
- n_cell: number of cells from the cell type
- n_ctrl: number of control gene sets
- assoc_mcp: MC p-value for cell type-disease association
- assoc_mcz: MC z-score for cell type-disease association
- hetero_mcp:  MC p-value for within-cell type heterogeneity in association with disease
- hetero_mcz:  MC z-score for within-cell type heterogeneity in association with disease

  |                 | n_cell | n_ctrl | assoc_mcp  | assoc_mcz  | hetero_mcp | hetero_mcz |
  | :-------------: | :----: | :----: | :--------: | :--------: | :--------: | :--------: |
  |   causal_cell   |  10.0  |  20.0  | 0.04761905 | 12.297529  |    1.0     |    1.0     |
  | non_causal_cell |  20.0  |  20.0  | 0.9047619  | -1.1364214 |    1.0     |    1.0     |

**{trait}.scdrs_var** file : cell-level variable analysis (association to disease)
- {trait} : trait name consistent with **{trait}.full_score.gz** file
- first column: all cell-level variables (specified in `--cell_variable`)
- n_ctrl: number of control gene sets
- corr_mcp: MC p-value for cell-level variable association with disease
- corr_mcz: MC z-score for cell-level variable association with disease

  |                     | n_ctrl |  corr_mcp  | corr_mcz  |
  | :-----------------: | :----: | :--------: | :-------: |
  |   causal_variable   |  20.0  | 0.04761905 | 3.4574268 |
  | non_causal_variable |  20.0  | 0.23809524 | 0.8974108 |
  |      covariate      |  20.0  | 0.1904762  | 1.1220891 |

**{trait}.scdrs_gene** file : correlation of gene expression with scDRS disease score
- {trait} : trait name consistent with **{trait}.full_score.gz** file
- first column: genes in `adata.var_names`
- CORR: correlation with scDRS disease score across all cells in adata
- RANK: correlation of correlation across genes

  |  index   | CORR  | RANK  |
  | :------: | :---: | :---: |
  | Serping1 | 0.314 |   0   |
  |   Lmna   | 0.278 |   1   |

-->
