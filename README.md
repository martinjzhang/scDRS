# scDRS
Single-cell disease-relevance score.

### Explore results
- Results for 74 diseases/traits and the TMS FACS data ([cellxgene visualization](https://scdrs-tms-facs.herokuapp.com/))
- Demo for 3 diseases/traits and 3 TMS FACS cell types ([cellxgene visualization](https://scdrs-demo.herokuapp.com/))

# Installation
- Go to the scDRS directory
- Run `pip install -r requirements.txt; pip install -e .`

# Usage 
- Compute scDRS results for a given scRNA-seq data set (.h5ad file) and multiple (GWAS) gene sets (.gs file): compute_score.py
```
python compute_score.py \
    --h5ad_file $H5AD_FILE\ # 
    --h5ad_species $SPECIES\
    --gs_file $GS_FILE\
    --gs_species $SPECIES\
    --flag_filter True\
    --flag_raw_count True\
    --n_ctrl 1000\
    --flag_return_ctrl_raw_score False\
    --flag_return_ctrl_norm_score True\
    --out_folder $OUT_FOLDER
```

- compute_downstream.py Script for scDRS-based downsteam analyses. 

## File formats
- .h5ad file (compatible with [scanpy](https://scanpy.readthedocs.io/en/stable/index.html))

- .gs file: .tsv file

    1. TRAIT: trait name
    2. GENESET: a comma-separated string of genes 

        Example:
    
        |           TRAIT           |         GENESET          |
        | :-----------------------: | :----------------------: |
        |        PASS_HbA1C         |   FN3KRP,FN3K,HK1,GCK    |
        | PASS_MedicationUse_Wu2019 | FTO,SEC16B,ADCY3,DNAJC27 |
            
- .cov file: .tsv file

    1. index: cell names, should be the same as adata.obs_names
    2. cov1: numerical-values covariate
    3. cov2: numerical-values covariate

        Example:
        |                 index                 | const | n_genes | sex_male |  age  |
        | :-----------------------------------: | :---: | :-----: | :------: | :---: |
        | A10_B000497_B009023_S10.mm10-plus-0-0 |   1   |  2706   |    1     |  18   |
        | A10_B000756_B007446_S10.mm10-plus-0-0 |   1   |  3212   |    1     |  18   |
  
- .score.gz file:
 
    1. index: cell names, should be the same as adata.obs_names
    2. raw_score: raw disease score
    3. norm_score: normalized disease score
    3. mc_pval: cell-level MC p-value
    3. pval: cell-level scDRS p-value
    3. nlog10_pval: -log10(pval)
    3. zscore: z-score converted from pval

        Example:
        |                 index                 | raw_score  | norm_score |  mc_pval   |     pval     | nlog10_pval |  zscore   |
        | :-----------------------------------: | :--------: | :--------: | :--------: | :----------: | :---------: | :-------: |
        | A10_B000497_B009023_S10.mm10-plus-0-0 | 0.7298449  | 7.0396357  | 0.04761905 | 0.0016638935 |  2.7788744  | 2.9357162 |
        | A10_B000756_B007446_S10.mm10-plus-0-0 | 0.72515404 |  7.300498  | 0.04761905 | 0.0016638935 |  2.7788744  | 2.9357162 |
