# scDRS
Single-cell disease-relevance score

# Installation
- Go to the scTRS directory
- Run `pip install -e .`

# Usage 
- compute_score.py
Script for batch processing scores. 

## File formats
- .gs file: .tsv file

    1. TRAIT: trait name
    2. GENESET: a comma-separated string of genes 

        Example:
    
        | TRAIT | GENESET |
        | :---: | :---: |
        | PASS_HbA1C | FN3KRP,FN3K,HK1,GCK |
        | PASS_MedicationUse_Wu2019 | FTO,SEC16B,ADCY3,DNAJC27 |
            
- .cov file: .tsv file

    1. index: cell names, should be the same as adata.obs_names
    2. cov1: numerical-values covariate
    3. cov2: numerical-values covariate

        Example:
        | index | const | n_genes | sex_male | age |
        | :---: | :---: | :---: | :---: | :---: |
        | A10_B000497_B009023_S10.mm10-plus-0-0 | 1 | 2706 | 1 | 18 |
        | A10_B000756_B007446_S10.mm10-plus-0-0 | 1 | 3212 | 1 | 18 |
  
- .score file:
