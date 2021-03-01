# scTRS
Trait-relevance score informed scRNA-seq analysis

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
    
            TRAIT                        GENESET
            PASS_HbA1C                   FN3KRP,FN3K,HK1,GCK
            PASS_MedicationUse_Wu2019    FTO,SEC16B,ADCY3,DNAJC27
  
- .score file:
