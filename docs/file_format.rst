File formats
============

.sumstats
~~~~~~~~~
GWAS summary statistics following the `LDSC format <https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format>`_.
    
.. csv-table:: Example .sumstats file
   :header: "GENE", "BMI", "HEIGHT"
   :delim: space
   
   SNP A1 A2 N CHISQ Z
   rs7899632 A G 59957 3.4299 -1.852
   rs3750595 A C 59957 3.3124 1.82

.h5ad
~~~~~

Single-cell data :code:`.h5ad` file as defined in `AnnData <https://anndata.readthedocs.io/en/latest/>`_ and `Scanpy <https://scanpy.readthedocs.io/en/stable/>`_.


pval_file,zscore_file
~~~~~~~~~~~~~~~~~~~~~
GWAS gene-level p-values / z-scores for different traits. A :code:`.tsv` file with first column corresponding to genes and other columns corresponding to p-values / z-scores of traits (one trait per column).
    
.. csv-table:: Example pval_file
   :header: "GENE", "BMI", "HEIGHT"
   :delim: space
   
   OR4F5   0.001  0.01
   DAZ3    0.01   0.001

.gs
~~~~

scDRS gene set file. A :code:`.tsv` file with two columns :code:`["TRAIT", "GENESET"]` and one line per trait. Can be generated using customized code or from p-value or z-score files using scDRS CLI :code:`scdrs munge-gs`.

TRAIT
    Trait (gene set) identifier.
GENESET
    Comma-separated list of gene-weight pairs with the form "gene1\:weight1,gene2\:weight2,..." 
    or "gene1,gene2,..." (meaning weights are 1). 


.. csv-table:: Example weighted .gs file
   :header: "TRAIT", "GENESET"
   :delim: space
   :align: center
   :width: 50%
   
   PASS_HbA1C FN3KRP:1.2,FN3K:2.3,HK1:4.7,GCK:5.2
   PASS_MedicationUse_Wu2019 FTO:3,SEC16B:0.6,ADCY3:1.5,DNAJC27:1.3

.. csv-table:: Example unweighted .gs file
   :header: "TRAIT", "GENESET"
   :delim: space
   :align: center
   :width: 50%
   
   PASS_HbA1C FN3KRP,FN3K,HK1,GCK
   PASS_MedicationUse_Wu2019 FTO,SEC16B,ADCY3,DNAJC27
  

.cov
~~~~

scDRS covariate file for the :code:`.h5ad` single-cell data. :code:`.tsv` file.

- First column: cell names, consistent with :code:`adata.obs_names`.
- Other comlumns: covariates with numerical values.

.. csv-table:: Example .cov file
   :header: "index", "const", "n_genes", "sex_male", "age"
   :align: center
   :width: 50%
   
   A10_B000497_B009023_S10, 1, 2706, 1, 18 
   A10_B000497_B009023_S10, 1, 2501, 0, 24 
   

<trait>.score.gz
~~~~~~~~~~~~~~~~

scDRS score file for a give trait. :code:`.tsv.gz` file.

- First column: cell names, should be the same as :code:`adata.obs_names`.
- raw_score: raw disease score.
- norm_score: normalized disease score.
- mc_pval: cell-level MC p-value.
- pval: cell-level scDRS p-value.
- nlog10_pval: -log10(pval).
- zscore: z-score converted from pval.

.. csv-table:: Example <trait>.score.gz file
   :header: "index", "raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval", "zscore"
  
   A10_B000497_B009023_S10, 0.730, 7.04, 0.0476, 0.00166, 2.78, 2.94
   A10_B000756_B007446_S10, 0.725, 7.30, 0.0476, 0.00166, 2.78, 2.94
   
        
<trait>.full_score.gz
~~~~~~~~~~~~~~~~~~~~~

scDRS full score file for a give trait. :code:`.tsv.gz` file.

- All columns of :code:`{trait}.score.gz` file.
- ctrl_raw_score_<i_ctrl> : raw control scores, specified by :code:`--flag_return_ctrl_raw_score True`.
- ctrl_norm_score_<i_ctrl> : normalized control scores, specified by :code:`--flag_return_ctrl_norm_score True`.


<trait>.scdrs_group.<annot>
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Results for scDRS group-level analysis for a give trait and a given cell-group annotation (e.g., cell type). :code:`.tsv` file.

- <trait> : trait name consistent with :code:`<trait>.full_score.gz` file.
- <annot> : cell-annotation in :code:`adata.obs.columns`, specified by :code:`group_analysis` in CLI.
- First column: different cell groups in :code:`adata.obs[<annot>]`.
- n_cell: number of cells from the cell group.
- n_ctrl: number of control gene sets.
- assoc_mcp: MC p-value for cell group-disease association.
- assoc_mcz: MC z-score for cell group-disease association.
- hetero_mcp:  MC p-value for within-cell group heterogeneity in association with disease.
- hetero_mcz:  MC z-score for within-cell group heterogeneity in association with disease.

.. csv-table:: Example <trait>.scdrs_group.<annot> file
   :header: "", "n_cell", "n_ctrl", "assoc_mcp", "assoc_mcz", "hetero_mcp", "hetero_mcz"
   
   causal_cell    , 10.0,   20.0, 0.04761905, 12.297529 , 1.0, 1.0
   non_causal_cell, 20.0,   20.0, 0.9047619 , -1.1364214, 1.0, 1.0


<trait>.scdrs_cell_corr
~~~~~~~~~~~~~~~~~~~~~~~

Results for scDRS cell-level correlation analysis for a given trait. :code:`.tsv` file.

- <trait> : trait name consistent with :code:`<trait>.full_score.gz` file.
- First column: all cell-level variables, specified by specified by :code:`corr_analysis` in CLI.
- n_ctrl: number of control gene sets.
- corr_mcp: MC p-value for cell-level variable association with disease.
- corr_mcz: MC z-score for cell-level variable association with disease.

.. csv-table:: Example <trait>.scdrs_cell_corr file
   :header: "", "n_cell", "corr_mcp", "corr_mcz"
   
   causal_variable    , 20.0, 0.04761905, 3.4574268
   non_causal_variable, 20.0, 0.23809524, 0.8974108
   covariate          , 20.0, 0.1904762 , 1.1220891

<trait>.scdrs_gene
~~~~~~~~~~~~~~~~~~

Results for scDRS gene-level correlation analysis for a given trait. :code:`.tsv` file.

- <trait> : trait name consistent with :code:`<trait>.full_score.gz` file.
- First column: genes in :code:`adata.var_names`.
- CORR: correlation with scDRS disease score across all cells in :code:`adata`.
- RANK: rank of correlation across genes (starting from 0).

.. csv-table:: Example <trait>.scdrs_gene file
   :header: "index", "CORR", "RANK"   

   Serping1, 0.314, 0
     Lmna  , 0.278, 1