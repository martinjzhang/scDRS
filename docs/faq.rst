FAQ
################################

Here are some frequently asked questions about scDRS.


Typical scDRS workflow
======================

1. Create scDRS gene set file (:code:`.gs`) from GWAS 

    a. Generate gene-level p-values/z-scores for a given trait (pval_file/zscore_file) using `MAGMA <https://ctg.cncr.nl/software/magma>`_ from GWAS summary statistics.
    b. Create :code:`.gs` file from pval_file/zscore_file using scDRS CLI :code:`scdrs munge-gs`.
   
2. Compute scDRS individual cell-level results.

    a. Use scDRS CLI :code:`scdrs compute-score`
    b. Inputs: scDRS gene set file (:code:`.gs`), scRNA-seq data (:code:`.h5ad`), covariates (:code:`.cov`). 
    c. Outputs: scDRS score file (:code:`<trait>.score.gz`), scDRS full score file (:code:`<trait>.full_score.gz`).
    
3. Perform scDRS downstream cell-group analyses.

    a. Use scDRS CLI :code:`scdrs perform-downstream` for a) cell group-trait association, b) within-cell group association heterogeneity, c) correlating scDRS disease score with a cell variable, d) correlating scDRS disease score with gene expression.
    b. Use customized code for other analyses:
    
        i. Compute a test statistic using the normalized disease score across a given set of cells.
        ii. Compute the same test statistic using each of the :code:`n_ctrl` sets of normalized control scores across the same given set of cells.
        iii. MC-pvalue = :code:`(# of control statistics exceeding disease statistic) / (n_ctrl+1)`
    b. Input: scDRS full score file (:code:`<trait>.full_score.gz`), cell annotations stored in :code:`adata.obs` of scRNA-seq :code:`.h5ad` file.
    c. Output: test statistics/p-values based on the scDRS MC tests.


How to create MAGMA gene sets?
==============================

Please see the `instructions <https://github.com/martinjzhang/scDRS/issues/2>`_


Use scDRS for other gene sets?
=====================================

Yes, you can use other gene sets instead of GWAS gene sets with scDRS to identify cells in scRNA-seq with excess expression of genes in the gene set.


Requirement of gene set size?
========================================

The gene set should have a moderate size (e.g., >50 genes and <20% of all genes) for the scDRS results to be statistically valid. In practice, we observe a reasonable performance as long as the gene set size is >=10. Please see details in Methods in Zhang & Hou et al. Nat Genet 2022. 


Computational complexity?
====================================

scDRS scales linearly with the number of cells and number of control gene sets for both computation and memory use. It takes around 3 hours and 60GB for a single-cell data set with a million cells). Please see details in Methods in Zhang & Hou et al. Nat Genet 2022. 


Use scDRS for other types of single-cell data?
====================================================

scDRS is tailored for scRNA-seq. Best practices for using scDRS on other data types and systematic comparisons with alternative methods remain interesting future directions.

We empirically observed that scDRS works for other types of RNA-seq data like spatial transcriptomics. 

We empirically observed that scDRS works for single-cell DNA methylation data. 

scDRS should work for single-cell ATAC-seq in principle, although you may need some customized Python codes. To do this,

1. Right after calling :code:`scdrs.preprocess`, create a categorical :code:`adata.obs['atac_match']` for your control gene matching criteria by dividing features (genes/peaks) into discrete bins. We recommend >20 features per bin. We recommend matching for mean accessibility and GC contents, as done in `gchromVAR <https://github.com/caleblareau/gchromVAR>`_.
2. When calling :code:`scdrs.score_cell`, tell scDRS to use this matching criteria by :code:`ctrl_match_key='atac_match'`.

Relevant works for individual cell-level associations for scATAC-seq: `Ulirsch et al. Nat Genet 2019 <https://www.nature.com/articles/s41588-019-0362-6>`_, `Chiou et al. Nat Genet 2021 <https://www.nature.com/articles/s41588-021-00823-0>`_, `Yu et al. Nat Biotechnol 2022 <https://www.nature.com/articles/s41587-022-01341-y>`_.




