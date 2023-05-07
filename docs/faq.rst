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
    
    
Which GWAS and scRNA-seq data to use?
======================================================

To ensure a reasonable number of scDRS discoveries, we recommend using GWAS data with a heritability z-score greater than 5 or a sample size greater than 100K. We also recommend using scRNA-seq data with a diverse set of cells potentially relevant to disease, although a smaller number of cells should not affect the scDRS power.


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


scDRS detected few significant cells (FDR<0.2)?
==================================================

scDRS may be underpowered for certain GWAS/scRNA-seq data sets. In these cases, the ensuing scDRS group analyses may still have sufficient power, because scDRS group analyses aggregate results of individual cells and hence have higher power than the scDRS individual cell-level analyses. To assess if scDRS has sufficient power, we suggest performing the `scDRS group analyses <https://martinjzhang.github.io/scDRS/reference_cli.html#perform-downstream>`_ to assess significance at an aggregated level. In addition, it is helpful to visually inspect the scDRS normalized disease score on the UMAP plot. Localized enrichments of high scDRS disease scores on the UMAP usually indicate that scDRS have detected interesting biological signals.


MC z-scores are much more significant than MC p-values in group analysis due to the MC limit? 
===========================================================================================

Increasing :code:`--n-ctrl` in `compute-score` will produce more control scores, which will be used in the group analysis to increase the number of MC samples for MC tests. Alternatively, you can compute a p-value from assoc_mcz when assoc_mcp is reasonably small. As mentioned in the Methods section of our manuscript: "We recommend using MC P values to determine statistical significance and using MC z-scores to further prioritize associations whose MC P values have reached the MC limit. "


Use scDRS for other types of single-cell data?
====================================================

scDRS is tailored for scRNA-seq. Best practices for using scDRS on other data types and systematic comparisons with alternative methods remain interesting future directions.

We empirically observed that scDRS works for other types of RNA-seq data like spatial transcriptomics. 

We empirically observed that scDRS works for single-cell DNA methylation data. 

scDRS should work for single-cell ATAC-seq in principle, although you may need some customized Python codes. To do this,

1. Right after calling :code:`scdrs.preprocess`, create a categorical :code:`adata.obs['atac_match']` for your control gene matching criteria by dividing features (genes/peaks) into discrete bins. We recommend >20 features per bin. We recommend matching for mean accessibility and GC contents, as done in `gchromVAR <https://github.com/caleblareau/gchromVAR>`_.
2. When calling :code:`scdrs.score_cell`, tell scDRS to use this matching criteria by :code:`ctrl_match_key='atac_match'`.

Relevant works for individual cell-level associations for scATAC-seq: `Ulirsch et al. Nat Genet 2019 <https://www.nature.com/articles/s41588-019-0362-6>`_, `Chiou et al. Nat Genet 2021 <https://www.nature.com/articles/s41588-021-00823-0>`_, `Yu et al. Nat Biotechnol 2022 <https://www.nature.com/articles/s41587-022-01341-y>`_.




