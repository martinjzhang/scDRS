CLI
===
scDRS CLI supports 3 main functions:

1. Munge gene sets :code:`scdrs munge-gs`, producing :code:`.gs` file.
   
2. Compute scores :code:`scdrs compute-score`, producing :code:`<trait>.score.gz` and :code:`<trait>.full_score.gz` files.
   
3. Perform downstream analyses :code:`scdrs perform-downstream`, producing :code:`<trait>.scdrs_group.<annot>`, :code:`<trait>.scdrs_cell_corr`, and :code:`<trait>.scdrs_gene` files.

See `scDRS file formats <file_format.html>`_. Use :code:`--help` to check out more CLI details, e.g., :code:`scdrs compute-score --help`. 

scDRS CLI does not distinguish between dash "-" and underscore "_". E.g., followings are equivalent:

.. code-block:: bash
    
    # All dashes
    scdrs compute-score --h5ad-file <h5ad_file> --out-folder <out_folder> ...
    # All underscores
    scdrs compute_score --h5ad_file <h5ad_file> --out_folder <out_folder> ...
    # Mixed
    scdrs compute_score --h5ad-file <h5ad_file> --out_folder <out_folder> ... 
    
    
munge-gs
~~~~~~~~

Convert a :code:`.tsv` GWAS gene statistics file to an scDRS :code:`.gs` file.

1. Read result from a :code:`.tsv` p-value or z-score file.
2. Select a subset of genes for each trait:

   - If both :code:`fdr` and :code:`fwer` are :code:`None`, select the top :code:`n_max` genes.
   - If :code:`fdr` is not :code:`None`, select genes for based on FDR (across all genes of a given trait) 
     and cap between :code:`n_min` and :code:`n_max`. 
   - If :code:`fwer` is not :code:`None`, select genes based on FWER (across all genes of a given trait) and 
     cap between :code:`n_min` and :code:`n_max`. 
   
3. Assign gene weights based on :code:`weight`.
4. Write the :code:`.gs` file to :code:`out_file`.

.. code-block:: bash
    
    # Select top 1,000 genes and use z-score weights
    scdrs munge-gs \
        --out-file <out_file> \
        --zscore-file <zscore_file> \
        --weight zscore \
        --n-max 1000
        
out_file : str
    Output scDRS :code:`.gs` file.
pval_file : str, optional
    P-value file. A .tsv file with first column corresponding to genes and other 
    columns corresponding to p-values of traits (one trait per column). 
    One of :code:`pval-file` and :code:`zscore-file` is expected. Default is :code:`None`. 
zscore_file : str, optional
    Z-score file. A .tsv file with first column corresponding to genes and other 
    columns corresponding to z-scores of traits (one trait per column). 
    One of :code:`pval-file` and :code:`zscore-file` is expected. Default is :code:`None`. 
weight : str, optional
    Gene weight options. One of :code:`zscore` or :code:`uniform`. Default is :code:`zscore`. 
fdr : float, optional
    FDR threshold. Default is :code:`None`. E.g., :code:`--fdr 0.05`
fwer : float, optional
    FWER threshold. Default is :code:`None`. E.g., :code:`--fwer 0.05`
n_min : int, optional
    Minimum number of genes for each gene set. Default is :code:`100`. E.g., :code:`--n-min 100`
n_max : int, optional
    Maximum number of genes for each gene set. Default is :code:`1000`. E.g., :code:`--n-min 1000`
    
Example p-value file::
    
    GENE    BMI    HEIGHT
    OR4F5   0.001  0.01
    DAZ3    0.01   0.001

    
compute-score
~~~~~~~~~~~~~
Compute scDRS scores. Generate :code:`.score.gz` and :code:`.full_score.gz` files for each trait.

.. code-block:: bash

    scdrs compute-score \
        --h5ad-file <h5ad_file>\
        --h5ad-species mouse\
        --gs-file <gs_file>\
        --gs-species human\
        --out-folder <out_folder>\
        --cov-file <cov_file>\
        --flag-filter-data True\
        --flag-raw-count True\
        --n-ctrl 1000\
        --flag-return-ctrl-raw-score False\
        --flag-return-ctrl-norm-score True\
        
h5ad_file : str
    Single-cell :code:`.h5ad` file.
h5ad_species : str
    Species of :code:`h5ad_file`. One of :code:`hsapiens`, :code:`human`, :code:`mmusculus`, or :code:`mouse`.
gs_file : str
    scDRS gene set :code:`.gs` file.
gs_species : str
    Species of :code:`gs_file`. One of :code:`hsapiens`, :code:`human`, :code:`mmusculus`, or :code:`mouse`.
out_folder : str
    Output folder. Save scDRS score files as :code:`<out_folder>/<trait>.score.gz` and 
    scDRS full score files as :code:`<out_folder>/<trait>.full_score.gz`, where trait identifier 
    :code:`<trait>` is from :code:`gs_file` file.
cov_file : str, optional
    scDRS covariate :code:`.cov` file. Default is :code:`None`.
flag_filter_data : bool, optional
    If to apply minimal cell and gene filtering to :code:`h5ad_file`. Default is :code:`True`.
flag_raw_count : bool, optional
    If to apply size-factor normalization and log1p-transformation to :code:`h5ad_file`. Default is :code:`True`.
n_ctrl : int, optional
    Number of control gene sets. Default is :code:`1000`.
flag_return_ctrl_raw_score : bool, optional
    If to return raw control scores. Default is :code:`False`.
flag_return_ctrl_norm_score : bool, optional
    If to return normalized control scores. Default is :code:`True`.


perform-downstream
~~~~~~~~~~~~~~~~~~

Perform scDRS downstream analyses based on precomputed scDRS :code:`.full_score.gz` files.
        
--group-analysis  For a given cell group-level annotation (e.g., tissue or cell type), assess cell 
    group-disease association (control-score-based MC tests using 5% quantile) and within-cell
    group disease-association heterogeneity (control-score-based MC tests using Geary's C).
--corr-analysis  For a given individual cell-level variable (e.g., T cell effectorness gradient),
    assess association between disease and the individual cell-level variable (control-score-based 
    MC tests using Pearson's correlation).
--gene-analysis  Compute Pearson's correlation between expression of each gene and the scDRS disease score. 

.. code-block:: bash

    scdrs perform-downstream \
        --h5ad-file <h5ad_file>\
        --score-file <score_file>\
        --out-folder <out_folder>\
        --group-analysis cell_type \
        --corr-analysis causal_variable,non_causal_variable,covariate\
        --gene-analysis\
        --flag-filter-data True\
        --flag-raw-count True
        
h5ad_file : str
    Single-cell :code:`.h5ad` file.
score_file : str
    scDRS :code:`.full_score.gz` file. Use "@" to specify multiple file names,
    e.g., :code:`<score_folder>/@.full_score.gz`. However, :code:`<score_folder>` 
    should not contain "@".
out_folder : str
    Output folder. 
group_analysis : str, optional
    Comma-seperated column names for cell group annotations in :code:`adata.obs.columns`, e.g., 
    cell types or tissues. Results are saved as :code:`<out_folder>/<trait>.scdrs_group.<annot>`, 
    one file per annotation. Default is :code:`None`. 
corr_analysis : str, optional
    Comma-seperated column names for continuous annotations in :code:`adata.obs.columns`,
    e.g., T cell effectorness gradient. Results are saved as 
    :code:`<out_folder>/<trait>.scdrs_cell_corr` for all variables. 
    Default is :code:`None`.
gene_analysis : str, optional
    Flag to perform the gene prioritization by correlating gene expression with scDRS
    scores. Specifying :code:`--gene-analysis` without any arguments. Results are saved as
    :code:`<out_folder>/<trait>.scdrs_gene` for all genes.  Default is :code:`None`.
flag_filter_data : bool, optional
    If to apply minimal cell and gene filtering to :code:`h5ad_file`. Default is :code:`True`.
flag_raw_count : bool, optional
    If to apply size-factor normalization and log1p-transformation to :code:`h5ad_file`. 
    Default is :code:`True`.
knn_n_neighbors : int, optional
    :code:`n_neighbors` parameter for computing KNN graph using :code:`sc.pp.neighbors`.
    Default is :code:`15` (consistent with the TMS pipeline).
knn_n_pcs : int, optional
    :code:`n_pcs` parameter for computing KNN graph using :code:`sc.pp.neighbors`.
    Default is :code:`20` (consistent with the TMS pipeline).