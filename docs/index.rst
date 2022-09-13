scDRS
=====

scDRS (single-cell disease-relevance score) is a method for associating individual cells in scRNA-seq data with disease GWASs, built on top of `AnnData <https://anndata.readthedocs.io/en/latest/>`_ and `Scanpy <https://scanpy.readthedocs.io/en/stable/>`_.

Check out our manuscript `Zhang*, Hou*, et al. "Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data <https://www.nature.com/articles/s41588-022-01167-z>`_.

Explore results for 74 diseases/traits and the TMS FACS data on `cellxgene <https://scdrs-tms-facs.herokuapp.com/>`_.


Installation
============

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS
   git checkout -b v102 v1.0.2
   pip install -e .
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_CLI.py -p no:warnings
   
Install via `PyPI <https://pypi.org/project/scdrs/1.0.1/#description>`_

.. code-block:: bash

    pip install scdrs==1.0.2
    
Quick test for PyPI installation: open Python (>=3.5) and run the code in the Usage section below.
   
`Install other versions <versions.html>`_


Usage
=====

Use `scDRS command-line interface (CLI) <reference_cli.html>`_ for standard analyses.

Use `scDRS Python API <reference.html>`_ for customized analyses. 

Here is a toy example for computing scDRS scores. 

.. code-block:: python

    import os
    import pandas as pd
    import scdrs

    DATA_PATH = scdrs.__path__[0]
    H5AD_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.h5ad")
    COV_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.cov")
    GS_FILE = os.path.join(DATA_PATH, "data/toydata_mouse.gs")

    # Load .h5ad file, .cov file, and .gs file
    adata = scdrs.util.load_h5ad(H5AD_FILE, flag_filter_data=False, flag_raw_count=False)
    df_cov = pd.read_csv(COV_FILE, sep="\t", index_col=0)
    df_gs = scdrs.util.load_gs(GS_FILE)

    # Preproecssing .h5ad data compute scDRS score
    scdrs.preprocess(adata, cov=df_cov)
    gene_list = df_gs['toydata_gs_mouse'][0]
    gene_weight = df_gs['toydata_gs_mouse'][1]
    df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=20)

    print(df_res.iloc[:4])
    
Expected results:

.. csv-table::
   :header: "index", "raw_score", "norm_score", "mc_pval", "pval", "nlog10_pval", "zscore"
   
   N1.MAA000586.3_8_M.1.1-1-1            , 4.741197 , 6.3260064 , 0.04761905, 0.0016638935, 2.7788744   , 2.9357162
   F10.D041911.3_8_M.1.1-1-1             , 4.739066 , 5.916272  , 0.04761905, 0.0016638935, 2.7788744   , 2.9357162
   A17_B002755_B007347_S17.mm10-plus-7-0 , 4.6366262, 5.5523157 , 0.04761905, 0.0016638935, 2.7788744   , 2.9357162
   C22_B003856_S298_L004.mus-2-0-1       , 4.6805663, 7.2986684 , 0.04761905, 0.0016638935, 2.7788744   , 2.9357162
   G12.B002765.3_38_F.1.1-1-1            , 4.640043 , 5.7792473 , 0.04761905, 0.0016638935, 2.7788744   , 2.9357162
   H5.B003278.3_38_F.1.1-1-1             , 4.4457436, -0.5613674, 0.7619048 , 0.687188    , 0.16292442  , -0.48789537
   O14.MAA000570.3_8_M.1.1-1-1           , 4.4552336, -1.5821338, 0.95238096, 0.9467554   , 0.023762206 , -1.6141763
   J21.B000634.3_56_F.1.1-1-1            , 4.4433637, -2.3119287, 1.0       , 0.9916805   , 0.0036282123, -2.3945906
   E5.B002765.3_38_F.1.1-1-1             , 4.4870768, 1.1566308 , 0.23809524, 0.13311148  , 0.87578446  , 1.1118028
   K20_B000268_B009896_S260.mm10-plus-4-0, 4.53548  , -3.1656132, 1.0       , 1.0         , -0.0        , -10.0

Examples
========
- `Tutorial on a mouse Cortex data set. <notebooks/quickstart.html>`_
- Coming soon


.. Citation
   ========
   If scDRS is useful for your research, consider citing:
   
   **Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data**
   
   <author list>

   *BioRxiv* 2021. <doi>


.. toctree::
   :maxdepth: 2
   :hidden:

   reference_cli
   reference
   file_format
   faq
   versions
   notebooks/quickstart.ipynb
   downloads
   

