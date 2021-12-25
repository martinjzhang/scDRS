scDRS
=====

scDRS (single-cell disease-relevance score) is a method for associating individual cells in scRNA-seq data with disease GWASs, built on top of `AnnData <https://anndata.readthedocs.io/en/latest/>`_ and `Scanpy <https://scanpy.readthedocs.io/en/stable/>`_.

Check out the  manuscript `Zhang*, Hou*, et al. "Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data <https://www.biorxiv.org/content/10.1101/2021.09.24.461597v1>`_.

Results for 74 diseases/traits and the TMS FACS data `(cellxgene visualization) <https://scdrs-tms-facs.herokuapp.com/>`_.

Demo for 3 diseases/traits and 3 TMS FACS cell types `(cellxgene visualization) <https://scdrs-demo.herokuapp.com/>`_. 


Installation
============

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS; pip install -e .
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_CLI.py -p no:warnings
   
`Install other versions of scDRS <versions>`_



Usage
=====

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


`Jump to an scDRS demo on Cortex data set. <notebooks/quickstart.html>`_


Citation
========
If scDRS is useful for your research, consider citing:

**Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data**

<author list>

*BioRxiv* 2021. <doi>


.. toctree::
   :maxdepth: 2
   :hidden:

   notebooks/quickstart.ipynb
   reference
   file_format
   versions
   examples
   downloads
   

