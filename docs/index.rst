scDRS
=====

scDRS (single-cell disease-relevance score) is a method for associating individual cells in scRNA-seq data with disease GWASs, built on top of `AnnData <https://anndata.readthedocs.io/en/latest/>`_ and `Scanpy <https://scanpy.readthedocs.io/en/stable/>`_.

Check out the bioRxiv manuscript `Zhang*, Hou*, et al. "Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data <XXX>`_.

Results for 74 diseases/traits and the TMS FACS data `(cellxgene visualization) <https://scdrs-tms-facs.herokuapp.com/>`_.

Demo for 3 diseases/traits and 3 TMS FACS cell types `(cellxgene visualization) <https://scdrs-demo.herokuapp.com/>`_. 


Installation
============

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS; pip install -e .


Usage
=====

.. code-block:: python

   import os
   import pandas as pd
   from anndata import read_h5ad
   import scdrs

   # Load data
   DATA_PATH = scdrs.__path__[0]
   adata = read_h5ad(os.path.join(DATA_PATH, "data/toydata_mouse.h5ad"))
   df_gs = pd.read_csv(os.path.join(DATA_PATH, "data/toydata_mouse.gs"), sep="\t")

   # Compute scDRS gene-level and cell-level statistics
   scdrs.method.compute_stats(adata)

   # Compute scDRS results
   gene_list = df_gs["GENESET"].values[0].split(",")
   df_res = scdrs.method.score_cell(adata, gene_list)
   print(df_res.iloc[:4])

`Jump to an scDRS demo on Cortex data set. <notebooks/quickstart.html>`_


Citation
========
If you use scDRS for published work, please cite our manuscript:

**Polygenic enrichment distinguishes disease associations of individual cells in single-cell RNA-seq data**

<author list>

*BioRxiv* 2021. <doi>


.. toctree::
   :maxdepth: 2
   :hidden:

   notebooks/quickstart.ipynb
   reference
   examples
   downloads
   

