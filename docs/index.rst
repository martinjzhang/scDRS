scDRS
=====
scDRS is a Python package for evaluating disease association at each single cell.

.. TODO: texts describing features of scDRS

.. TODO: vignettes demonstraing functions of scDRS

Results for 74 diseases/traits and the TMS FACS data `(cellxgene visualization) <https://scdrs-tms-facs.herokuapp.com/>`_.

Demo for 3 diseases/traits and 3 TMS FACS cell types `(cellxgene visualization) <https://scdrs-demo.herokuapp.com/>`_. 


Installation
============

.. code-block:: bash

   git clone git@github.com:martinjzhang/scDRS.git & cd scDRS
   pip install -r requirements.txt; pip install -e .

Usage
=====

.. code-block:: bash

   # Script for batch processing scores.
   compute_score.py 
   

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
   

