.. module:: scdrs
.. automodule:: scdrs
   :noindex:

Documentation
=============

Import scDRS as::

   import scdrs

Compute scDRS score
~~~~~~~~~~~~~~~~~~~

.. autosummary:: 
   :toctree: reference

   scdrs.preprocess
   scdrs.score_cell

scDRS downstream analyses
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary:: 
   :toctree: reference

   scdrs.method.downstream_group_analysis
   scdrs.method.downstream_corr_analysis
   scdrs.method.downstream_gene_analysis
   
Data loading
~~~~~~~~~~~~

.. autosummary:: 
   :toctree: reference
    
   scdrs.util.load_h5ad
   scdrs.util.load_scdrs_score
   scdrs.util.load_gs
   scdrs.util.load_homolog_mapping
   
Utils
~~~~~

.. autosummary:: 
   :toctree: reference
    
   scdrs.method.test_gearysc
   scdrs.method.gearys_c
   scdrs.method._pearson_corr
   scdrs.method._pearson_corr_sparse
   scdrs.pp.compute_stats
   scdrs.pp.reg_out
   scdrs.pp._get_mean_var
   scdrs.pp._get_mean_var_implicit_cov_corr
