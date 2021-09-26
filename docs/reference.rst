.. module:: scdrs
.. automodule:: scdrs
   :noindex:

Documentation
=============

Import scDRS as::

   import scdrs

Individual cell disease-relevance score
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary:: 
   :toctree: reference

   scdrs.method.score_cell
   scdrs.method.compute_stats

Downstream analyses
~~~~~~~~~~~~~~~~~~~

.. autosummary:: 
   :toctree: reference

   scdrs.util.test_gearysc
   scdrs.util.convert_gs_species

.. Data loading
.. ~~~~~~~~~~~~
.. utility function for downloading data, we can have a function which load commonly used GWAS datasets