Versions
========

v1.0.1: publication (072222)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS 
   git checkout -b pub v1.0.1
   pip install -e .
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_CLI.py -p no:warnings
   
v1.0.0: revision 1 (030622)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS
   git checkout -b rv1 v1.0.0
   pip install -e .
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_CLI.py -p no:warnings

v0.1: initial submission (093022)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS
   git checkout -b initial_submission v0.1 
   pip install -e .
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_scdrs.py -p no:warnings