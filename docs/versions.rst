Versions
========

1.0.0: current version
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS; pip install -e .
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_CLI.py -p no:warnings

beta: initial submission version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/martinjzhang/scDRS.git
   cd scDRS; pip install -e .
   git checkout -b initial_submission v0.1 
   
Quick test:

.. code-block:: bash

   python -m pytest tests/test_scdrs.py -p no:warnings