.. _install:

Installation
============

Dependencies
------------

These packages are required:

 * `numpy`_ >= 1.8
 * `scipy`_
 * `netCDF4 <netcdf4-py_>`_

Manual installation
-------------------

With the above packages ready, installation of Inelastica is derformed with the command

.. code-block:: bash

    python setup.py install --prefix=<prefix>
    # or
    python setup.py install --home=<my-python-home>

One may also wish to set the following environment variables

.. code-block:: bash

    export PYTHONPATH=<my-python-home>/lib/python/
    export PATH=$PATH:<my-python-home>/bin/
