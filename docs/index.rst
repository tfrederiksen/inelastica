.. Inelastica documentation master file, created by
   sphinx-quickstart on Fri Feb 23 13:48:57 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|license|_
|buildstatus|_


Welcome to Inelastica's documentation!
======================================

Inelastica is a `Python`_ package for the `SIESTA`_/`TranSIESTA`_ DFT codes
as well as a script to compute inelastic transport (and inelastica 
is the corresponding repository name).

Inelastica is hosted here http://github.com/tfrederiksen/inelastica.


Features
--------

Inelastica contains a number of scripts such as:

 * **geom2geom**: Geometry conversion between different file formats
 * **Bandstructures**: Computation of electron and phonon band structures
 * **pyTBT**: A Python version of tbtrans for elastic electron transport
 * **EigenChannels**: Eigenchannel analysis and generation of real-space scattering state wave functions
 * **Phonons**: Vibration modes/frequencies and electron-vibration couplings
 * **Inelastica**: Inelastic transport (IETS)

Installation
------------

Dependencies
~~~~~~~~~~~~

These packages are required:

 * `numpy`_
 * `scipy`_
 * `netCDF4-python <netcdf4-py_>`_

Contributions, issues and bugs
------------------------------

Contributions are highly appreciated.

If you find any bugs please form a `bug report/issue <issue_>`_.

If you have a fix please consider adding a `pull request <pr_>`_.

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. toctree::
   :hidden:
   :maxdepth: 2

   cite
   more



.. |buildstatus| image:: https://travis-ci.org/tfrederiksen/inelastica.svg
.. _buildstatus: https://travis-ci.org/tfrederiksen/inelastica

.. |license| image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
.. _license: https://www.gnu.org/licenses/lgpl-3.0

