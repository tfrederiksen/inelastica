.. Inelastica documentation master file, created by
   sphinx-quickstart on Fri Feb 23 13:48:57 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


|license|_
|buildstatus|_
|codacy|_


Welcome to Inelastica's documentation!
======================================

Inelastica is both the name of this whole `Python`_ package as well as that of an included script to compute inelastic transport characteristics.
Inelastica is based on the `SIESTA`_/`TranSIESTA`_ DFT codes.

The latest release can obtained `here <releases_>`_ and the development version through:

.. code-block:: bash

    git clone https://github.com/tfrederiksen/inelastica.git


History
-------
The project was initiated around 2003-2005 by Thomas Frederiksen and Magnus Paulsson while they worked in the group of Mads Brandbyge at the Technical University of Denmark.
Inelastica was originally hosted at `SourceForge <inelastica-old_>`_ but moved to `GitHub <inelastica_>`_ in February 2018.

Features
--------

Inelastica contains a number of scripts such as:

 * **geom2geom**: Geometry conversion between different file formats
 * **Bandstructures**: Computation of electron and phonon band structures
 * **pyTBT**: A Python version of tbtrans for elastic electron transport
 * **EigenChannels**: Eigenchannel analysis and generation of real-space scattering state wave functions
 * **Phonons**: Vibration modes/frequencies and electron-vibration couplings
 * **Inelastica**: Inelastic transport (IETS)


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

   install
   examples
   faq
   more

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Publications

   cite
   publications

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Reference documentation

   api


.. |buildstatus| image:: https://travis-ci.org/tfrederiksen/inelastica.svg
.. _buildstatus: https://travis-ci.org/tfrederiksen/inelastica

.. |license| image:: https://img.shields.io/badge/License-LGPL%20v3-blue.svg
.. _license: https://www.gnu.org/licenses/lgpl-3.0

.. |codacy| image:: https://api.codacy.com/project/badge/Grade/013fe70aa6564ea1bec9df0b3831c834
.. _codacy: https://www.codacy.com/app/brandimarte/inelastica?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tfrederiksen/inelastica&amp;utm_campaign=Badge_Grade
