.. _au_fcc:

Au FCC
------

Phonon band structure for Au (see `Examples/Phonon_bands/Au_FCC`).

Lattice parameter a = 4.08 Å and a 6x6x6 repetition of the primitive cell. Main computational settings:

.. code-block:: bash

    PAO.EnergyShift 0.01 Ry
    PAO.BasisSize   DZP
    XC.functional   GGA
    XC.authors      PBE
    MeshCutoff      500. Ry
    MD.FCDispl      0.02 Ang


Electrons
~~~~~~~~~

**Band structure**

.. image:: results/Au_FCC_Electrons.png
   :scale: 90 %
   :alt: electronic band structure
   :align: center

**Desity of states**

Density of states (DOS) sampled on a grid of 64x64x64 k-points:

.. image:: results/ElectronDOS-AuFCC.png
   :scale: 90 %
   :alt: electronic density of states
   :align: center


Phonons
~~~~~~~

**Band structure**

Phonon band structure computed with different force cutoff radii r = 5.0, 7.0, 9.0 Å:

.. image:: results/AuFCC_r5.0_Phonons.png
   :scale: 80 %
   :alt: phonon band structure with force cutoff radii of 5.0 angstroms
   :align: center
.. image:: results/AuFCC_r7.0_Phonons.png
   :scale: 80 %
   :alt: phonon band structure with force cutoff radii of 7.0 angstroms
   :align: center
.. image:: results/AuFCC_r9.0_Phonons.png
   :scale: 80 %
   :alt: phonon band structure with force cutoff radii of 9.0 angstroms
   :align: center

**Desity of states**

Density of states (DOS) sampled on a grid of 64x64x64 k-points:

.. image:: results/PhononDOS-AuFCC.png
   :scale: 90 %
   :alt: phonon density of states 
   :align: center


Reference results
~~~~~~~~~~~~~~~~~

Lynn, J. W.; Smith, H. G. & Nicklow, R. M. `Lattice-dynamics of gold <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.8.3493>`_, Phys. Rev B **8**, 3493-3499 (1973).

.. image:: results/Lynn_1973.png
   :scale: 90 %
   :alt: reference phonon density of statesband structure 
   :align: center
