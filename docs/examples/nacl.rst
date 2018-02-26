.. _nacl:

NaCl
----

Phonon band structure for NaCl (see `Examples/Phonon_bands/NaCl`).

Lattice parameter a = 5.64 Å and a 6x6x6 repetition of the primitive cell. Main computational settings:

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

.. image:: results/NaCl_Electrons.png
   :scale: 90 %
   :alt: electronic band structure
   :align: center

**Desity of states**

Density of states (DOS) sampled on a grid of 64x64x64 k-points:

.. image:: results/ElectronDOS-NaCl.png
   :scale: 90 %
   :alt: electronic density of states
   :align: center


Phonons
~~~~~~~

**Band structure**

Phonon band structure computed with different force cutoff radii r = 5.0, 7.0, 9.0 Å:

.. image:: results/NaCl_r5.0_Phonons.png
   :scale: 80 %
   :alt: phonon band structure with force cutoff radii of 5.0 angstroms
   :align: center
.. image:: results/NaCl_r7.0_Phonons.png
   :scale: 80 %
   :alt: phonon band structure with force cutoff radii of 7.0 angstroms
   :align: center
.. image:: results/NaCl_r9.0_Phonons.png
   :scale: 80 %
   :alt: phonon band structure with force cutoff radii of 9.0 angstroms
   :align: center

**Desity of states**

Density of states (DOS) sampled on a grid of 64x64x64 k-points:

.. image:: results/PhononDOS-NaCl.png
   :scale: 90 %
   :alt: phonon density of states 
   :align: center


Reference results
~~~~~~~~~~~~~~~~~

One important difference is the modes at *Gamma*. With `Inelastica <docs_>`_ the optical modes are found to be degenerate. The missing LO-TO splitting is likely the effect of macroscopic polarization in polar materials. The long-range electrostatic interactions are not corrected in the current code revisions.

.. image:: results/Raunio_1970.png
   :scale: 90 %
   :alt: reference phonon density of statesband structure 
   :align: center

Raunio, G. & Rolandson, S. `Lattice Dynamics of NaCl, KCl, RbCl, and RbF <http://journals.aps.org/prb/abstract/10.1103/PhysRevB.2.2098>`_, Phys. Rev. B **2**, 2098 (1970).
