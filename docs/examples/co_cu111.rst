.. _co_cu111:

CO/Cu(111)
----------

In this example we describe how to obtain STM images of a CO molecule adsorbed on a Cu(111) surface consisting of 3x3 Cu atoms with lattice constant 2.57 Å, using 5x5 k points along the surface plane.

.. code-block:: bash

   #TranSiesta:
   XC.functional   GGA
   XC.authors      PBE
   MeshCutoff      200. Ry        
   PAO.EnergyShift 0.001 Ry       
   PAO.BasisSize   SZP#DZP for C and O

   #STM command:
   STM -F 64 -L 87 -n 1 -x 5 -y 5 STM5x5

STM image
~~~~~~

**k averaged STM image**

The k-averaged STM image and its cross-section:

.. image:: results/STMimage_CO_at_Cu111.png
   :scale: 80 %
   :alt: 
   :align: center

**Images of individual k points**

STM images of individual k points along the surface plane:

.. image:: results/STMimage_CO_at_Cu111_kspace.png
   :scale: 80 %
   :alt: 
   :align: center


Convergence in k space
~~~~~~

Cross-section profiles using different number of k points in the surface plane:

.. image:: results/STMimage_CO_at_Cu111_conv.png
   :scale: 80 %
   :alt: 
   :align: center
	  
