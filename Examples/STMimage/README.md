STM-image simulation
====================

This example provides an STM-image calculation for a CO molecule adsorbed on a Cu(111) surface,
scanned by a four-atom pyramidal Cu tip.

Note 1
------
STM calculations require a minimum vacuum gap of 5 Ang, i.e., the distance along
the direction of transport between the most protruding substrate atom and the tip-apex atom.

Note 2
------
A relatively large lateral surface plane is in general required, as well as a dense `k` grid,
in order to obtain reliable results. The minimal example included here however qualitatively
reproduces the expected conductance dip over the CO molecule. 

Note 3
------
The calculated images are displayed with the default tip position at its center.

Needs
-----
To run the `STM` script one needs first to obtain the following files from `transiesta`:

        TotalPotential.grid.nc
        Rho.grid.nc

Add the following lines in `RUN.fdf` for the `transiesta` calculation:

        SaveTotalPotential true
        SaveRho true

Partial STM images are saved as `FDcurr[0-999].nc` in respective `k`-point library `[0-999]` in the `./DestDir`.
Wave functions and their derivatives, etc, are saved in the same way as `FD[0-999].nc`.
 
The `k`-averaged STM image is saved in `./DestDir/STMimage.nc`.

Mathematica
-----------
For Mathematica users a Mathematica notebook (`AnalyzeSTMimage.nb`) is included
to visualize the results saved in `STMimage.nc`.
