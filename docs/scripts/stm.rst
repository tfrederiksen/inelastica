.. _stm:

STM
===

Script that calculates STM images using the Bardeen approximation outlined in
PRB 93 115434 (2016) and PRB 96 085415 (2017). The script is divided into 3
parts:

1. Calculation of the scattering states at the Fermi-energy on the same
real space grid as TranSiesta (real-space cutoff). These are saved in
DestDir/SystemLabel.A[LR][0-99].nc files and are reused if found. NEEDS:
TranSiesta calculation.

2. Propagation of the scattering states from a
surface (defined by a constant charge density) out into the vacuum region.
After the x-y plane, where the average potential of the slice is maximum (the
separation plane), is found, the potential is ascribed a constant value at
this average. Saves the propagated wavefunctions at the separation plane in
DestDir/[kpoint]/FD[kpoint].nc. NEEDS: TotalPotential.grid.nc and Rho.grid.nc.

3. Conductance calculation where the tip/substrate wavefunctions are displaced
to simulate the conductance at different tip-positions. The k averaged STM
image and the STM images of individual k points are saved in
DestDir/STMimage.nc.

Usage:
  STM [options] DestinationDirectory

Options:
  -h, --help            show this help message and exit
  -F DEVICEFIRST, --DeviceFirst=DEVICEFIRST
                        First device atom (SIESTA numbering) [TS.TBT.PDOSFrom]
  -L DEVICELAST, --DeviceLast=DEVICELAST
                        Last device atom (SIESTA numbering) [TS.TBT.PDOSTo]
  -e ENERGY, --Energy=ENERGY
                        Energy where scattering states are evaluated [0.0 eV]
  --eta=ETA             Imaginary part added to all energies (device and
                        leads) [1e-06 eV]
  -l ETALEAD, --etaLead=ETALEAD
                        Additional imaginary part added ONLY in the leads
                        (surface GF) [0.0 eV]
  -f FN, --fdf=FN       Input fdf-file for TranSIESTA calculations [./RUN.fdf]
  -s ISPIN, --iSpin=ISPIN
                        Spin channel [0]
  -p, --savePOS         Save the individual solutions as .pos files
  --shift               Shift current 1/2 cell in x, y directions
  --bulk                Use bulk in electrodes. The Hamiltonian from the
                        electrode calculation is inserted into the electrode
                        region in the TranSIESTA cell [TS.UseBulkInElectrodes]
  --nobulk              Use only self-energies in the electrodes. The full
                        Hamiltonian of the TranSIESTA cell is used in
                        combination with self-energies for the electrodes
                        [TS.UseBulkInElectrodes]
  --scaleSigL=SCALESIGL
                        Scale factor applied to Sigma_L [default=1.0]
  --scaleSigR=SCALESIGR
                        Scale factor applied to Sigma_R [default=1.0]
  -u, --useSigNC        Use SigNCfiles [False]
  --SpectralCutoff=SPECTRALCUTOFF
                        Cutoff value for SpectralMatrix functions (for
                        ordinary matrix representation set cutoff<=0.0)
                        [default=0.0]
  -n NCPU, --nCPU=NCPU  Number of processors [1]
  -x NK1, --Nk1=NK1     k-points Nk1 along a1 [1]
  -y NK2, --Nk2=NK2     k-points Nk2 along a2 [1]
  -r RHOISO, --rhoiso=RHOISO
                        Density at the isosurface from which the localized-
                        basis wave functions are propagated [default=0.001
                        Bohr^-3Ry^-1]
  --ssp=SHIFTSEPARATIONPLANE
                        Manually shift the separation plane (>0 means away
                        from the substrate) [default=0 Ang]
  --sc=SAMPLINGSCALE    Sampling scale of wave functions in lateral plane. 1
                        means same real-space resolution as used in
                        TranSiesta, while 2 means doubling of the lateral
                        lattice constants, etc. [default=2]
  -t TOLSOLVE, --tolsolve=TOLSOLVE
                        Tolerance for the iterative linear solver
                        (scipy.linalg.isolve.gmres) [default=1e-06]
  --savelocwfs          Save localized-basis wave functions.
