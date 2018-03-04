.. _phonons:

Phonons
=======

Methods to calculate vibrations and e-ph couplings from SIESTA output

usage:
  Phonons [-h] [-c] [-r] [--CheckPointNetCDF CHECKPOINTNETCDF] [-s] [-F DEVICEFIRST] [-L DEVICELAST] [--FCfirst FCFIRST] [--FClast FCLAST] [--EPHfirst EPHFIRST] [--EPHlast EPHLAST] [--PBCFirst PBCFIRST] [--PBCLast PBCLAST] [--FCwildcard FCWILDCARD] [--OSdir ONLYSDIR] [-a] [-i ISOTOPES] [-x K1] [-y K2] [-z K3] [-g] DestDir

positional arguments:
  DestDir               Destination directory

optional arguments:
  -h, --help            show this help message and exit
  -c, --CalcCoupl       Calculate e-ph couplings [default=False]
  -r, --Restart         Restart from a previous run [default=False]
  --CheckPointNetCDF CHECKPOINTNETCDF
                        Old NetCDF file used for restart [default=None]
  -s, --SinglePrec      Calculate e-ph couplings using single precision arrays
                        [default=False]
  -F DEVICEFIRST, --DeviceFirst DEVICEFIRST
                        First device atom index (in the electronic basis)
                        [default=1]
  -L DEVICELAST, --DeviceLast DEVICELAST
                        Last device atom index (in the electronic basis)
                        [default=1000]
  --FCfirst FCFIRST     First FC atom index [default=1]
  --FClast FCLAST       Last FC atom index [default=1000]
  --EPHfirst EPHFIRST   First atom index for which the e-ph. couplings are
                        evaluated [default=FCfirst]
  --EPHlast EPHLAST     Last atom index for which the e-ph. couplings are
                        evaluated [default=FClast]
  --PBCFirst PBCFIRST   For eliminating interactions through periodic boundary
                        conditions in z-direction [default=1]
  --PBCLast PBCLAST     For eliminating interactions through periodic boundary
                        conditions in z-direction [default=1000]
  --FCwildcard FCWILDCARD
                        Wildcard for FC directories [default=./FC*]
  --OSdir ONLYSDIR      Location of OnlyS directory [default=./OSrun]
  -a, --AbsoluteEnergyReference
                        Use an absolute energy reference (Fermi energy of
                        equilibrium structure) for displaced Hamiltonians
                        (e.g., when eF is not well-defined) instead of the
                        instantaneous Fermi energy for the displaced
                        geometries, cf. Eq.(17) in PRB 75, 205413 (2007)
                        [default=False]
  -i ISOTOPES, --Isotopes ISOTOPES
                        String, formatted as a list [[i1,m1],...], where the
                        mass of atom index i1 (SIESTA numbering) will be set
                        to m1. Alternatively, the argument can be a file with
                        the string [default=[]]
  -x K1, --k1 K1        k-point along a1 where e-ph couplings are evaluated
                        [0.0]
  -y K2, --k2 K2        k-point along a2 where e-ph couplings are evaluated
                        [0.0]
  -z K3, --k3 K3        k-point along a3 where e-ph couplings are evaluated
                        [0.0]
  -g, --WriteGradients  Write real-space gradients dH/dR to NetCDF
                        [default=False]
