.. _computedos:

ComputeDOS
==========

Computing density of states (DOS) from a band calculation on its underlying k-mesh (assuming equal weights to each k-point).

usage:
  ComputeDOS [-h] [-a EMIN] [-b EMAX] [-n PTS] [-s SMEAR] NC XMG

positional arguments:
  NC                    Input netCDF file from SupercellPhonons (containing eigenvalues and k-grid)
  XMG                   Output xmgrace filename (to be created)

optional arguments:
  -h, --help            show this help message and exit
  -a EMIN, --emin EMIN  Energy minimum (default: 0.0 eV)
  -b EMAX, --emax EMAX  Energy maximum (default: 1.0 eV)
  -n PTS, --pts PTS     Points on energy grid (default: 1001 eV)
-s SMEAR, --smear SMEAR
                        Gaussian smearing of eigenvalues (default: 0.0001 eV)
