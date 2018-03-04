.. _bandstructures:

Bandstructures
==============

Methods to calculate electron and phonon band structures from finite-difference calculations

usage:
  Bandstructures [-h] [--FCwildcard FCWILDCARD] [--OSdir ONLYSDIR] [-r RADIUS] [--AtomicMass ATOMICMASS] [-k KFILE] [-q QFILE] [-s STEPS] [--mesh MESH] [--sort] [--TSdir ONLYTSDIR] [--nbands NBANDS] DestDir

positional arguments:
  DestDir               Destination directory

optional arguments:
  -h, --help            show this help message and exit
  --FCwildcard FCWILDCARD
                        Wildcard for FC directories [default=./FC*]
  --OSdir ONLYSDIR      Location of OnlyS directory [default=./OSrun]
  -r RADIUS, --radius RADIUS
                        Force cutoff radius in Angstroms [default=0.0]
  --AtomicMass ATOMICMASS
                        Option to add to (or override!) existing dictionary of
                        atomic masses. Format is a list
                        [[anr1,mass1(,label)],...] [default=[]]
  -k KFILE, --kpointfile KFILE
                        Input file with electronic k-points to be evaluated
                        [default=None]
  -q QFILE, --qpointfile QFILE
                        Input file with phonon q-points to be evaluated
                        [default=None]
  -s STEPS, --steps STEPS
                        Number of points on path between high-symmetry
                        k-points [default=100]
  --mesh MESH           Mesh sampling over one BZ (powers of 2)
                        [default=[0,0,0]]
  --sort                Sort eigenvalues along k-mesh for nice plots?
                        [default=False]
  --TSdir ONLYTSDIR     Location of TranSIESTA calculation directory (will
                        ignore FC and OnlyS directories) [default=None]
  --nbands NBANDS       Number of electronic bands to be included in netCDF
                        output (lower energy bands) [default=None]
