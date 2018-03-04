.. _geom2zmat:

geom2zmat
=========

Tool for converting geometries to zmatrix FDF format. Currently supports XV, xyz, fdf, ANI, and mkl. (fdf support limited to Ang, cartesian coord)

usage:
  geom2zmat [-h] [-F FIRST] [-L LAST] [-x A1] [-y A2] [-z A3] GIN GOUT

positional arguments:
  GIN                   Input geometry file
  GOUT                  Output zmat-fdf filename

optional arguments:
  -h, --help            show this help message and exit
  -F FIRST, --first FIRST
                        First atom in Z-matrix range (default: 1)
  -L LAST, --last LAST  Last atom in Z-matrix range (default: 1)
  -x A1, --A1 A1        Repeat unitcell along lattice vector A1
  -y A2, --A2 A2        Repeat unitcell along lattice vector A2
  -z A3, --A3 A3        Repeat unitcell along lattice vector A3
