.. _geom2geom:

geom2geom
=========

Tool for converting geometries between different files. Currently supports XV, xyz, fdf, ANI, and mkl. Repeat unitcell from XV and fdf. (fdf support limited to Ang, cartesian coord)

usage:
  geom2geom [-h] [-x A1] [-y A2] [-z A3] GIN GOUT

positional arguments:
  GIN             Input geometry file
  GOUT            Output geometry file

optional arguments:
  -h, --help      show this help message and exit
  -x A1, --A1 A1  Repeat unitcell along lattice vector A1
  -y A2, --A2 A2  Repeat unitcell along lattice vector A2
  -z A3, --A3 A3  Repeat unitcell along lattice vector A3
