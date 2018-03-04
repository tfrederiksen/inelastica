.. _average-gridfunc:

average-gridfunc
================

Tool for 2D averages of SIESTA 3D gridfunctions

usage:
  average-gridfunc [-h] [-a AXIS] [-s SPIN] [-B] ncfile [ncfile ...]

positional arguments:
  ncfile                Input netcdf gridfunction file

optional arguments:
  -h, --help            show this help message and exit
  -a AXIS, --axis AXIS  Axis along which the transverse 2D average should be
                        performed [default: 3]
  -s SPIN, --spin SPIN  Spin index [default: 0]
  -B, --Bohr            Bohr units [default: False]
