.. _kaverage-tbt:

kaverage-TBT
============

Utility to perform k-averages on relevant quantities in netCDF4 files from TBTrans.

usage:
  kaverage-TBT [-h] [-w] [-v VAR [VAR ...]] [-t TERM [TERM ...]] [-s SKIP [SKIP ...]] [-e EXT] NC [NC ...]

positional arguments:
  NC                    A netCDF4 file from TBTrans

optional arguments:
  -h, --help            show this help message and exit
  -w                    Write k-averaged quantities also in ascii format
  -v VAR [VAR ...], --var VAR [VAR ...]
                        Variable names for k-averaging
  -t TERM [TERM ...], --term TERM [TERM ...]
                        String terminations of variables for k-averaging
  -s SKIP [SKIP ...], --skipvar SKIP [SKIP ...]
                        Skip variables matching this list [default:
  -e EXT, --ext EXT     Output filename termination
