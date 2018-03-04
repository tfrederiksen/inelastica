.. _setupfcrun:

setupFCrun
==========

Tool for generating a FCrun folder starting from CGrun

usage:
  setupFCrun [-h] [-F FCFIRST] [-L FCLAST] [-d DISPLACEMENT] CG FC

positional arguments:
  CG                    Input CGrun (or TSrun) directory name
  FC                    Output FCrun directory name (to be created)

optional arguments:
  -h, --help            show this help message and exit
  -F FCFIRST, --FCfirst FCFIRST
                        First dynamic atom (default: 1)
  -L FCLAST, --FClast FCLAST
                        Last dynamic atom (default: 1)
  -d DISPLACEMENT, --displacement DISPLACEMENT
                        Finite displacement amplitude (default: 0.02 Ang
