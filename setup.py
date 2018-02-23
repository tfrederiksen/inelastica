#!/usr/bin/env python

import sys, time

# TODO: Compile F90helpers
# TODO: Check prerequireies
# TODO: Testcalculations

def test_prereq():

    try:
        import numpy as N
        import numpy.linalg as LA
    except:
        print "#### ERROR ####"
        print "Inelastica needs the package 'numpy' to run."
        print "#### ERROR ####"
        raise NameError('numpy package not found')

    try:
        import numpy.distutils
        import numpy.distutils.extension
    except:
        print "#### ERROR ####"
        print "Inelastica requires the f2py extension of numpy."
        print "#### ERROR ####"
        raise NameError('numpy f2py package not found')

    try:
        import netCDF4 as NC4
    except:
        print "#### ERROR ####"
        print "Inelastica requires netCDF4 (1.2.7 or newer recommended)"
        print "https://pypi.python.org/pypi/netCDF4"
        print "#### ERROR ####"
        raise NameError('netCDF4 package not found')

    # Make sure that numpy is compiled with optimized LAPACK/BLAS
    st = time.time()

    # For release 600!
    a = N.ones((600,600),N.complex)
    b = N.dot(a,a)
    c,d = LA.eigh(b)
    en = time.time()
    if en - st > 4.0:
        print "#### Warning ####"
        print "A minimal test showed that your system takes %3.2f s"%(en-st)
        print "numpy was compiled with a slow versions of BLAS/LAPACK."
        print "  (normal Xeon5430/ifort/mkl10 takes ~ 1 s)"
        print "Please see http://dipc.ehu.es/frederiksen/inelastica/index.php"
        print "#### Warning ####"

    try:
        import scipy
        import scipy.linalg as SLA
        import scipy.special as SS
    except:
        print "#### Warning ####"
        print 'Some modules will not work without the scipy package'
        print '(needed for solving generalized eigenvalue problems'
        print 'and spherical harmonics)'
        print "#### Warning ####"

test_prereq()

from numpy.distutils.core import setup
from numpy.distutils.system_info import get_info, NotFoundError
import numpy.distutils.extension as Next

# Fortran helper files
F90ext = Next.Extension('Inelastica.F90helpers',
                        ['Inelastica/fortran/expansion_SE.f90',
                         'Inelastica/fortran/readTSHS.f90',
                         'Inelastica/fortran/removeUnitCellXij.f90',
                         'Inelastica/fortran/setkpointhelper.f90'],
                        )

# Retrieve the LAPACK-library...
lapack_opt = get_info('lapack_opt')
if not lapack_opt:
    raise NotFoundError('No LAPACK/BLAS resources found')
F90extLapack = Next.Extension('Inelastica.F90_lapack',
                              ['Inelastica/fortran/surfaceGreen.f90'],
                              **lapack_opt)

# Main setup of python modules
setup(name='Inelastica',
      version='1.2-rc',
      # Define the requirements for Inelastica
      # These probably needs to be adjusted... 
      requires = ['python (>=2.7)','numpy (>=1.8)','netCDF4 (>=1.2.7)'],
      description='Python tools for SIESTA/TranSIESTA', 
      author='Magnus Paulsson and Thomas Frederiksen', 
      author_email='magnus.paulsson@lnu.se / thomas_frederiksen@ehu.es',  
      url='https://github.com/tfrederiksen/inelastica', 
      license='GPL', 
      package_dir={'Inelastica': 'Inelastica'},
      scripts  = ['Inelastica/scripts/agr2pdf',
                  'Inelastica/scripts/Inelastica',
                  'Inelastica/scripts/EigenChannels',
                  'Inelastica/scripts/pyTBT',
                  'Inelastica/scripts/geom2geom',
                  'Inelastica/scripts/geom2zmat',
                  'Inelastica/scripts/Bandstructures',
                  'Inelastica/scripts/ComputeDOS',
                  'Inelastica/scripts/siesta_cleanup',
                  'Inelastica/scripts/Vasp2Siesta',
                  'Inelastica/scripts/bands2xmgr',
                  'Inelastica/scripts/Phonons',
                  'Inelastica/scripts/NEB',
                  'Inelastica/scripts/grid2grid',
                  'Inelastica/scripts/setupFCrun',
                  'Inelastica/scripts/setupOSrun',
                  'Inelastica/scripts/kaverage-TBT',
                  'Inelastica/scripts/STM',
                  'Inelastica/scripts/kaverage-IETS',
                  'Inelastica/scripts/average-gridfunc',
                  'Inelastica/scripts/WriteWavefunctions'
                  ],
      packages=['Inelastica'],
      ext_modules=[F90ext,F90extLapack],
      )

