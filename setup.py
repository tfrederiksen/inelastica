import sys, time

# TODO: Compile F90helpers
# TODO: Check prerequireies
# TODO: Testcalculations

def test_prereq():
    # Check for numpy.distutils.
    print "Testing : numpy f2py."
    try:
        import numpy as N
        import numpy.linalg as LA
    except:
        print "#### ERROR ####"
        print "Inelastica needs the package 'numpy' to run."
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        sys.exit(1)

    try:
        import numpy.distutils
        import numpy.distutils.extension
    except:
        print "#### ERROR ####"
        print "Inelastica requires the f2py extension of numpy."
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        sys.exit(1)

    # Check for ScientificPython including netCDF.
    print "Testing : ScientificPython"
    try:
        import Scientific.IO.NetCDF as nc
    except:
        print "#### ERROR ####"
        print "Inelastica requires ScientificPython with NetCDF extensions."
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        sys.exit(1)

    print "Testing : numpy speed."
    # Make sure that numpy is compiled with optimized LAPACK/BLAS
    st = time.time()

    # For release 600!
    a = N.ones((600,600),N.complex)
    b = N.dot(a,a)
    c,d = LA.eigh(b)
    en = time.time()
    print "A minimal test showed that your system takes %3.2f s"%(en-st)
    if en-st>4.0:
        print "#### Warning ####"
        print "numpy was compiled with a slow versions of BLAS/LAPACK."
        print "  (normal Xeon5430/ifort/mkl10 takes ~ 1 s)"
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        tmp = raw_input("Press [enter] to continue.")

    print "Testing passed!"

test_prereq()

from numpy.distutils.core import setup
from numpy.distutils.system_info import get_info, NotFoundError
import numpy.distutils.extension as Next

# Fortran helper files
F90ext = Next.Extension('Inelastica.F90helpers',\
                            ['package/F90/expansion_SE.f90',
                             'package/F90/readTSHS.f90',
                             'package/F90/removeUnitCellXij.f90',
                             'package/F90/setkpointhelper.f90'],
                        )

# Retrieve the LAPACK-library...
lapack_opt = get_info('lapack_opt')
if not lapack_opt:
    raise NotFoundError('No LAPACK/BLAS resources found')
F90extLapack = Next.Extension('Inelastica.F90_lapack',\
                                  ['package/F90/surfaceGreen.f90'],
                              **lapack_opt)

# Main setup of python modules
setup(name='Inelastica',
      version='1.2',
      # Define the requirements for Inelastica
      # These probably needs to be adjusted... 
      install_requires = ['python>=2.5','numpy>=1.6','ScientificPython>=2.6'],
      requires = ['python (>=2.5)','numpy (>=1.6)','ScientificPython (>=2.6)'],
      description='Python tools for SIESTA/TranSIESTA', 
      long_description="""
Provides:
1:	File format conversions for geometry, try: geom2geom --help
2:	Phonon calculations (including e-ph coupling)
3:	Transmission calculations, try: pyTBT --help
4:	Eigenchannels analysis, try: EigenChannels --help
5:	IETS calculations, try: Inelastica --help
6:      grid2grid, try: grid2grid --help
7:	Scripts to set up the above type of calculations.""",
      author='Magnus Paulsson and Thomas Frederiksen', 
      author_email='magnus.paulsson@lnu.se / thomas_frederiksen@ehu.es',  
      url='https://sourceforge.net/apps/mediawiki/inelastica', 
      license='GPL. Please cite: Frederiksen et al., PRB 75, 205413 (2007)', 
      package_dir={'Inelastica': 'package'},
      scripts  = ['scripts/Inelastica',
                  'scripts/EigenChannels',
                  'scripts/pyTBT',
                  'scripts/geom2geom',
                  'scripts/geom2zmat',
                  'scripts/BandStruct',
                  'scripts/siesta_cleanup',
                  'scripts/STMimage',
                  'scripts/Vasp2Siesta',
                  'scripts/bands2xmgr',
                  'scripts/Phonons',
                  'scripts/NEB',
                  'scripts/grid2grid'
                  ],
      packages=['Inelastica'],
      ext_modules=[F90ext,F90extLapack],
      data_files=[('Inelastica/PBS', ['PBS/RUN.OS.pbs','PBS/RUN.py.pbs', \
                           'PBS/RUN.TS.pbs'],)]
      )

