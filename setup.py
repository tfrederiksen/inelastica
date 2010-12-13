import sys, time

# TODO: Compile F90helpers
# TODO: Check prerequireies
# TODO: Testcalculations

def test_prereq():
    # Check for numpy.distutils.
    print "# Testing : numpy f2py."
    try:
        import numpy as N
        import numpy.linalg as LA
    except:
        print "#### ERROR ####"
        print "Inelastica needs the package 'numpy' to run."
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        sys.exit([1])

    try:
        import numpy.distutils
        import numpy.distutils.extension
    except:
        print "#### ERROR ####"
        print "Inelastica requires the f2py extension of numpy."
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        sys.exit([1])

    # Check for ScientificPython including netCDF.
    print "# Testing : ScientificPython."
    try:
        import Scientific.IO.NetCDF as nc
    except:
        print "#### ERROR ####"
        print "Inelastica requires ScientificPython with NetCDF extensions."
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        sys.exit([1])

    print "# Testing : numpy speed."
    # Make sure that numpy is compiled with optimized LAPACK/BLAS
    st=time.time()

    # For release 600!
    a=N.ones((600,600),N.complex)
    #a=N.ones((60,60),N.complex)
    for ii in range(2):
        b=N.dot(a,a)
    c,d = LA.eigh(b)
    en=time.time()
    if en-st>4.0:
        print "#### Warning ####"
        print "numpy was compiled with a slow versions of BLAS/LAPACK."
        print "A minimal test showed that your system takes %3.2f s"%(en-st)
        print "  (normal Xeon5430/ifort/mkl10 takes ~ 1 s)"
        print "Please see http://sourceforge.net/apps/mediawiki/inelastica/"
        tmp = raw_input("Press [enter] to continue.")
 
    print "# Testing passed!"

test_prereq()

from numpy.distutils.core import setup
import numpy.distutils.extension as Next

# Fortran helper files
F90ext = Next.Extension('Inelastica.F90helpers',\
                        ['package/F90/distributegs.F',
                         'package/F90/readNewTSHS.F90',
                         'package/F90/removeUnitCellXij.F',
                         'package/F90/setkpointhelper.F'])

# Main setup of python modules
setup(name='Inelastica',
      version='1.1',
      description='Python tools for SIESTA/TranSIESTA', 
      long_description="""
Provides:
1:	File format conversions for geometry, try: geom2geom --help
2:	Phonon calculations (including e-ph coupling)
3:	Transmission calculations, try: pyTBT --help
4:	Eigenchannels analyzis, try: Eigenchannels --help
5:	IETS calculations, try: Inelastica --help
6:	Scripts to set up the above type of calculations.""",
      author='Magnus Paulsson and Thomas Frederiksen', 
      author_email='magnus.paulsson@hik.se',  
      url='https://sourceforge.net/apps/mediawiki/inelastica', 
      license='GPL. Please cite: Frederiksen et al., PRB 75, 205413 (2007)', 
      package_dir={'Inelastica': 'package'},
      scripts  = ['scripts/Inelastica',
                  'scripts/EigenChannels',
                  'scripts/pyTBT',
                  'scripts/geom2geom',
                  'scripts/geom2zmat',
                  'scripts/SDOS',
                  'scripts/BandStruct',
                  'scripts/siesta_cleanup',
                  'scripts/STMimage',
                  'scripts/Vasp2Siesta',
                  'scripts/bands2xmgr',
                  'scripts/Phonons'
                  ],
      packages=['Inelastica'],
      ext_modules=[F90ext],
      data_files=[('Inelastica/PBS', ['PBS/RUN.OS.pbs','PBS/RUN.py.pbs', \
                           'PBS/RUN.TS.pbs'],)]
      )

