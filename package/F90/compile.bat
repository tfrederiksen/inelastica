f2py -c -m F90helpers distributegs.f90 setkpointhelper.f90 removeUnitCellXij.f90 readTSHS.f90 surfaceGreen.f90 --fcompiler=intelem --compiler=intelem
cp F90helpers.so ..
