f2py -c -m F90helpers distributegs.f90 setkpointhelper.f90 removeUnitCellXij.f90 readTSHS.f90 --fcompiler=intelem --compiler=intelem
cp F90helpers.so ..

f2py -c -m F90_lapack surfaceGreen.f90 --fcompiler=intelem --compiler=intelem
cp F90_lapack.so ..
