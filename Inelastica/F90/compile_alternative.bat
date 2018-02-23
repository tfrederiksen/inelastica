f2py -c -m F90helpers expansion_SE.f90 setkpointhelper.f90 removeUnitCellXij.f90 readTSHS.f90 --fcompiler=gfortran
cp F90helpers.so ..

f2py -c -m F90_lapack surfaceGreen.f90 --fcompiler=gfortran  --link-lapack_opt
cp F90_lapack.so ..
