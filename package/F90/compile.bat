f2py -c -m F90helpers distributegs.F setkpointhelper.F removeUnitCellXij.F readNewTSHS.F90 --fcompiler=intelem --compiler=intelem
cp F90helpers.so ..
