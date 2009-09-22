f2py -c -m F90helpers  distributegs.F setkpointhelper.F removeUnitCellXij.F readNewTSHS.F90 --fcompiler=gnu95
cp F90helpers.so ..
