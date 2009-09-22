import sys, profile
sys.path+=['..']

import pyTBT
import SiestaIO as SIO
import numpy as N
import numpy.random as RA

def main():
    if not pyTBT.F90imported:
        print "To test the F90 routines you better compile them first"
        kuk

    # Test readNewTSHS routines
    
    HS1=SIO.HS('../TestCalculations/Self-energy-FCC111/ELEC-1x1//Au3D_BCA.TSHS',UseF90helpers=False)
    HS2=SIO.HS('../TestCalculations/Self-energy-FCC111/ELEC-1x1//Au3D_BCA.TSHS',UseF90helpers=True)
    
    # Check that we read the same TSHS data
    if HS1.gamma!=HS2.gamma or HS1.onlyS!=HS2.onlyS or HS1.nuo!=HS2.nuo or HS1.no!=HS2.no or HS1.nspin!=HS2.nspin or HS1.maxnh!=HS2.maxnh or HS1.qtot!=HS2.qtot or HS1.temp!=HS2.temp or HS1.nua!=HS2.nua or HS1.ef!=HS2.ef or HS1.istep!=HS2.istep or HS1.ia1!=HS2.ia1:
        print "readNesTSHS failed to give same results between python and fortran code for simple variables!!!"
        kuk
    if N.sum(N.abs(HS1.cell-HS2.cell))+N.sum(N.abs(HS1.lasto-HS2.lasto))+N.sum(N.abs(HS1.numh-HS2.numh))+N.sum(N.abs(HS1.listh-HS2.listh))+N.sum(N.abs(HS1.indxuo-HS2.indxuo))+N.sum(N.abs(HS1.xa-HS2.xa))+N.sum(N.abs(HS1.isa-HS2.isa))+N.sum(N.abs(HS1.Ssparse-HS2.Ssparse))+N.sum(N.abs(HS1.Hsparse-HS2.Hsparse))+N.sum(N.abs(HS1.xij-HS2.xij))>1e-9:
        print "readNesTSHS failed to give same results between python and fortran code for arrays!!!"
        kuk

    # Test removeUnitCellXij
    elec1=pyTBT.surfaceGF('../TestCalculations/Self-energy-FCC111/ELEC-1x1//Au3D_BCA.TSHS',3,3,UseF90helpers=False)
    elec2=pyTBT.surfaceGF('../TestCalculations/Self-energy-FCC111/ELEC-1x1//Au3D_BCA.TSHS',3,3,UseF90helpers=True)
    maxerr=N.max(abs(elec1.HS.xij-elec2.HS.xij))
    print "Maximum difference between Xij :",maxerr    

    elec=pyTBT.surfaceGF('../TestCalculations/Self-energy-FCC111/ELEC-1x1//Au3D_BCA.TSHS',3,3)

    for ii in range(10):
        k=N.array(RA.random(3),N.float)
        print " Checking k-point: %f,%f,%f"%(k[0],k[1],k[2])
        elec.HS.kpoint=N.zeros((3,),N.float) # To ensure it is calculated!
        elec.HS.setkpoint(k,UseF90helpers=True)
        H1, S1 = elec.HS.H.copy(), elec.HS.S.copy()
        elec.HS.kpoint=N.zeros((3,),N.float) # To ensure it is calculated!
        elec.HS.setkpoint(k,UseF90helpers=False)
        H2, S2 = elec.HS.H.copy(), elec.HS.S.copy()
        tmp1=N.max(abs(H1-H2))
        tmp2=N.max(abs(S1-S2))
        print "Max difference kpointhelper: ",max(tmp1,tmp2)

        ee=RA.random(1)+0.0001j
        SGFf90=elec.getSig(ee,k[0:2].copy(),UseF90helpers=True)
        SGF=elec.getSig(ee,k[0:2].copy(),UseF90helpers=False)
        SGFerr=N.max(abs(SGFf90-SGF))
        print "Max difference for self-energy: ",SGFerr
        maxerr=max(tmp1,tmp2,SGFerr,maxerr)

    print "Maximum error : ",maxerr
    print 
    print "###############################################################"
    print "###############################################################"
    if maxerr>1e-9:
        print "ERROR!"
        kuk
    else:
        print "Test passed for remove Xij, distributeg0, readnewtshs, and setkpoint!"
    print "###############################################################"
    print "###############################################################"
    print 

main()
