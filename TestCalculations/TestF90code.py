import sys, profile
sys.path+=['..']

import Inelastica.pyTBT as pyTBT
import Inelastica.SiestaIO as SIO
import Inelastica.NEGF as NEGF
import numpy as N
import numpy.random as RA

tol = 1e-9

def compare_H(HS1,HS2,not_checks={}):
    # Check that we read the same TSHS data
    if HS1.gamma!=HS2.gamma or HS1.onlyS!=HS2.onlyS or HS1.nuo!=HS2.nuo \
            or HS1.nspin!=HS2.nspin or HS1.maxnh!=HS2.maxnh or HS1.qtot!=HS2.qtot \
            or HS1.temp!=HS2.temp or HS1.nua!=HS2.nua or HS1.ef!=HS2.ef or HS1.istep!=HS2.istep or HS1.ia1!=HS2.ia1:
        print "Failed to give same results between python and fortran code for simple variables!!!"
        print HS1.gamma,HS2.gamma
        print HS1.onlyS,HS2.onlyS
        print HS1.nuo,HS2.nuo
        print HS1.no,HS2.no
        print HS1.nspin,HS2.nspin
        print HS1.maxnh,HS2.maxnh
        print HS1.qtot,HS2.qtot
        print HS1.temp,HS2.temp
        print HS1.nua,HS2.nua
        print HS1.ef,HS2.ef
        kuk

    if N.sum(N.abs(HS1.cell-HS2.cell)) > tol:
        print "Failed to give same results between python and fortran code for cell!!!"
        kuk
    if N.sum(N.abs(HS1.lasto-HS2.lasto)) +N.sum(N.abs(HS1.numh-HS2.numh)) \
            +N.sum(N.abs(HS1.listh%HS1.nuo-HS2.listh%HS2.nuo)) > tol:
        print "Failed to give same results between python and fortran code for index arrays!!!"
        kuk
    if N.sum(N.abs(HS1.xa-HS2.xa)) > tol:
        print "Failed to give same results between python and fortran code for coordinates!!!"
        kuk
    if N.sum(N.abs(HS1.Ssparse-HS2.Ssparse)) > tol:
        print "Failed to give same results between python and fortran code for overlap!!!"
        kuk
    if N.sum(N.abs(HS1.xij-HS2.xij))>tol and not 'xij' in not_checks:
        print "Failed to give same results between python and fortran code for xij!!!"
        kuk
    if not HS1.onlyS: # Not reading onlyS files
        if N.sum(N.abs(HS1.Hsparse-HS2.Hsparse))>tol:
            print "Failed to give same results between python and fortran code for Hamiltonian!!!"
            kuk

        if False:
            for ii in range(2):
                k = N.array(RA.random(3),N.float)
                print " Checking k-point: %f,%f,%f"%(k[0],k[1],k[2])
                HS1.kpoint = N.zeros((3,),N.float) # To ensure it is calculated!
                HS1.setkpoint(k)
                H1, S1 = HS1.H.copy(), HS1.S.copy()
                HS2.kpoint = N.zeros((3,),N.float) # To ensure it is calculated!
                HS2.setkpoint(k)
                H2, S2 = HS2.H.copy(), HS2.S.copy()
                tmp1 = N.max(abs(H1-H2))
                tmp2 = N.max(abs(S1-S2))
                print "Max difference kpointhelper: ",max(tmp1,tmp2)
                if max(tmp1,tmp2)>1e-9:
                    print "ERROR!"
                    kuk

    print 'Reading of file %s vs. %s'%(HS1.fn,HS2.fn)
    print 'PASSED!\n'

def main():
    if not SIO.F90imported:
        print "To test the F90 routines you better compile them first"
        kuk

    print 'TESTING reading routines:\n'

    # Test readTSHS routines
    for file in ['Self-energy-FCC111/ELEC-1x1/Au3D_BCA.TSHS','sample.onlyS']:
        HS1=SIO.HS('../TestCalculations/'+file,UseF90helpers=False)
        HS2=SIO.HS('../TestCalculations/'+file,UseF90helpers=True)
        compare_H(HS1,HS2)
    
    # Test readTSHS routines new vs. old
    for file in ['Self-energy-FCC100/ELEC-1x1/ABAB.TSHS',
                 'Self-energy-FCC111/ELEC-1x1/Au3D_BCA.TSHS']:      
        HS1=SIO.HS('../TestCalculations/' + file,UseF90helpers=True)
        HS2=SIO.HS('../TestCalculations/' + file.replace('.TSHS','_NEW.TSHS'),UseF90helpers=True)
        compare_H(HS1,HS2,not_checks={'xij':False})

    # Test removeUnitCellXij
    print 'TESTING removeUnitCellXij method'
    elec1 = NEGF.ElectrodeSelfEnergy('../TestCalculations/Self-energy-FCC111/ELEC-1x1/Au3D_BCA.TSHS',1,1,UseF90helpers=True)
    elec2 = NEGF.ElectrodeSelfEnergy('../TestCalculations/Self-energy-FCC111/ELEC-1x1/Au3D_BCA.TSHS',1,1,UseF90helpers=False)
    maxerr=N.max(abs(elec1.HS.xij-elec2.HS.xij))
    print "Maximum difference between Xij :",maxerr 

    for ii in range(10):
        k=N.array(RA.random(3),N.float)
        print " Checking k-point: %f,%f,%f"%(k[0],k[1],k[2])
        elec1.HS.kpoint=N.zeros((3,),N.float) # To ensure it is calculated!
        elec1.HS.setkpoint(k,UseF90helpers=True)
        H1, S1 = elec1.HS.H.copy(), elec1.HS.S.copy()
        elec2.HS.kpoint=N.zeros((3,),N.float) # To ensure it is calculated!
        elec2.HS.setkpoint(k,UseF90helpers=False)
        H2, S2 = elec2.HS.H.copy(), elec2.HS.S.copy()
        tmp1 = N.max(abs(H1-H2))
        tmp2 = N.max(abs(S1-S2))
        print "Max difference kpointhelper: ",max(tmp1,tmp2)
        if maxerr>1e-9:
            print "ERROR!"
            kuk

        ee = RA.random(1)+0.0001j
        elec1.semiinf = 2
        elec2.semiinf = 2
        SGF1 = elec1.getSig(ee,k[0:2].copy(),UseF90helpers=True)
        SGF2 = elec2.getSig(ee,k[0:2].copy(),UseF90helpers=False)
        SGFerr = N.max(abs(SGF1-SGF2))
        print "Max difference for self-energy: ",SGFerr
        maxerr = max(tmp1,tmp2,SGFerr,maxerr)
        if maxerr>1e-9:
            print "ERROR!"
            kuk

    # SurfaceGF
    print '\nTESTING pyTBT.surfaceGF method (1x1 vs 3x3)'
    elecF90 = NEGF.ElectrodeSelfEnergy('../TestCalculations/Self-energy-FCC111/ELEC-1x1/Au3D_BCA.TSHS',3,3,UseF90helpers=True)
    elecNoF90 = NEGF.ElectrodeSelfEnergy('../TestCalculations/Self-energy-FCC111/ELEC-1x1/Au3D_BCA.TSHS',3,3,UseF90helpers=False)

    for ii in range(10):
        k=N.array(RA.random(3),N.float)
        print " Checking k-point: %f,%f,%f"%(k[0],k[1],k[2])
        elecF90.HS.kpoint=N.zeros((3,),N.float) # To ensure it is calculated!
        elecF90.HS.setkpoint(k,UseF90helpers=True)
        H1, S1 = elecF90.HS.H.copy(), elecF90.HS.S.copy()
        elecNoF90.HS.kpoint=N.zeros((3,),N.float) # To ensure it is calculated!
        elecNoF90.HS.setkpoint(k,UseF90helpers=False)
        H2, S2 = elecNoF90.HS.H.copy(), elecNoF90.HS.S.copy()
        tmp1=N.max(abs(H1-H2))
        tmp2=N.max(abs(S1-S2))
        print "Max difference kpointhelper: ",max(tmp1,tmp2)
        if maxerr>1e-9:
            print "ERROR!"
            kuk

        ee=RA.random(1)+0.0001j
        elecF90.semiinf = 2
        elecNoF90.semiinf = 2
        SGFf90=elecF90.getSig(ee,k[0:2].copy(),UseF90helpers=True)
        SGF=elecNoF90.getSig(ee,k[0:2].copy(),UseF90helpers=False)
        SGFerr=N.max(abs(SGFf90-SGF))
        print "Max difference for self-energy: ",SGFerr
        maxerr=max(tmp1,tmp2,SGFerr,maxerr)
        if maxerr>1e-9:
            print "ERROR!"
            kuk

    print "Maximum error : ",maxerr
    print 
    print "###############################################################"
    print "###############################################################"
    if maxerr>1e-9:
        print "ERROR!"
        kuk
    else:
        print "Tests passed for remove Xij, expansion_SE, readTSHS, and setkpoint!"
    print "###############################################################"
    print "###############################################################"
    print 

main()
