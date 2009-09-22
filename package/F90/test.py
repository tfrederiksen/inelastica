import sys, profile
sys.path+=['..']

import pyTBT
import numpy as N
import numpy.random as RA

def main():
    elec=pyTBT.surfaceGF('../TestCalculations/Self-energy-FCC111/ELEC-1x1//Au3D_BCA.TSHS',3,3)
    maxerr=0.0

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
        print "Max difference : ",max(tmp1,tmp2)
        maxerr=max(tmp1,tmp2,maxerr)

    print "Maximum error : ",maxerr
    if maxerr>1e-9:
        print "ERROR!"
    else:
        print "Test passed!"
profile.run('main()')
