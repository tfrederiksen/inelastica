import sys, profile
sys.path+=['../..']

import Inelastica.NEGF as NEGF
import numpy as N
import numpy.random as RA

def err_SE_k(ee,E1,E2,k,**kwargs):
    E1.semiinf = 2
    E2.semiinf = 2
    Sig1 = E1.getSig(ee,k,**kwargs)
    Sig2 = E2.getSig(ee,k,**kwargs)
    return N.max(abs(Sig1-Sig2))

def comp_SE(E1,E2):
    maxerr = 0.
    for ii in range(5):
        ee = RA.random()*10-5+RA.random()*5.0j
        k = RA.random(2)
        print " Checking energy: %f+%fi k-point: %f, %f"%(ee.real,ee.imag,k[0],k[1])
        tmp = err_SE_k(ee,E1,E2,k,left=True,Bulk=True)
        maxerr = max(tmp,maxerr)
        print 'Left, Bulk : max(abs(Sig1-Sig2)) ',tmp 
        tmp = err_SE_k(ee,E1,E2,k,left=False,Bulk=True)
        maxerr = max(tmp,maxerr)
        print 'Right, Bulk : max(abs(Sig1-Sig2)) ',tmp
        tmp = err_SE_k(ee,E1,E2,k,left=True,Bulk=False)
        maxerr = max(tmp,maxerr)
        print 'Left, no-Bulk : max(abs(Sig1-Sig2)) ',tmp
        tmp = err_SE_k(ee,E1,E2,k,left=False,Bulk=False)
        maxerr = max(tmp,maxerr)
        print 'Right, no-Bulk : max(abs(Sig1-Sig2)) ',tmp
    return maxerr
    

def main():
    maxerr = 0.0

    elec1 = NEGF.ElectrodeSelfEnergy('Self-energy-FCC100/ELEC-1x1/ABAB.TSHS',3,3)
    elec3 = NEGF.ElectrodeSelfEnergy('Self-energy-FCC100/ELEC-3x3/ABAB.TSHS',1,1)

    maxerr = comp_SE(elec1,elec3)
    del elec1,elec3
    
    print
    print "#########################################"
    print "#########################################"
    print "Maximum error : ",maxerr
    if maxerr>0.001:
        print "ERROR!"
        kuk
    else:
        print "Test passed!"
    print "#########################################"
    print "#########################################"
    print

main()
