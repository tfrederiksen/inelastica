import sys, profile
sys.path+=['../..']

import Inelastica.pyTBT as pyTBT
import numpy as N
import numpy.random as RA

def main():
    elec1=pyTBT.surfaceGF('Self-energy-FCC100/ELEC-1x1//ABAB.TSHS',3,3)
    elec3=pyTBT.surfaceGF('Self-energy-FCC100/ELEC-3x3//ABAB.TSHS',1,1)
    maxerr=0.0

    for ii in range(10):
        ee=RA.random()*10-5+RA.random()*5.0j
        k=RA.random(2)
        print " Checking energy: %f+%fi k-point: %f,%f"%(ee.real,ee.imag,k[0],k[1])
        
        Sig1=elec1.getSig(ee,k,left=True,Bulk=True)
        Sig3=elec3.getSig(ee,k,left=True,Bulk=True)
        tmp=N.max(abs(Sig1-Sig3))
        maxerr=max(tmp,maxerr)
        print 'Left, Bulk : max(abs(Sig1-Sig3)) ',tmp 
        Sig1=elec1.getSig(ee,k,left=False,Bulk=True)
        Sig3=elec3.getSig(ee,k,left=False,Bulk=True)
        tmp=N.max(abs(Sig1-Sig3))
        maxerr=max(tmp,maxerr)
        print 'Right, Bulk : max(abs(Sig1-Sig3)) ',tmp
        Sig1=elec1.getSig(ee,k,left=True,Bulk=False)
        Sig3=elec3.getSig(ee,k,left=True,Bulk=False)
        tmp=N.max(abs(Sig1-Sig3))
        maxerr=max(tmp,maxerr)
        print 'Left, no-Bulk : max(abs(Sig1-Sig3)) ',tmp
        Sig1=elec1.getSig(ee,k,left=False,Bulk=False)
        Sig3=elec3.getSig(ee,k,left=False,Bulk=False)
        tmp=N.max(abs(Sig1-Sig3))
        maxerr=max(tmp,maxerr)
        print 'Right, no-Bulk : max(abs(Sig1-Sig3)) ',tmp

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
