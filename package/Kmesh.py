# Kmesh module by T. Frederiksen
import MiscMath as MM
import numpy as N

# This function is used in the class below for each of the
# three spatial directions (in momentum space)
def generatelinmesh(Nk):
    "Generate Nk points sampling linearly the interval [-0.5;0.5]"
    kpts = [(ii*1.0+0.5)/Nk-0.5 for ii in range(Nk)]
    wgts = [1.0/Nk for ii in range(Nk)]
    return N.array(kpts), N.array(wgts)

class kmesh:

    def __init__(self,Nk1,Nk2,Nk3,meshtype=['LIN','LIN','LIN'],invsymmetry=True):
        self.Nk = N.array([Nk1,Nk2,Nk3])
        self.type = meshtype
        self.genklists()
        self.genkmesh()
        if invsymmetry:
            self.SymmetryReduce()

    def genklists(self):
        self.kpts = []
        self.wgts = []
        self.errorw = []
        for i in range(3): # loop over the three k-components
            if self.type[i].upper() == 'GK' or self.type[i].upper() == 'GAUSSKRONROD':
                self.type[i] = 'GK'
                if self.Nk[i]>1: # GK-method fails with fewer points
                    kpts, wgts, errorw = MM.GaussKronrod(self.Nk[i])
                    self.kpts.append(kpts)
                    self.wgts.append(wgts)
                    self.errorw.append(errorw)
                else:
                    print 'Kmesh.py: GK method requires Nk=%i>1'%(self.Nk[i])
                    kuk
            elif self.type[i].upper() == 'LIN' or self.type[i].upper() == 'LINEAR':
                self.type[i] = 'LIN'
                kpts, wgts = generatelinmesh(self.Nk[i])
                self.kpts.append(kpts)
                self.wgts.append(wgts)
                self.errorw.append(wgts)
            else:
                print 'Kmesh.py: Unknown meshtype:', self.type[i].upper()
            self.Nk[i] = len(self.kpts[i])
        self.NNk = N.prod(self.Nk)
        self.NGK = len(self.errorw) # Number of GK axes
        print 'Kmesh.py: Generating mesh:'
        print ' ... type = ', self.type
        print ' ... Nk = ', self.Nk
        
    def genkmesh(self):
        # repete out in 3D
        kpts = N.zeros((self.NNk,3)) # Array of k-points
        wgts = N.ones((4,self.NNk)) # (wgts, errorw1, errorw2, errorw3)
        nn = 0
        for i in range(self.Nk[0]):
            for j in range(self.Nk[1]):
                for k in range(self.Nk[2]):
                    kpts[nn,:] = [self.kpts[0][i],self.kpts[1][j],self.kpts[2][k]]
                    wgts[0,nn] = self.wgts[0][i]*self.wgts[1][j]*self.wgts[2][k]
                    wgts[1,nn] = self.errorw[0][i]*self.wgts[1][j]*self.wgts[2][k]
                    wgts[2,nn] = self.wgts[0][i]*self.errorw[1][j]*self.wgts[2][k]
                    wgts[3,nn] = self.wgts[0][i]*self.wgts[1][j]*self.errorw[2][k]
                    nn += 1
        print ' ... NNk = %i, sum(wgts) = %.8f'%(self.NNk,N.sum(wgts[0]))
        print ' ... sum(errorw) = (%.8f,%.8f,%.8f)'%tuple(N.sum(wgts[i+1]) for i in range(3))
        self.kpts = kpts
        self.wgts = wgts

    def SymmetryReduce(self):
        # Remove duplicates for symmetry
        # INVERSION SYMMETRY:
        # If the Bloch function
        #    \psi(k) = exp(ikr)u(k),
        # with crystal momentum k, is an eigenstate of the Schroedinger equation then also
        #    \psi^\dagger(k) = exp(-ikr)u^\dagger(k)
        # with crystal momentum -k, is an eigenstate with same eigenvalue.
        # Hence E(k) = E(-k).
        # TIME REVERSAL SYMMETRY:
        # t,\psi(r,t) --> -t,\psi^\dagger(r,-t). T(k) = T(-k).
        # (Elastic) propagation from L to R is always identical to propagation from R to L.
        print ' ... Applying inversion symmetry (the simple way)'
        # No brute force (and therefore terribly slow) pairing here
        indx = [[ii,2] for ii in range(self.NNk/2,self.NNk)] # Keep the last half of the k-points with double weight
        k0 = self.kpts[self.NNk/2]
        if N.dot(k0,k0)==0: # gamma in the kptsist
            indx[0] = [self.NNk/2,1] # lower weight to one
        indx, weight = N.array([ii[0] for ii in indx]), N.array([ii[1] for ii in indx])
        kpts, wgts = self.kpts[indx], self.wgts[:,indx]*weight
        self.kpts = kpts
        self.NNk = len(kpts)
        self.wgts = wgts
        print ' ... NNk = %i, sum(wgts) = %.8f'%(self.NNk,N.sum(wgts[0]))
        print ' ... sum(errorw) = (%.8f,%.8f,%.8f)'%tuple(N.sum(wgts[i+1]) for i in range(3))

# Test example
if __name__ == '__main__':
    k = kmesh(4,3,1,meshtype=['LIN','GK','LIN'])
    print 'ki wi'
    for i in range(len(k.kpts)):
        print k.kpts[i],k.wgts[0,i]
