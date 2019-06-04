"""

:mod:`Inelastica.physics.mesh`
==============================

Written by Thomas Frederiksen.

Classes
-------

.. autosummary::
   :toctree:

   kmesh

.. currentmodule:: Inelastica.physics.mesh

"""
from __future__ import print_function

import Inelastica.math as MM
import numpy as N

# This function is used in the class below for each of the
# three spatial directions (in momentum space)


def generatelinmesh(Nk):
    "Generate Nk points sampling linearly the interval [-0.5;0.5]"
    kpts = [(ii*1.0+0.5)/Nk-0.5 for ii in range(Nk)]
    wgts = [1.0/Nk for ii in range(Nk)]
    return N.array(kpts), N.array(wgts)


class kmesh(object):

    """
    Create a k-mesh samling where each of the three components
    sample the range [-0.5,0.5]). They are not in reciprocal space.
    Instead they correspond to the mathematical orthogonal space
    that is Fourier transformed.

    Variables:
    k    : Array of 3D k-points (shape = [NNk,3])
    w    : Array of integration weights (shape = [4,NNk]), where
           w[0] are the full weights and
           w[i=1,2,3] are weights with a reduced number of GK points
           for the three axes. These are used for error estimates
           within the GK method.
    NNk  : Total number of k-points explicitly in the mesh
    Nk   : Array with sampling for each axis
    type : List of sampling type for each axis
    """

    def __init__(self, Nk1, Nk2, Nk3, meshtype=['LIN', 'LIN', 'LIN'], invsymmetry=True):
        """
        Returns an instance of a k-mesh with each of the k-vector axes
        sampled either linearly (LIN) or using a Gauss-Kronrod (GK) scheme.

        An axis i sampled by (Nki,LIN) generates Nk points, while
        (Nki,GK) returns 2*Nk+1 points.

        By applying inversion symmetry (k=-k) the number of points in the mesh
        can be reduced by a factor two (except when the Gamma point is included
        which have no partner).
        """
        self.Nk = N.array([Nk1, Nk2, Nk3])
        self.type = meshtype
        self.genkmesh()
        self.invsymmetry = invsymmetry
        if invsymmetry:
            self.SymmetryReduce()

    def genkmesh(self):
        "Generate mesh for a given sampling of the axes"
        self.k = []
        self.w = []
        errorw = []
        for i in range(3): # loop over the three k-components
            if self.type[i].upper() == 'GK' or self.type[i].upper() == 'GAUSSKRONROD':
                self.type[i] = 'GK'
                if self.Nk[i] > 1: # GK-method fails with fewer points
                    kpts, wgts, ew = MM.GaussKronrod(self.Nk[i])
                    self.k.append(kpts)
                    self.w.append(wgts)
                    errorw.append(ew)
                else:
                    print('Kmesh.py: GK method requires Nk=%i>1'%(self.Nk[i]))
                    sys.exit(1)
            elif self.type[i].upper() == 'LIN' or self.type[i].upper() == 'LINEAR':
                self.type[i] = 'LIN'
                kpts, wgts = generatelinmesh(self.Nk[i])
                self.k.append(kpts)
                self.w.append(wgts)
                errorw.append(wgts)
            else:
                print('Kmesh.py: Unknown meshtype:', self.type[i].upper())
            self.Nk[i] = len(self.k[i])
        self.NNk = N.prod(self.Nk)
        print('Kmesh.py: Generating mesh:')
        print(' ... type = ', self.type)
        print(' ... Nk = ', self.Nk)
        # repete out in 3D
        kpts = N.zeros((self.NNk, 3)) # Array of k-points
        wgts = N.ones((4, self.NNk)) # (wgts, errorw1, errorw2, errorw3)
        nn = 0
        for i in range(self.Nk[0]):
            for j in range(self.Nk[1]):
                for k in range(self.Nk[2]):
                    kpts[nn, :] = [self.k[0][i], self.k[1][j], self.k[2][k]]
                    wgts[0, nn] = self.w[0][i]*self.w[1][j]*self.w[2][k]
                    wgts[1, nn] = errorw[0][i]*self.w[1][j]*self.w[2][k]
                    wgts[2, nn] = self.w[0][i]*errorw[1][j]*self.w[2][k]
                    wgts[3, nn] = self.w[0][i]*self.w[1][j]*errorw[2][k]
                    nn += 1
        print(' ... NNk = %i, sum(wgts) = %.8f'%(self.NNk, N.sum(wgts[0])))
        print(' ... sum(errorw) = (%.8f,%.8f,%.8f)'%tuple(N.sum(wgts[i+1]) for i in range(3)))
        self.k = kpts
        self.w = wgts

    def SymmetryReduce(self):
        """
        Remove duplicates for symmetry
        INVERSION SYMMETRY:
        If the Bloch function
           \psi(k) = exp(ikr)u(k),
        with crystal momentum k, is an eigenstate of the Schroedinger equation then also
           \psi^\dagger(k) = exp(-ikr)u^\dagger(k)
        with crystal momentum -k, is an eigenstate with same eigenvalue.
        Hence E(k) = E(-k).
        TIME REVERSAL SYMMETRY:
        t,\psi(r,t) --> -t,\psi^\dagger(r,-t). T(k) = T(-k).
        (Elastic) propagation from L to R is always identical to propagation from R to L.
        """
        print(' ... Applying inversion symmetry (the simple way)')
        # No brute force (and therefore terribly slow) pairing here
        indx = [[ii, 2] for ii in range(self.NNk // 2, self.NNk)]  # Keep the last half of the k-points with double weight
        k0 = self.k[self.NNk // 2]
        if N.dot(k0, k0) == 0: # Gamma in the kptlist
            indx[0] = [self.NNk // 2, 1] # lower weight to one
        indx, weight = N.array([ii[0] for ii in indx]), N.array([ii[1] for ii in indx])
        kpts, wgts = self.k[indx], self.w[:, indx]*weight
        self.k = kpts
        self.NNk = len(kpts)
        self.w = wgts
        print(' ... NNk = %i, sum(wgts) = %.8f'%(self.NNk, N.sum(wgts[0])))
        print(' ... sum(errorw) = (%.8f,%.8f,%.8f)'%tuple(N.sum(wgts[i+1]) for i in range(3)))

    def mesh2file(self, fn):
        "Writes the k-mesh to file"
        f = open(fn, 'w')
        f.write('# k1 k2 k3 w0 w1 w2 w3\n')
        for i in range(self.NNk):
            s = ''
            for j in range(3):
                s += '%.8f '%self.k[i, j]
            for j in range(4):
                s += '%.8e '%self.w[j, i]
            f.write(s+'\n')
        f.close()


def test():
    """
    Test function
    """
    mesh = kmesh(4, 3, 3, meshtype=['LIN', 'GK', 'LIN'], invsymmetry=False)
    keys = list(mesh.__dict__.keys())
    print(keys)
    print(mesh.Nk)
    print(mesh.NNk)
    print(mesh.type)
    print(mesh.invsymmetry)
    #print 'ki wi'
    #for i in range(len(mesh.k)):
    #    print mesh.k[i], mesh.w[0, i]
    mesh.mesh2file('mesh-test.dat')

    print('Integrate some simple functions over [-0.5,0.5]:')
    print('   f(x,y,z)=1 => \int f dxdydz =', N.sum(mesh.w[0]))
    for i, s in enumerate(['x', 'y', 'z']):
        f = mesh.k[:, i]
        print('   f(x,y,z)=%s => \int f dxdydz ='%s, N.sum(f*mesh.w[0]))
    f = mesh.k[:, 0]*mesh.k[:, 1]*mesh.k[:, 2]
    print('   f(x,y,z)=x*y*z => \int f dxdydz =', N.sum(f*mesh.w[0]))
    f = (mesh.k[:, 0]+1)*(mesh.k[:, 1]+1)*(mesh.k[:, 2]+1)
    print('   f(x,y,z)=(x+1)*(y+1)*(z+1) => \int f dxdydz =', N.sum(f*mesh.w[0]))
    for i, s in enumerate(['x', 'y', 'z']):
        f = N.cos(2*mesh.k[:, i])
        print('   f(x,y,z)=cos(2%s) => \int f dxdydz ='%s, N.sum(f*mesh.w[0]), '[exact: sin(1) ~ 0.841470984807897]')

if __name__ == '__main__':
    test()
