"""

Symmetry (:mod:`Inelastica.Symmetry`)
=====================================

Symmetry operations of lattice and basis.

.. currentmodule:: Inelastica.Symmetry

classes
-------

.. autosummary::
   :toctree:

   Symmetry

"""

import Inelastica.io.siesta as SIO
import Inelastica.MiscMath as MM
import numpy as N
import sys
mm = MM.mm


class Symmetry:
    """
    Classify symmetry of lattice and basis.
    Class contain:
    # Basic quantities
    NN  : Number of atoms
    xyz : xyz coordinates [NN,3]
    snr : Siesta numbering of atoms [NN]
    anr : atomic number [NN]
    pbc : periodic boundary conditions [3,3] (first index is basis vector,
                                              second is x/y/z)
    # Derived quantities
    a1..3    : Lattice vectors of unitcell
    b1..3    : Reciprocal lattice vectors
    bp1..bp3 : Reciprocal lattice vectors to PBC.
    latticeType : string CUBIC/FCC/BCC/HEX/TETRAGONAL/ORTHORHOMBIC
    basis.xyz : xyz of basis inside a1..a3
    basis.snr, .anr, .NN : Siesta/atomic number and number of atoms

    pointU33 : List of lattice point symmetries [?,3,3], i.e., the
             rotation/mirror symmetry operators in space

    # Full set of point operations around rotation center
    U33      : List of symmetry operations for lattice + basis, i.e.,
               subset of pointU33 [?,3,3]
    origo    : List of rotation origin for the operations in fullU33

    # Irreducible symmetry operations not needed anymore

    All arrays are numpy type.

    Use:
       sym = Symmetry.Symmetry(fn='filename.XV',accuracy=1e-4)
       Where accuracy is the fudge factor in Ang.
    """

    def __init__(self, fn=None, accuracy=1e-4, onlyLatticeSym=False):
        self.accuracy=accuracy
        if fn!=None:
            self.readXV(fn, onlyLatticeSym)

    ###########################################################
    # Init

    def readXV(self, fn, onlyLatticeSym=False):
        # Get geometry
        self.pbc, self.snr, self.anr, self.xyz = SIO.ReadXVFile(fn)
        self.pbc, self.snr, self.anr, self.xyz = N.array(self.pbc), N.array(self.snr), N.array(self.anr), N.array(self.xyz)
        self.NN = len(self.anr)
        self.findSymmetry(onlyLatticeSym)

    def setupGeom(self, pbc, snr, anr, xyz, onlyLatticeSym=False):
        self.pbc, self.snr, self.anr, self.xyz = N.array(pbc), N.array(snr), N.array(anr), N.array(xyz)
        self.NN = len(self.anr)
        self.findSymmetry(onlyLatticeSym)

    ###########################################################
    # Main symmetry finding

    def findSymmetry(self, onlyLatticeSym=False):
        # Reciprocal lattice for PBC
        pbc = self.pbc
        bp1, bp2, bp3 = N.cross(pbc[1], pbc[2]), N.cross(pbc[0], pbc[2]),\
            N.cross(pbc[0], pbc[1])
        bp1, bp2, bp3 = bp1/mm(pbc[0], bp1), bp2/mm(pbc[1], bp2),\
            bp3/mm(pbc[2], bp3)
        self.bp1, self.bp2, self.bp3 = bp1, bp2, bp3

        # Find minimal unit cell
        self.findLattice()

        # Calculate point group symmetry operators for lattice system
        self.latticeGroup()

        if not onlyLatticeSym:
            # Point group including basis (lower than for lattice system)
            self.pointGroup()

            ## Reduce to smallest possible subgroup which generate the full group
            #self.findIrreducible()

        return

    ###########################################################
    # Make the dynamical matrix symmetric

    def symmetrizeFC(self, FC, FCfirst, FClast, radi=0.0):
        NN, NU, NFC = self.NN, len(self.U33), FClast-FCfirst+1
        if self.basis.NN>NFC:
            print "Phonons: ERROR: FCfirst/last do contain all atoms in the basis (%i)."%self.basis.NN
            sys.exit('Symmetry error:')

        # Rearrange to FC_ia,jb, force from atom j, axis b to atom i, axis a
        if len(FC.shape)==3:
            FCreshape = True
            FCn = N.zeros((NFC, 3, self.NN, 3))
            for ii in range(NFC):
                for jj in range(3):
                    FCn[ii, jj, :, :] = FC[ii*3+jj, :, :]
            FC = FCn
        else:
            FCreshape = False

        # Rearange basis to fit FCfirst...FClast order in Siesta FC file
        basisxyz = moveIntoCell(self.xyz[FCfirst-1:FClast],\
                                self.a1, self.a2, self.a3, self.accuracy)
        ipiv = []
        for elem in basisxyz[0:FClast-FCfirst+1]:
            for jj, comp in enumerate(self.basis.xyz):
                if N.allclose(elem, comp, atol=self.accuracy):
                    ipiv+=[jj]
                    break
        if len(N.unique(ipiv))!=self.basis.NN:
            print "Symmetry: Some confusion about the order of the basis of atoms and the FCfirst/last region."
            sys.exit('Symmetry: This should never happen')
        self.basis.xyz = self.basis.xyz[ipiv, :]
        self.basis.snr, self.basis.anr = N.array(self.basis.snr)[ipiv], N.array(self.basis.anr)[ipiv]

        # Find out which basis atom corresponds to each atom
        xyz = moveIntoCell(self.xyz, self.a1, self.a2, self.a3, self.accuracy)
        self.basisatom = N.zeros((self.NN), N.int)
        for ii in range(self.basis.NN):
            indx = N.where(N.sum(N.abs(xyz-self.basis.xyz[ii, :]), axis=1)<self.accuracy)
            self.basisatom[indx[0]]=ii

        # Symmetry operations are complicated by the periodic structure
        # in the Siesta calculation. So ... limit range of interaction.

        # Find radie and center for biggest sphere fitting into pbc
        if radi==0.0:
            radi = findRadi(self.pbc[0], self.pbc[1], self.pbc[2])*\
                (1-5*self.accuracy)
        print "Symmetry: Longest range of forces %f Ang."%radi

        # Symmetry operation on FC is composed of:
        #     FC_ia,jb = U^t_aa' PL^j_ii' FC_i'a',j'b' PR_j'j UR_b'b
        # (where 'ed are summed). I.e., normal matrix mult for the L/R
        # unitary matrices. The permutation matrices are complicated.
        # PR : U may change which atom in the basis is moved.
        # PL : Should point to the atom closest to the atom moved and
        #      thus may depend also on j, i.e., the atom moved

        UL = [N.transpose(self.U33[ii]) for ii in range(NU)]
        UR = [self.U33[ii] for ii in range(NU)]
        PL  = N.zeros((NU, NFC, NN, NN), N.int)
        PR  = N.zeros((NU, NFC), N.int)
        for iU in range(NU):
            for ii in range(NFC):
                SIO.printDone(iU*NFC+ii, NU*NFC, 'Symmetrizing')

                # Figure out which atoms are connected by symmetry op.
                xyz = self.xyz.copy()-self.origo[iU]
                FCxyz = xyz[FCfirst-1+ii, :].copy()
                xyz = moveIntoClosest(N.array(xyz)-FCxyz,\
                                    self.pbc[0], self.pbc[1], self.pbc[2])+FCxyz
                nxyz = N.transpose(mm(self.U33[iU], N.transpose(xyz)))
                # Did the moving atom move between unit cells?
                nFCxyz, iix = nxyz[FCfirst-1+ii, :].copy(),  10
                for ix in range(-3, 4):
                    for iy in range(-3, 4):
                        for iz in range(-3, 4):
                            if N.any(N.sum(N.abs(nFCxyz+ix*self.a1+iy*self.a2+iz*self.a3-xyz[FCfirst-1:FClast]), axis=1)<self.accuracy):
                                iix, iiy, iiz = ix, iy, iz
                if iix==10: sys.exit('Symmetry Error')
                tFCxyz = nFCxyz+iix*self.a1+iiy*self.a2+iiz*self.a3
                shift = tFCxyz-nFCxyz # Shift all new coordinates
                nxyz = moveIntoClosest(N.array(nxyz)+shift-FCxyz,\
                                    self.pbc[0], self.pbc[1], self.pbc[2])+FCxyz

                # Which atoms are within radi?
                indx = N.where(distance(xyz-FCxyz)<radi)[0]

                # Find the target atom
                diff = N.sum(N.abs(nxyz[indx].reshape((-1, 1, 3))-\
                         xyz.reshape((1, -1, 3))), axis=2)<self.accuracy
                tindx = N.where(diff)
                tindx2 = tindx[1]
                indx3= indx[tindx[0]]
                if len(indx3)!=len(indx):
                    for ix in range(-1, 2):
                        for iy in range(-1, 2):
                            for iz in range(-1, 2):
                                if ix!=0 or iy!=0 or iz!=0:
                                    tindxs = N.where(N.sum(N.abs(nxyz[indx].reshape((-1, 1, 3))+\
                                                                  self.pbc[0]*ix+self.pbc[1]*iy+self.pbc[2]*iz-\
                                                                  xyz.reshape((1, -1, 3))), axis=2)<self.accuracy)
                                    indx3 = N.concatenate((indx3, indx[tindxs[0]]))
                                    tindx2 = N.concatenate((tindx2, tindxs[1]))
                if len(indx3)!=len(indx):
                    sys.exit('Symmetry error')
                # Make permutation matrix
                PL[iU, ii, tindx2, indx3] = 1
                indx2 = N.where(distance(xyz-tFCxyz)<self.accuracy)
                PR[iU, ii] = self.basisatom[indx2[0]]

                # TF: Check that the symmetry operations apply to FC??
                for kk in range(len(indx3)):
                    oFC = FC[ii, :, indx3[kk], :].reshape((3, 3))
                    sFC = mm(UL[iU], FC[indx2[0], :, tindx2[kk], :].reshape((3, 3)), UR[iU])
                    if N.max(N.abs(oFC-sFC))>0.5:
                        print oFC
                        print sFC

        def applySym(FC):
            FCn = N.tensordot(UR[iU], FC, ((1, 1)))
            FCn = N.swapaxes(FCn, 0, 1)
            FCn = N.tensordot(FCn, UL[iU], ((3, 0)))
            FCn2=0*FCn
            for ii in range(NFC):
                tmp = N.tensordot(PL[iU, ii, :, :], \
                                      FCn[ii, :, :, :], ((1, 1)))
                tmp = N.swapaxes(tmp, 0, 1)
                FCn2[PR[iU, ii], :, :, :] = tmp
            return FCn2

        # Symmetrize dynamical matrix
        print "Symmetry: Iterative application of symmetries"
        niter, change = 0, 10
        FCo = FC.copy()
        # Uncomment the two lines below to skip the application of symmetries
        #niter, change = 10, 10
        #FCs = FC.copy()
        while niter<10 and change>1e-10:
            FCs = FCo
            for iU in range(NU):
                FCn = FCs
                FCs = 0
                for jj in range(self.rankU33[iU]):
                    FCn = applySym(FCn)
                    FCs += FCn
                FCs = FCs/self.rankU33[iU]

            change = N.max(N.abs(FCs-FCo))
            print "Symmetry: max change between old and symmetrized = %e"%change
            print "Symmetry: relative change = %e\n"%(N.max(N.abs(FCs-FCo))/N.max(abs(FCo)))
            FCo = FCs
        if N.max(N.abs(FCs-FC))/N.max(abs(FC))>0.05:
            print "Symmetry: WARNING: large relative difference"

        # Change format back ...
        if FCreshape:
            FC = N.zeros((NFC*3, self.NN, 3))
            for ii in range(NFC):
                for jj in range(3):
                    FC[ii*3+jj, :, :] = FCs[ii, jj, :, :]
            FCs = FC

        return FCs

    ###########################################################
    # NOT USED currently!
    # Find irreducible symmetry operations and print
    #
    def findIrreducible(self):
        # Just print for the lattice group
        ipiv, sign, rotn = self.reduce(self.pointU33)

        # Go through the different rotation centers
        lasto, nU, nO, nR = 0, [], [], []
        for ii in range(1, len(self.origo)):
            if not N.allclose(self.origo[ii-1], self.origo[ii], atol=self.accuracy):
                ipiv, sign, rotn = self.reduce(N.concatenate((self.U33[lasto:ii], [N.eye(3)])))
                if len(rotn)>1 or rotn[0]!=1:
                    for jj, kk in enumerate(lasto+N.array(ipiv)):
                        nU += [self.U33[kk]]
                        nO += [self.origo[kk]]
                        nR += [rotn[jj]]
                    print "Symmetry: Irreducible point operations around %f %f %f"%(self.origo[ii-1][0], self.origo[ii-1][1], self.origo[ii-1][2])
                    print "Ranks*determinant : ", rotn*sign
                lasto = ii

        if len(self.U33)!=0:
            ipiv, sign, rotn = self.reduce(N.concatenate((self.U33[lasto:], [N.eye(3)])))
            for jj, kk in enumerate(lasto+N.array(ipiv)):
                nU += [self.U33[kk]]
                nO += [self.origo[kk]]
                nR += [rotn[jj]]

            print "Symmetry: Irreducible point operations around %f %f %f"%(self.origo[ii-1][0], self.origo[ii-1][1], self.origo[ii-1][2])
            print "Ranks*determinant : ", rotn*sign

        self.fullU33, self.fullOrigo = self.U33, self.origo
        self.U33, self.origo, self.rankU33 = nU, nO, nR
        return

    def reduce(self, Ulist):
        # Find least number of symmetry operations that generate the whole group
        keep = [ii for ii in N.array(Ulist).copy()]
        ChT = N.zeros((len(keep), len(keep)), N.int) # Character table
        for ii in range(len(keep)):
            for jj in range(len(keep)):
                tmp = mm(keep[ii], keep[jj])
                indx=N.where(N.sum(N.sum(N.abs(keep-tmp), axis=2), axis=1)<self.accuracy)
                ChT[ii, jj]=indx[0][0]

        # Calculate rank
        rotn = []
        for ii in keep:
            a, M = ii, 1
            while not N.allclose(N.eye(3), a, atol=self.accuracy):
                a, M =mm(a, ii), M+1
            rotn+=[M]

        def pick(x, y):
            tmp=[]
            for ii in y:
                tmp=N.concatenate((tmp, x[ii, N.array(y+self.accuracy, N.int)]))
            return N.array(tmp).reshape((-1, ))

        def grow(x):
            # See which points in table that can be reached from the points in x
            NN = -1
            while len(x)!=NN:
                NN=len(x)
                x= N.unique(N.concatenate((x, pick(ChT, x))))
            return x

        # Sort and start with largest rank
        order = N.argsort(-abs(N.array(rotn))).reshape((-1, ))
        irr = grow(N.array([order[0]], N.int))
        basis = [order[0]]
        for ii in order:
            if N.sum(irr==ii)==0:
                basis += [ii]
                irr = grow(N.unique(N.concatenate((basis, irr))))

        sign = N.array([N.linalg.det(keep[ii]) for ii in basis])
        keep = N.array([keep[ii] for ii in basis])

        rotn = []
        for ii in keep:
            a, M = ii, 1
            while not N.allclose(N.eye(3), a, atol=self.accuracy):
                a, M =mm(a, ii), M+1
            rotn+=[M]

        return basis, sign, rotn

    ###########################################################
    # point group of lattice+basis

    def pointGroup(self):
        from itertools import product as iprod
        import numpy.linalg as LA
        # Find the symmetry operators and centers where the basis is unchanged

        # Use the least common siesta atom type to find possible centers
        abundance = [N.sum(N.array(self.basis.snr)==snr) for snr in self.basis.snr]
        whichsnr  = self.basis.snr[N.where(N.array(abundance)==N.min(abundance))[0][0]]

        Ulist, a1, a2, a3 = self.pointU33, self.a1, self.a2, self.a3
        pointU,  pointO = [],  []
        for iU, U in enumerate(Ulist):
            SIO.printDone(iU, len(Ulist), 'Looking for point group')

            xyz= self.basis.xyz[N.where(self.basis.snr==whichsnr)[0]]
            # Find centers by solving U(x_i-o)==x_j-o +R where o is center, R is lattice vector
            centers = []
            A = U-N.eye(3)
            for i1, i2, i3 in iprod(range(-1, 2), repeat=3):
                for ii, jj in iprod(range(len(xyz)), repeat=2):
                    sol = LA.lstsq(A, mm(U, xyz[ii, :].transpose())-xyz[jj, :].transpose()+(a1*i1+a2*i2+a3*i3))[0]
                    correct = not N.any(N.abs(mm(A, sol)-mm(U, xyz[ii, :].transpose())+xyz[jj, :].transpose()-(a1*i1+a2*i2+a3*i3))>1e-6)
                    if correct: centers += [sol]
            centers = myUnique2(moveIntoCell(N.array(centers), a1, a2, a3, self.accuracy), self.accuracy)
            # All centers are not correct since we only looked at one atom type and
            # remove repeats of the same symmetry by looking at ipvi, i.e., which atom moves onto which
            origin, ipiv = [], []
            for o in centers:
                xyz, nxyz = self.basis.xyz, mm(U, self.basis.xyz.transpose()-o.reshape((3, 1))).transpose()+o
                xyz, nxyz = moveIntoCell(xyz, a1, a2, a3, self.accuracy), moveIntoCell(nxyz, a1, a2, a3, self.accuracy)
                thisipiv = []
                for x in xyz:
                    for iy, y in enumerate(nxyz):
                        if N.allclose(x, y, atol=self.accuracy):
                            thisipiv+=[iy]

                trueSym = len(thisipiv)==self.basis.NN and N.allclose(N.array(self.basis.snr)[thisipiv], self.basis.snr) and N.allclose(N.sort(thisipiv), N.arange(self.basis.NN))
                if trueSym and len(myIntersect(N.array(ipiv), N.array([thisipiv]), self.accuracy))==0:
                    origin, ipiv = origin+[o], ipiv+[thisipiv]
            if len(origin)>0:
                pointU+=[U]
                pointO+=[origin]

        self.U33, self.origo = pointU, pointO
        print("Symmetry: %i point symmetry operations (with rotation centers) found for lattice+basis."%len(pointU))

        # Calculate rank of operations
        self.rankU33, sign = [], []
        for U in self.U33:
            tmp, ii = U, 1
            while ii<7 and not N.allclose(N.eye(3), tmp, atol=self.accuracy):
                tmp, ii = mm(tmp, U), ii+1
            if ii > 6: sys.exit('Symmetry error: rank >6 !!')
            self.rankU33 += [ii]
            sign += [N.linalg.det(U)]
        print "Symmetry: rank*det"
        print N.array(self.rankU33)*sign
        return

    ###########################################################
    # point group of lattice

    def latticeGroup(self):
        from itertools import product as iprod
        # Generate all unitary matrices that take a1...a3 to
        #   lattice points

        a1, a2, a3= self.a1, self.a2, self.a3
        b1, b2, b3= self.b1, self.b2, self.b3
        points, dist = [], []
        # Generate lattice points
        for i1, i2, i3 in iprod(range(-2, 3), repeat=3):
            points +=[i1*a1+i2*a2+i3*a3]
            dist+=[distance(points[-1])]
        points,  dist = N.array(points),  N.array(dist)
        a = [a1, a2, a3]
        targets = [[], [], []]
        # Choose all possible target points for a1...a3 with same length
        for ii in range(3):
            targets[ii]=points[N.where(abs(dist-distance(a[ii]))<self.accuracy)]
        # Make all possible combinations of targets for a1..a3
        possible = myPermute(targets)
        # Generate unitary transformation from a1, a2, a3 -> possible
        Ulist = []
        for ii in possible:
            # Generate symmetry operations
            U = ii[0].reshape(3, 1)*b1.reshape(1, 3)+\
                ii[1].reshape(3, 1)*b2.reshape(1, 3)+\
                ii[2].reshape(3, 1)*b3.reshape(1, 3)
            if N.allclose(mm(U, N.transpose(U)), N.eye(3), atol=self.accuracy):
                Ulist+=[U]

        self.pointU33 = Ulist
        print "Symmetry: %i point symmetry operations found for lattice"%len(Ulist)
        sign = N.array([N.linalg.det(ii) for ii in Ulist])

        rotn = []
        for ii in Ulist:
            a, M = ii, 1
            while not N.allclose(N.eye(3), a, atol=self.accuracy):
                a, M =mm(a, ii), M+1
            rotn+=[M]

        self.pointRank = rotn
        print "Symmetry: rank * determinant"
        print sign*self.pointRank

        return

    ###########################################################
    # Minimum lattice and supercell

    def findLattice(self):
        """
        Find lattice vectors a1..3, reciprocal vectors b1..3
        such that a_i b_j=delta_ij
        """
        # Use least abundant atom type
        # Use the least common siesta atom type to find possible centers
        abundance = [N.sum(N.array(self.snr)==snr) for snr in self.snr]
        whichsnr  = self.snr[N.where(N.array(abundance)==N.min(abundance))[0][0]]

        # Brute force ... find vectors connecting same type of atoms,
        # try and see if they are lattice vectors starting with shortest
        sameatoms = N.argwhere(self.snr==whichsnr)
        samexyz = self.xyz[sameatoms]
        samexyz = samexyz.reshape((-1, 3))
        samexyz = moveIntoCell(samexyz, self.pbc[0, :], self.pbc[1, :],\
                                   self.pbc[2, :], self.accuracy)
        # Make matrix of difference between coordinates
        possible = samexyz.reshape((1, -1, 3))-samexyz.reshape((-1, 1, 3))
        # Reshape into list of vectors
        possible = possible.reshape((-1, 3))
        # Add the periodic boundary conditions in case the unit
        # cell is the whole cell.
        possible = N.concatenate((possible, self.pbc))
        # remove duplicates etc
        possible = myUnique(possible, self.accuracy)

        # Go through combinations and check if possible
        i1, i2, i3, done, NP = 0, 1, 2, False, len(possible)

        def increase(i1, i2, i3, NP, incindx):
            # Keep track of indices
            if incindx<2:
                i3=0
            if incindx<1:
                i2=0
            if incindx==2:
                i3=(i3+1)%NP
                if i3==0:
                    i3=i2+2
                    incindx=1
            if incindx==1:
                i2=(i2+1)%NP
                if i2==0:
                    i2=i1+2
                    incindx=0
            if incindx==0: i1+=1
            i2 = max(i1+1, i2)
            i3 = max(i2+1, i3)
            #if incindx<=1: print i1, i2, i3
            return i1, i2, i3

        # Step through possibilities
        while i1<NP-2 and not done:
            a1, a2 = possible[i1], possible[i2]
            if distance(N.cross(a1, a2))>self.accuracy: # Non-parallell
                a3 = possible[i3]
                incindx = self.checkLattice(a1, a2, a3, samexyz)
                if incindx!=-1:
                    i1, i2, i3 = increase(i1, i2, i3, NP, incindx)
                else:
                    done = True
            else:
                i1, i2, i3 = increase(i1, i2, i3, NP, 1)
            print 'AQUI', NP, i1, i2, i3, incindx

        if not done:
            # Should not be necessary since PBC is added to possible!
            print("Symmerty: WARNING Could not find smaller unitcell")
            sys.exit('Symmetry error')

        # Rearange to fit with PBC
        a1, a2, a3 = self.makeHumanReadable(a1, a2, a3)

        b1, b2, b3 = N.cross(a2, a3), N.cross(a3, a1), N.cross(a1, a2)
        b1, b2, b3 = b1/mm(a1, b1), b2/mm(a2, b2), b3/mm(a3, b3)

        self.a1, self.a2, self.a3 = a1, a2, a3
        self.b1, self.b2, self.b3 = b1, b2, b3

        indx = N.round(abs(mm(self.pbc, N.transpose(N.array([b1, b2, b3])))))
        N1,  N2,  N3 =indx[0, 0],  indx[1, 1],  indx[2, 2]
        print "Symmetry: Lattice structure"
        print "%i atoms in the basis"%(self.NNbasis)
        print "a1 = (%f,%f,%f), N1=%i"%(a1[0], a1[1], a1[2], N1)
        print "a2 = (%f,%f,%f), N2=%i"%(a2[0], a2[1], a2[2], N2)
        print "a3 = (%f,%f,%f), N3=%i"%(a3[0], a3[1], a3[2], N3)
        print
        print "b1 = (%f,%f,%f)"%(b1[0], b1[1], b1[2])
        print "b2 = (%f,%f,%f)"%(b2[0], b2[1], b2[2])
        print "b3 = (%f,%f,%f)"%(b3[0], b3[1], b3[2])

        # Shift to more convenient origin
        xyz = moveIntoCell(self.xyz, a1, a2, a3, self.accuracy)
        tmp=mm(N.array([b1, b2, b3]), N.transpose(xyz))
        mx, my, mz = min(tmp[0, :]), min(tmp[1, :]), min(tmp[2, :])
        shift = mx*a1+my*a2+mz*a3
        print "Shifting coordinates by : ", -shift
        self.xyz = self.xyz-shift

        # Find basis
        xyz = moveIntoCell(self.xyz, a1, a2, a3, self.accuracy)
        class basis:
            pass
        basis.xyz, basis.snr, basis.anr = [], [], []
        for ii in range(len(xyz)):
            for jj in range(len(basis.xyz)):
                if N.allclose(xyz[ii], basis.xyz[jj], atol=self.accuracy):
                    break
            else:
                basis.xyz+=[xyz[ii]]
                basis.snr+=[self.snr[ii]]
                basis.anr+=[self.anr[ii]]
        basis.xyz = N.array(basis.xyz)
        basis.NN = len(basis.xyz)
        if len(basis.anr)!=self.NNbasis:
            print "Symmetry: ERROR! Inconsistent number of basis atoms. Probably bug in program."
            sys.exit('Symmetry error')
        self.basis = basis

        # Find out which basis atom corresponds to each atom
        xyz = moveIntoCell(self.xyz, self.a1, self.a2, self.a3, self.accuracy)
        self.basisatom = N.zeros((self.NN), N.int)
        for ii in range(self.basis.NN):
            indx = N.where(N.sum(N.abs(xyz-self.basis.xyz[ii, :]), axis=1)<self.accuracy)
            self.basisatom[indx[0]]=ii

        # Determine types of lattices from the angles / lengths
        # Sort on length
        length = N.array([distance(ii) for ii in [a1, a2, a3]])
        ipiv = N.argsort(length)
        length = N.array(length)[ipiv]
        a123 = N.array([a1, a2, a3])[ipiv]

        # Calculate angles in degrees
        angles = []
        for ii, a in enumerate([a123[0], a123[1]]):
            for jj, b in enumerate([a123[1], a123[2]][ii:]):
                angles += [N.round(N.arccos(mm(a, b)/distance(a)/distance(b))*\
                                   180/N.pi)]
        # Move angels into range [0,90]
        angles = [min(ii, 180-ii) for ii in angles]
        print 'Lattice angles  =', angles
        print 'Lattice lengths =', length
        latticetype = 'UNKNOWN'
        if N.allclose(length[0], length[1], atol=self.accuracy) and \
           N.allclose(length[1], length[2], atol=self.accuracy):
            # All lattice vectors with same length |a1|=|a2|=|a3|:
            if N.sum(N.abs(N.array(angles)-90))<2:
                latticetype = 'CUBIC'
                # a1 = [a,0,0]
                # a2 = [0,a,0]
                # a3 = [0,0,a]
            elif N.sum(N.abs(N.array(angles)-60))<2:
                latticetype = 'FCC'
                # a1 = [0,a/2,a/2]
                # a2 = [a/2,0,a/2]
                # a3 = [a/2,a/2,0]
            elif N.sum(N.abs(N.array(angles)-70.5))<4:
                latticetype = 'BCC'
                # a1 = [-a/2,a/2,a/2]
                # a2 = [a/2,-a/2,a/2]
                # a3 = [a/2,a/2,-a/2]
        elif N.allclose(length[0], length[1], atol=self.accuracy):
            # |a1|=|a2|:
            if N.abs(angles[0]-60)<2 and N.sum(N.abs(N.array(angles[1:])-90))<2:
                latticetype="HEX"
                # a1 = [a/2,-a*3**.5/2,0]
                # a2 = [a/2,+a*3**.5/2,0]
                # a3 = [0,0,c]
            elif N.sum(N.abs(N.array(angles)-90))<2:
                latticetype = 'TETRAGONAL'
                # a1 = [a,0,0]
                # a2 = [0,a,0]
                # a3 = [0,0,c]
        elif N.sum(N.abs(N.array(angles)-90))<2:
            latticetype="ORTHORHOMBIC"
                # a1 = [a,0,0]
                # a2 = [0,b,0]
                # a3 = [0,0,c]
        print "Symmetry: Lattice = %s"%latticetype
        self.latticeType = latticetype

    def makeHumanReadable(self, a1, a2, a3):
        from itertools import product as iprod
        # TODO! Does not give "Human" readbility ...
        # Choose vectors that give the smallest angles between a1..a3

        # Sort on length
        ipiv = N.argsort(distance(N.array([a1, a2, a3])))
        na1 = [a1, a2, a3][ipiv[0]]
        na2 = [a1, a2, a3][ipiv[1]]
        na3 = [a1, a2, a3][ipiv[2]]
        a1,  a2,  a3 = na1,  na2,  na3

        # Make possible n a1 + m a2 + l a3
        poss = []
        for i1, i2, i3 in iprod(range(-1, 2), repeat=3):
            poss += [a1*i1+a2*i2+a3*i3]
        poss = N.array(poss)

        # possiblities for a1, a2, a3 based on length
        poss1 = poss[N.where(N.abs(distance(poss)-distance(a1))<self.accuracy)]
        poss2 = poss[N.where(N.abs(distance(poss)-distance(a2))<self.accuracy)]
        poss3 = poss[N.where(N.abs(distance(poss)-distance(a3))<self.accuracy)]
        poss = N.concatenate((poss1, poss2, poss3)) # Keep duplicates

        # Choose the one fitting the PBC best
        vec = self.pbc[0, :]
        scal = mm(vec, N.transpose(poss))
        a1=poss[N.where(N.sort(scal)[-1]==scal)[0][0]]

        # Remove linear dependent
        poss=poss[N.where(N.abs(N.abs(mm(poss, a1))-\
                                    N.abs(distance(poss)*distance(a1)))>self.accuracy)]

        # Choose the one fitting the PBC best
        vec = self.pbc[1, :]
        scal = mm(vec, N.transpose(poss))
        a2=poss[N.where(N.sort(scal)[-1]==scal)[0][0]]

        # Remove linear dependent
        poss=poss[N.where(N.abs(N.abs(mm(poss, a2))-\
                                    N.abs(distance(poss)*distance(a2)))>self.accuracy)]

        # Choose the one fitting the PBC best
        vec = self.pbc[2, :]
        scal = mm(vec, N.transpose(poss))
        a3=poss[N.where(N.sort(scal)[-1]==scal)[0][0]]

        if N.linalg.det([a1, a2, a3])<-self.accuracy:
            a3=-a3

        return a1, a2, a3

    def checkLattice(self, a1, a2, a3, samexyz):
        # Check that a1..3 span lattice, return -1 for success
        #  otherwise index to increase

        # Test of PBC
        b1, b2, b3 = N.cross(a2, a3), N.cross(a1, a3), N.cross(a1, a2)
        if distance(b1)<self.accuracy or distance(b2)<self.accuracy or\
                abs(mm(a1, N.cross(a2, a3)))<self.accuracy:
            return 2
        b1, b2, b3 = b1/mm(a1, b1), b2/mm(a2, b2), b3/mm(a3, b3)

        # Try to rule out wrong a1 and a2 quickly
        if not self.quickCheck(a1, samexyz):
            return 0
        if not self.quickCheck(a2, samexyz):
            return 1

        # Should be integers
        # PBC comensurate with a1...3
        Ints=mm(N.array([b1, b2, b3]), N.transpose(self.pbc))
        if N.max(N.abs(N.round(Ints)-Ints))>self.accuracy:
            return 2

        # Do volume check
        pbcvol = abs(mm(self.pbc[0], N.cross(self.pbc[1], self.pbc[2])))
        vol = abs(mm(a1, N.cross(a2, a3)))
        numcells = pbcvol/vol
        basisatoms = self.NN/(1.0*numcells)

        # numcells and basisatoms should be integers!
        if abs(N.round(numcells)-numcells)>self.accuracy or\
                abs(N.round(basisatoms)-basisatoms)>self.accuracy:
            return 2

        self.NNbasis = N.round(basisatoms)

        # Now check all atoms
        if not self.fullCheck(a1):
            return 0
        if not self.fullCheck(a2):
            return 1
        if not self.fullCheck(a3):
            return 2

        # Success
        return -1

    def quickCheck(self, a, samexyz):
        # Check that first atom repeats over a
        pbc = self.pbc
        bp1, bp2, bp3 = self.bp1, self.bp2, self.bp3

        diff = samexyz[0:]-samexyz[0].reshape((1, 3))

        for ii in range(1, 6):
            ndiff = diff-a*ii # One atom should have ndiff = m_i pbc_i
            indx = N.floor(self.accuracy+mm(ndiff, N.transpose(N.array([bp1, bp2, bp3]))))
            ndiff -= mm(indx, pbc)
            dist = N.sum(ndiff*ndiff, axis=1)
            if min(dist)>self.accuracy:
                return False

        return True

    def fullCheck(self, a):
        # Check that atoms repeats over a
        atomtypes = N.unique(self.snr)
        passed = True
        for atomtype in atomtypes:
            sameatoms = N.argwhere(self.snr==atomtype)
            samexyz = self.xyz[sameatoms]
            samexyz = samexyz.reshape((-1, 3))
            shifted = samexyz+a.reshape((1, 3))

            # Move inside PBC
            samexyz = moveIntoCell(samexyz, self.pbc[0, :], self.pbc[1, :], self.pbc[2, :], self.accuracy)
            shifted = moveIntoCell(shifted, self.pbc[0, :], self.pbc[1, :], self.pbc[2, :], self.accuracy)

            # Should be the same if sorted!
            ipiv = N.lexsort(N.round(N.transpose(samexyz)/self.accuracy))
            samexyz = samexyz[ipiv, :]
            ipiv = N.lexsort(N.round(N.transpose(shifted)/self.accuracy))
            shifted = shifted[ipiv, :]

            if not N.allclose(samexyz, shifted, atol=self.accuracy):
                passed = False

        return passed

    def what(self):
        # Returns directions to calculate bandstructure
        b1, b2, b3= self.b1, self.b2, self.b3
        # Here we adopt the solid-state physics definition
        # a_i b_j=2*pi*delta_ij
        b1, b2, b3= 2*N.pi*b1, 2*N.pi*b2, 2*N.pi*b3
        G = N.zeros(3, N.float)
        # Symmetry k-points from this publication:
        # https://dx.doi.org/10.1016%2Fj.commatsci.2010.05.010
        if self.latticeType == 'CUBIC':
            # Table 2 - Figure 1
            X = b1*0/1+b2*1/2+b3*0/1
            M = b1*1/2+b2*1/2+b3*0/1
            R = b1*1/2+b2*1/2+b3*1/2
            self.path = [[G, 'G'], [X, 'X'], [M, 'M'], [G, 'G'], [R, 'R'], [X, 'X']]
        elif self.latticeType == 'FCC':
            # Table 3 - Figure 2
            X = b1*1/2+b2*0/1+b3*1/2
            K = b1*3/8+b2*3/8+b3*3/4
            L = b1*1/2+b2*1/2+b3*1/2
            U = b1*5/8+b2*1/4+b3*5/8
            W = b1*1/2+b2*1/4+b3*3/4
            self.path = [[G, 'G'], [X, 'X'], [W, 'W'], [K, 'K'], [G, 'G'], [L, 'L'], [U, 'U']]
        elif self.latticeType == 'BCC':
            # Table 4 - Figure 3
            P = b1*1/4+b2*1/4+b3*1/4
            NN= b1*0/1+b2*0/1+b3*1/2
            H = b1*1/2-b2*1/2+b3*1/2
            self.path = [[G, 'G'], [H, 'H'], [P, 'P'], [G, 'G'], [NN, 'N'], [H, 'H']]
        elif self.latticeType == 'HEX':
            # Determine ortogonal array c3
            if N.allclose(N.dot(b1, b3), (0., 0., 0.)):
                if N.allclose(N.dot(b2, b3), (0., 0., 0.)):
                    c1, c2, c3 = b1, b2, b3
                else:
                    c1, c2, c3 = b3, b2, b1
            else:
                c1, c2, c3 = b1, b3, b2
            # Determine angle between c1 and c2
            a = N.arccos(mm(c1, c2)/distance(c1)/distance(c2))*180/N.pi
            if a > 90:
                # angle is 120deg, should be 60deg
                c2 = -c2
#               i.e.:
#               K = c1*1/3-c2*1/3+c3*0/1
#               H = c1*1/3-c2*1/3+c3*1/2
            A = c1*0/1+c2*0/1+c3*1/2
            H = c1*1/3+c2*1/3+c3*1/2
            K = c1*1/3+c2*1/3+c3*0/1
            L = c1*1/2+c2*0/1+c3*1/2
            M = c1*1/2+c2*0/1+c3*0/1
            self.path = [[G, 'G'], [M, 'M'], [K, 'K'], [G, 'G'], [A, 'A'], [L, 'L'], [H, 'H'], [A, 'A']]
        elif self.latticeType == 'TETRAGONAL':
            # Table 5 - Figure 4
            A = b1*1/2+b2*1/2+b3*1/2
            M = b1*1/2+b2*1/2+b3*0/1
            R = b1*0/1+b2*1/2+b3*1/2
            X = b1*0/1+b2*1/2+b3*0/1
            Z = b1*0/1+b2*0/1+b3*1/2
            self.path = [[G, 'G'], [X, 'X'], [M, 'M'], [G, 'G'], [Z, 'Z'], [R, 'R'], [A, 'A'], [Z, 'Z']]
        elif self.latticeType == 'ORTHORHOMBIC':
            X = b1*1/2+b2*0/2+b3*0/2
            Y = b1*0/2+b2*1/2+b3*0/2
            S = b1*1/2+b2*1/2+b3*0/2
            Z = b1*0/2+b2*0/2+b3*1/2
            self.path = [[G, 'G'], [X, 'X'], [S, 'S'], [Y, 'Y'], [G, 'G'], [Z, 'Z']]
        else:
            print "Symmetry: ERROR. Do not know what directions to calculate the phonon bandstructure for lattice %s."%self.latticeType
            sys.exit('Error: unknown lattice type')
        return self.path


###########################################################
# Mathematical helpers

def myUnique(myset, accuracy):
    # Union, sort, remove duplicates and change signs to get positive x

    # Change sign to get possitive x-coordinate or if zero y>0
    changex = myset[:, 0]<-accuracy
    zerox = N.abs(myset[:, 0])<=accuracy
    changey = zerox*myset[:, 1]<-accuracy
    zeroy = N.abs(myset[:, 1])<=accuracy
    changez = zerox*zeroy*myset[:, 2]<-accuracy
    change = -2*((changex+changey+changez)>0)+1
    myset = (N.transpose(N.array([[1], [1], [1]])*change))*myset

    # Sort on z,y,x
    ipiv = N.lexsort(N.round(N.transpose(myset)/accuracy))
    myset = myset[ipiv, :]

    # Remove duplicates
    newset=myset[N.where(N.sum(N.abs(myset[1:, :]-myset[:-1, :]), axis=1)>accuracy)]
    newset = N.concatenate((newset, [myset[-1, :]]))
    myset = newset

    # Sort by length
    dist = N.sum(myset*myset, 1) # Sort by length
    ipiv = dist.argsort()
    myset = myset[ipiv, :]

    # remove zero length
    if abs(distance(myset[0]))<accuracy:
        myset=myset[1:]

    return myset


def myIntersect(x, y, accuracy):
    # Intersection

    if len(x)*len(y)!=0:
        xx, yy, z = N.round(x/accuracy), N.round(y/accuracy), []
        for ii in xx:
            for jj in yy:
                if N.allclose(ii, jj, atol=accuracy):
                    z+=[x[N.where(N.all(ii==xx, 1))[0][0]]]
        if len(z)==0: return N.array([])
        else: return myUnique2(N.array(z), accuracy)
    else:
        return []


def myUnique2(myset, accuracy):
    # Union, sort, remove duplicates

    if len(myset)==0:
        return []
    else:
        # Sort on z,y,x
        ipiv = N.lexsort(N.round(N.transpose(myset)/accuracy))
        myset = myset[ipiv, :]

        # Remove duplicates
        newset=myset[N.where(N.sum(N.abs(myset[1:, :]-myset[:-1, :]), axis=1)>accuracy)]
        newset = N.concatenate((newset, [myset[-1, :]]))
        myset = newset

        # Sort by length
        dist = N.sum(myset*myset, axis=1) # Sort by length
        ipiv = dist.argsort()
        myset = myset[ipiv, :]

        return myset


def distance(x):
    return N.sqrt(N.sum(x*x, axis=-1))


def moveIntoCell(xyz, a1, a2, a3, accuracy):
    b1, b2, b3 = N.cross(a2, a3), N.cross(a1, a3), N.cross(a1, a2)
    b1, b2, b3 = b1/mm(a1, b1), b2/mm(a2, b2), b3/mm(a3, b3)
    copy = xyz.copy()
    indx = N.floor(accuracy+mm(copy, N.transpose(N.array([b1, b2, b3]))))
    copy -= mm(indx, N.array([a1, a2, a3]))
    return copy


def moveIntoClosest(xyz, a1, a2, a3):
    done = False
    while not done:
        done = True
        for ix in range(-2, 3):
            for iy in range(-2, 3):
                for iz in range(-2, 3):
                    indx = N.where(distance(xyz)>distance(xyz+ix*a1+iy*a2+iz*a3))
                    for ii in indx[0]:
                        xyz[ii] = xyz[ii]+ix*a1+iy*a2+iz*a3
                        done = False
    return xyz


def myPermute(plist):
    if len(plist)==0: return[]
    nlist=[]
    low = myPermute(plist[1:])
    for ii in range(len(plist[0])):
        if len(low)>0:
            for jj in low:
                nlist += [[plist[0][ii]]+jj]
        else:
            nlist += [[plist[0][ii]]]
    return nlist


def findRadi(a1, a2, a3):
    # Find sphere that fits in a1..3
    poss = []
    for i1 in range(-1, 2):
        for i2 in range(-1, 2):
            for i3 in range(-1, 2):
                poss += [a1*i1+a2*i2+a3*i3]
    dist = N.sort(distance(N.array(poss)))
    return dist[1]/2.0

if __name__ == '__main__':
    import numpy.linalg as LA

    # Tourture test FCC lattice
    a1=N.array([0, 1, 1])
    a2=N.array([1, 0, 1])
    a3=N.array([1, 1, 0])

    for iitest in range(200000):
        print '\nRunning test no.', iitest
        N1, N2, N3 = int(N.random.rand()*10+1), int(N.random.rand()*10+1), int(N.random.rand()*10+1)
        NB = 2# int(N.random.rand()*3+1)
        basis = (N.random.rand(NB, 3)-.5)*10-6
        basis = N.array([basis[0, :], basis[0, :]+a1/2.0])
        xyz, anr = [], []
        for ix in range(N1):
            for iy in range(N2):
                for iz in range(N3):
                    for ib in range(NB):
                        xyz+=[ix*a1+iy*a2+iz*a3+basis[ib, :]+(N.random.rand(3)-1)*1e-6]
                        anr += [ib]
        u1, u2 = N.random.rand(3)-0.5, N.random.rand(3)-0.5
        u3 = N.cross(u1, u2)
        u1, u3 = u1/LA.norm(u1), u3/LA.norm(u3)
        u2 = N.cross(u1, u3)

        U = N.array([u1, u2, u3])
        xyz = N.array([N.dot(U, x) for x in xyz])
        NN=len(xyz)
        tmp, ipiv = N.arange(NN), []
        for ii in range(NN):
            jj=int(N.random.rand()*len(tmp))
            ipiv += [tmp[jj]]
            tmp = N.hstack((tmp[:jj], tmp[jj+1:]))
        anr, ipiv = N.array(anr), N.array(ipiv)
        sym = Symmetry(accuracy=0.001)
        sym.setupGeom([N.dot(U, a1*N1), N.dot(U, a2*N2), N.dot(U, a3*N3)],
                      anr[ipiv], anr[ipiv], xyz[ipiv], onlyLatticeSym=False)
        if len(sym.pointU33)!=48 or sym.basis.NN!=NB or len(sym.U33)!=8:
            print(N1, N2, N3, NB)
            sys.exit('Failed in tourture test of symmetry')
