version = "SVN $Id$"
print version

import SiestaIO as SIO
import MiscMath as MM
import numpy as N

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
    latticeType : string FCC/BCC/CUBIC/HEX/GRAPHENE/POLYMER/TETRAGONAL
    basis.xyz : xyz of basis inside a1..a3
    basis.snr, .anr, .NN : Siesta/atomic number and number of atoms

    pointU33 : List of lattice point symmetries [?,3,3], i.e., the 
             rotation/mirror symmetry operators in space
    
    # Full set of point operations around rotation center
    fullU33      : List of symmetry operations for lattice + basis, i.e.,
               subset of pointU33 [?,3,3]
    fullOrigo    : List of rotation origin for the operations in fullU33

    # Irreducible symmetry operations for lattice + basis:         
    U33      : List of symmetry operations for lattice + basis, i.e.,
               subset of pointU33 [?,3,3]
    origo    : List of rotation origin for the operations in U33
    rankU33  : Rank, i.e., N; U^N=I

    All arrays are numpy type.

    Use:
       sym = Symmetry.Symmetry(fn='filename.XV',accuracy=1e-4)
       Where accuracy is the fudge factor in Ang.
    """
    def __init__(self,fn=None,accuracy=1e-4,onlyLatticeSym=False):
        self.accuracy=accuracy
        if fn!=None:
            self.readXV(fn,onlyLatticeSym)
        pass
    
    ###########################################################
    # Init

    def readXV(self,fn,onlyLatticeSym=False):
        # Get geometry
        self.pbc, self.snr, self.anr, self.xyz = SIO.ReadXVFile(fn)
        self.pbc, self.snr, self.anr, self.xyz = N.array(self.pbc), N.array(self.snr), N.array(self.anr), N.array(self.xyz)
        self.NN = len(self.anr)
        self.findSymmetry(onlyLatticeSym)

    def setupGeom(self,pbc, snr, anr, xyz,onlyLatticeSym=False):
        self.pbc, self.snr, self.anr, self.xyz = N.array(pbc), N.array(snr), N.array(anr), N.array(xyz)
        self.NN = len(self.anr)
        self.findSymmetry(onlyLatticeSym)
        
    ###########################################################
    # Main symmetry finding

    def findSymmetry(self,onlyLatticeSym=False):
        # Reciprocal lattice for PBC
        pbc = self.pbc
        bp1, bp2, bp3 = N.cross(pbc[1],pbc[2]), N.cross(pbc[0],pbc[2]),\
            N.cross(pbc[0],pbc[1])
        bp1, bp2, bp3 = bp1/mm(pbc[0],bp1), bp2/mm(pbc[1],bp2),\
            bp3/mm(pbc[2],bp3)
        self.bp1, self.bp2, self.bp3 = bp1, bp2, bp3 

        # Find minimal unit cell
        self.findLattice()
        
        # Calculate point group symmetry operators for lattice system 
        self.latticeGroup()
        
        if not onlyLatticeSym:
            # Point group including basis (lower than for lattice system)
            self.pointGroup()

            # Reduce to smallest possible subgroup which generate the full group
            self.findIrreducible()

        return

    ###########################################################
    # Make the dynamical matrix symmetric 

    def symmetrizeFC(self,FC,FCfirst,FClast,radi=0.0):
        NN, NU, NFC = self.NN, len(self.U33), FClast-FCfirst+1
        if self.basis.NN>NFC:
            print "Phonons: ERROR: FCfirst/last do contain all atoms in the basis (%i)."%self.basis.NN
            kuk

        # Rearrange to FC_ia,jb, force from atom j, axis b to atom i, axis a
        FCn = N.zeros((NFC, 3, self.NN, 3))
        for ii in range(NFC):
            for jj in range(3):
                FCn[ii, jj, :, :] = FC[ii*3+jj, :, :]
        FC = FCn
        
        # Rearange basis to fit FCfirst...FClast order in Siesta FC file
        basisxyz = moveIntoCell(self.xyz[FCfirst-1:FClast],\
                                self.a1,self.a2,self.a3,self.accuracy)
        ipiv = []
        for elem in basisxyz[0:FClast-FCfirst+1]:
            for jj, cmp in enumerate(self.basis.xyz):
                if N.allclose(elem,cmp,atol=self.accuracy):
                    ipiv+=[jj]
                    break
        if len(N.unique(ipiv))!=self.basis.NN:
            print "Symmetry: Some confusion about the order of the basis of atoms and the FCfirst/last region."
            kuk
        self.basis.xyz = self.basis.xyz[ipiv,:]
        self.basis.snr, self.basis.anr = N.array(self.basis.snr)[ipiv], N.array(self.basis.anr)[ipiv]

        # Find out which basis atom corresponds to each atom
        xyz = moveIntoCell(self.xyz,self.a1,self.a2,self.a3,self.accuracy)
        self.basisatom = N.zeros((self.NN))
        for ii in range(self.basis.NN):
            indx = N.where(N.sum(N.abs(xyz-self.basis.xyz[ii,:]), axis=1)<self.accuracy)
            self.basisatom[indx[0]]=ii

        # Symmetry operations are complicated by the periodic structure 
        # in the Siesta calculation. So ... limit range of interaction.

        # Find radie and center for biggest sphere fitting into pbc
        if radi==0.0:
            radi = findRadi(self.pbc[0],self.pbc[1],self.pbc[2])*\
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
        PL  = N.zeros((NU, NFC, NN, NN),N.int)
        PR  = N.zeros((NU, NFC),N.int)
        for iU in range(NU):
            for ii in range(NFC):
                SIO.printDone(iU*NFC+ii, NU*NFC,'Symmetrizing')

                # Figure out which atoms are connected by symmetry op.
                xyz = self.xyz.copy()-self.origo[iU]
                FCxyz = xyz[FCfirst-1+ii, :].copy()
                xyz = moveIntoClosest(N.array(xyz)-FCxyz,\
                                    self.pbc[0],self.pbc[1],self.pbc[2])+FCxyz
                nxyz = N.transpose(mm(self.U33[iU],N.transpose(xyz)))
                # Did the moving atom move between unit cells?
                nFCxyz, iix = nxyz[FCfirst-1+ii, :].copy(),  10
                for ix in range(-2, 3): 
                    for iy in range(-2, 3): 
                        for iz in range(-2, 3):
                            if N.any(N.sum(N.abs(nFCxyz+ix*self.a1+iy*self.a2+iz*self.a3-xyz[FCfirst-1:FClast]), axis=1)<self.accuracy):
                                iix, iiy, iiz = ix, iy, iz
                if iix==10: kuk
                tFCxyz = nFCxyz+iix*self.a1+iiy*self.a2+iiz*self.a3
                shift = tFCxyz-nFCxyz # Shift all new coordinates
                nxyz = moveIntoClosest(N.array(nxyz)+shift-FCxyz,\
                                    self.pbc[0],self.pbc[1],self.pbc[2])+FCxyz
    
                # Which atoms are within radi?
                indx = N.where(distance(xyz-FCxyz)<radi)[0]
    
                # Find the target atom
                diff = N.sum(N.abs(nxyz[indx].reshape((-1,1,3))-\
                         xyz.reshape((1,-1,3))),axis=2)<self.accuracy
                tindx = N.where(diff)
                tindx2 = tindx[1]
                indx3= indx[tindx[0]]
                if len(indx3)!=len(indx):
                    for ix in range(-1,2):
                        for iy in range(-1,2):
                            for iz in range(-1,2):
                                if ix!=0 or iy!=0 or iz!=0:
                                    tindxs = N.where(N.sum(N.abs(nxyz[indx].reshape((-1,1,3))+\
                                                                  self.pbc[0]*ix+self.pbc[1]*iy+self.pbc[2]*iz-\
                                                                  xyz.reshape((1,-1,3))),axis=2)<self.accuracy)
                                    indx3 = N.concatenate((indx3, indx[tindxs[0]]))
                                    tindx2 = N.concatenate((tindx2, tindxs[1]))
                if len(indx3)!=len(indx):
                    print "kuk"
                # Make permutation matrix
                PL[iU,ii, tindx2,indx3] = 1
                indx2 = N.where(distance(xyz-tFCxyz)<self.accuracy)
                PR[iU, ii] = self.basisatom[indx2[0]]

                # TF: Check that the symmetry operations apply to FC??
                for kk in range(len(indx3)):
                    oFC = FC[ii,:,indx3[kk],:].reshape((3,3))
                    sFC = mm(UL[iU],FC[indx2[0],:,tindx2[kk],:].reshape((3,3)),UR[iU])
                    if N.max(N.abs(oFC-sFC))>0.5:
                        print oFC
                        print sFC
        
        def applySym(FC):
            FCn = N.tensordot(UR[iU], FC, ((1, 1)))
            FCn = N.swapaxes(FCn, 0, 1)
            FCn = N.tensordot(FCn, UL[iU], ((3, 0)))
            FCn2=0*FCn
            for ii in range(NFC):   
                tmp = N.tensordot(PL[iU,ii, :, :], \
                                      FCn[ii, :, :, :], ((1, 1)))
                tmp = N.swapaxes(tmp, 0, 1)
                FCn2[PR[iU, ii], :, :, :] = tmp
            return FCn2

        # Symmetrize dynamical matrix 
        print "Symmetry: Iterative application of symmetries"
        iter, change = 0, 10
        FCo = FC.copy()
        # Uncomment the two lines below to skip the application of symmetries
        #iter, change = 10, 10
        #FCs = FC.copy()
        while iter<10 and change>1e-10:
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
        FC = N.zeros((NFC*3, self.NN, 3))
        for ii in range(NFC):
            for jj in range(3):
                FC[ii*3+jj, :, :] = FCs[ii, jj, :, :] 
        
        return FC

    ###########################################################
    # Find irreducible symmetry operations and print
    # 
    def findIrreducible(self):
        # Just print for the lattice group
        ipiv, sign, rotn = self.reduce(self.pointU33)

        # Go through the different rotation centers
        lasto, nU, nO, nR, nS = 0, [], [], [], []
        for ii in range(1,len(self.origo)):
            if not N.allclose(self.origo[ii-1],self.origo[ii],self.accuracy):
                ipiv, sign, rotn = self.reduce(N.concatenate((self.U33[lasto:ii], [N.eye(3)])))
                if len(rotn)>1 or rotn[0]!=1: 
                    for jj, kk in enumerate(lasto+N.array(ipiv)):
                        nU += [self.U33[kk]]
                        nO += [self.origo[kk]]
                        nR += [rotn[jj]]
                    print "Symmetry: Irreducible point operations around %f %f %f"%(self.origo[ii-1][0],self.origo[ii-1][1],self.origo[ii-1][2])
                    print "Ranks*determinant : ",rotn*sign
                lasto = ii

        if len(self.U33)!=0:
            ipiv, sign, rotn = self.reduce(N.concatenate((self.U33[lasto:], [N.eye(3)])))
            for jj, kk in enumerate(lasto+N.array(ipiv)):
                nU += [self.U33[kk]]
                nO += [self.origo[kk]]
                nR += [rotn[jj]]

            print "Symmetry: Irreducible point operations around %f %f %f"%(self.origo[ii-1][0],self.origo[ii-1][1],self.origo[ii-1][2])
            print "Ranks*determinant : ",rotn*sign

        self.fullU33, self.fullOrigo = self.U33, self.origo
        self.U33, self.origo, self.rankU33 = nU, nO, nR
        return

    def reduce(self, Ulist):
        # Find least number of symmetry operations that generate the whole group
        keep = [ii for ii in N.array(Ulist).copy()]
        ChT = N.zeros((len(keep),len(keep)),N.int) # Character table
        for ii in range(len(keep)):
            for jj in range(len(keep)):
                tmp = mm(keep[ii],keep[jj])
                indx=N.where(N.sum(N.sum(N.abs(keep-tmp), axis=2), axis=1)<self.accuracy)
                ChT[ii,jj]=indx[0][0]

        # Calculate rank
        rotn = []
        for ii in keep:
            a, M = ii, 1
            while not N.allclose(N.eye(3),a,atol=self.accuracy):
                a, M =mm(a,ii), M+1
            rotn+=[M]
        fullRank = rotn[:]
        
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
                x= N.unique(N.concatenate((x,pick(ChT, x))))
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
            while not N.allclose(N.eye(3),a,atol=self.accuracy):
                a, M =mm(a,ii), M+1
            rotn+=[M]
        
        return basis, sign, rotn

    ###########################################################
    # point group of lattice+basis

    def pointGroup(self):
        # Check which of the point operations that leave the lattice unchanged
        # Guessing at usefull origins for the rotation
        
        atomtypes, Ulist, origo = N.unique(self.basis.snr), [], []

        # Generate possible rotation centers ...
        # Hope that midpoint between 2 or 3 atoms will suffice
        
        basisxyz = self.basis.xyz
        poss = [basisxyz[ii,:] for ii in range(self.basis.NN)]
        if self.basis.NN>1:
            for ii in basisxyz:
                for jj in basisxyz:
                    poss += [(ii+jj)/2.0]
        if self.basis.NN>2:
            for ii in basisxyz:
                for jj in basisxyz:
                    for kk in basisxyz:
                        poss += [(ii+jj+kk)/3.0]

        # Remove duplicates
        poss = myUnique2(N.array(poss),self.accuracy)
        poss = moveIntoCell(poss,  self.a1,  self.a2,  self.a3, 0)
        poss = myUnique2(N.array(poss),self.accuracy)

        # Loop over possible origins
        for iposs, x0 in enumerate([poss[ii,:] for ii in range(len(poss))]):
            for iU, U in enumerate(self.pointU33):
                SIO.printDone(iposs*len(self.pointU33)+iU,\
                                  len(self.pointU33)*len(poss),\
                                  'Checking basis symmetries')
                passed=True
                for atomtype in atomtypes:
                    xyz = self.basis.xyz[N.argwhere(self.basis.snr==atomtype)]
                    xyz = xyz.reshape((-1, 3))-x0
                    # Check operations
                    xyz =moveIntoCell(xyz,self.a1,self.a2,self.a3,self.accuracy)
                    ipiv = N.lexsort(N.round(N.transpose(xyz)/self.accuracy))
                    xyz = xyz[ipiv,:]
                    
                    res = N.transpose(mm(N.transpose(U), N.transpose(xyz)))
                    res =moveIntoCell(res,self.a1,self.a2,self.a3,self.accuracy)
                    ipiv = N.lexsort(N.round(N.transpose(res)/self.accuracy))
                    res = res[ipiv,:]
                    if not N.allclose(xyz,res,atol=self.accuracy):
                        passed = False 
                        break
                if passed and not N.allclose(U, N.eye(3), atol=self.accuracy):
                    Ulist += [U]
                    origo  += [x0]
        self.U33, self.origo = Ulist, origo
        print "Symmetry: %i point symmetry operations (with rotation centers) found for lattice+basis (some are unnecessary \"duplicates\")."%len(Ulist)

        # Calculate rank of operations
        self.rankU33, sign = [], []
        for U in self.U33:
            tmp, ii = U, 1
            while ii<7 and not N.allclose(N.eye(3),tmp,self.accuracy):
                tmp, ii = mm(tmp,U), ii+1
            if ii > 6: kuk
            self.rankU33 += [ii]
            sign += [N.linalg.det(U)]
        print "Symmetry: rank*det"
        print N.array(self.rankU33)*sign
        return


    ###########################################################
    # point group of lattice

    def latticeGroup(self):
        # Generate all unitary matrices that take a1...a3 to 
        #   lattice points
        
        a1, a2, a3= self.a1, self.a2, self.a3
        b1, b2, b3= self.b1, self.b2, self.b3
        points, dist = [], []        
        # Generate lattice points
        for i1 in range(-2,3): # -1,0,1 enough ? use 2 for safety
            for i2 in range(-2,3): 
                for i3 in range(-2,3):
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
            U = ii[0].reshape(3,1)*b1.reshape(1,3)+\
                ii[1].reshape(3,1)*b2.reshape(1,3)+\
                ii[2].reshape(3,1)*b3.reshape(1,3)
            if N.allclose(mm(U,N.transpose(U)),N.eye(3), self.accuracy):
                Ulist+=[U]

        self.pointU33 = Ulist
        print "Symmetry: %i point symmetry operations found for lattice"%len(Ulist)
        sign = N.array([N.linalg.det(ii) for ii in Ulist])
        
        rotn = []
        for ii in Ulist:
            a, M = ii, 1
            while not N.allclose(N.eye(3),a,atol=self.accuracy):
                a, M =mm(a,ii), M+1
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

        # Brute force ... find vectors connecting same type of atoms, 
        # try and see if they are lattice vectors starting with shortest
        sameatoms = N.argwhere(self.snr==self.snr[0])
        samexyz = self.xyz[sameatoms]
        samexyz = samexyz.reshape((-1, 3))
        samexyz = moveIntoCell(samexyz,self.pbc[0,:],self.pbc[1,:],\
                                   self.pbc[2,:],self.accuracy)
        # Make matrix of difference between coordinates
        possible = samexyz.reshape((1,-1,3))-samexyz.reshape((-1,1,3))
        # Reshape into list of vectors
        possible = possible.reshape((-1,3))
        # Add the periodic boundary conditions in case the unit 
        # cell is the whole cell.
        possible = N.concatenate((possible,self.pbc )) 
        # remove duplicates etc
        possible = myUnique(possible,self.accuracy) 

        # Go through combinations and check if possible
        i1, i2, i3, done, NP = 0, 1, 2, False, len(possible)

        def increase(i1,i2,i3,NP,incindx):
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
            i2 = max(i1+1,i2)
            i3 = max(i2+1,i3)
            #if incindx<=1: print i1, i2, i3
            return i1,i2,i3
        
        # Step through possibilities
        while i1<NP-2 and not done:
            a1, a2 = possible[i1], possible[i2]
            if distance(N.cross(a1,a2))>self.accuracy: # Non-parallell
                a3 = possible[i3]
                incindx = self.checkLattice(a1,a2,a3,samexyz)
                if incindx!=-1:
                    i1, i2, i3 = increase(i1,i2,i3,NP,incindx)
                else:
                    done = True
            else:
                i1, i2, i3 = increase(i1,i2,i3,NP,1)
        
        if not done:
            # Should not be necessary since PBC is added to possible!
            print "Symmerty: WARNING Could not find smaller unitcell"
            print "  ...   : Continuing with Siesta periodic cell"
            a1, a2, a3 = self.pbc[0,:], self.pbc[1,:], self.pbc[2,:]
            self.NNbasis = self.NN
            kuk

        # Rearange to fit with PBC
        a1, a2, a3 = self.makeHumanReadable(a1,a2,a3)

        b1, b2, b3 = N.cross(a2,a3), N.cross(a1,a3), N.cross(a1,a2)
        b1, b2, b3 = b1/mm(a1,b1), b2/mm(a2,b2), b3/mm(a3,b3)

        self.a1, self.a2, self.a3 = a1, a2, a3
        self.b1, self.b2, self.b3 = b1, b2, b3

        indx = N.round(abs(mm(self.pbc,N.transpose(N.array([b1,b2,b3])))))
        N1,  N2,  N3 =indx[0, 0],  indx[1, 1],  indx[2, 2]
        print "Symmetry: Lattice structure"
        print "%i atoms in the basis"%(self.NNbasis)
        print "a1 = (%f,%f,%f), N1=%i"%(a1[0],a1[1],a1[2],N1)
        print "a2 = (%f,%f,%f), N2=%i"%(a2[0],a2[1],a2[2],N2)
        print "a3 = (%f,%f,%f), N3=%i"%(a3[0],a3[1],a3[2],N3)

        # Shift to more convenient origin
        xyz = moveIntoCell(self.xyz, a1, a2, a3, self.accuracy)
        tmp=mm(N.array([b1,b2,b3]),N.transpose(xyz))
        mx, my, mz = min(tmp[0,:]),min(tmp[1,:]),min(tmp[2,:])
        shift = mx*a1+my*a2+mz*a3
        print "Shifting coordinates by : ",-shift
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
            kuk
        self.basis = basis

        # Determine types of lattices from the angles / lengths
        # Sort on length
        length = N.array([distance(ii) for ii in [a1,a2,a3]])
        ipiv = N.argsort(length)
        length = N.array(length)[ipiv]
        a123 = N.array([a1,a2,a3])[ipiv]

        # Calculate angles in degrees
        angles = []
        for ii,a in enumerate([a123[0],a123[1]]):
            for jj,b in enumerate([a123[1],a123[2]][ii:]):
                angles += [N.round(N.arccos(mm(a,b)/distance(a)/distance(b))*\
                                   180/N.pi)]
        # Move angels into range [0,90] 
        angles = [min(ii,180-ii) for ii in angles] 
        print 'Lattice angles  =', angles
        print 'Lattice lengths =', length
        latticetype = 'UNKNOWN'
        if N.allclose(length[0],length[1],atol=self.accuracy) and\
                N.allclose(length[1],length[2],atol=self.accuracy):
            # All lengths same
            if N.sum(N.abs(N.array(angles)-90))<2:
                latticetype = 'CUBIC'
            if N.sum(N.abs(N.array(angles)-60))<2:
                latticetype = 'FCC'
            if N.sum(N.abs(N.array(angles)-70.5))<4:
                latticetype = 'BCC'
        elif N.allclose(length[0],length[1],atol=self.accuracy) and \
                N.sum(N.abs(N.array(angles)-90))<2:
            latticetype = 'TETRAGONAL'
        elif N.allclose(length[1],length[2],atol=self.accuracy):
            if N.sum(N.abs(N.array(angles)-90))<2:
                latticetype = 'POLYMER'
        elif N.allclose(length[0],length[1],atol=self.accuracy):
            if N.abs(angles[0]-60)<2:
                latticetype="GRAPHENE"
        print "Symmetry: Lattice  = %s"%latticetype
        self.latticeType = latticetype

    def makeHumanReadable(self,a1,a2,a3):
        # TODO! Does not give "Human" readbility ...
        # Choose vectors that give the smallest angles between a1..a3

        # Sort on length
        ipiv = N.argsort(distance(N.array([a1,a2,a3])))
        na1 = [a1,a2,a3][ipiv[0]]
        na2 = [a1,a2,a3][ipiv[1]]
        na3 = [a1,a2,a3][ipiv[2]]
        a1,  a2,  a3 = na1,  na2,  na3
        
        # Make possible n a1 + m a2 + l a3 
        poss = []
        for i1 in range(-2,3):
            for i2 in range(-2,3):
                for i3 in range(-2,3):
                    poss += [a1*i1+a2*i2+a3*i3]
        poss = N.array(poss)

        # possiblities for a1, a2, a3 based on length
        poss1 = poss[N.where(N.abs(distance(poss)-distance(a1))<self.accuracy)]
        poss2 = poss[N.where(N.abs(distance(poss)-distance(a2))<self.accuracy)]
        poss3 = poss[N.where(N.abs(distance(poss)-distance(a3))<self.accuracy)]
        poss = N.concatenate((poss1,poss2,poss3)) # Keep duplicates
        
        # Choose the one fitting the PBC best 
        vec = self.pbc[0,:]
        scal = mm(vec,N.transpose(poss))
        a1=poss[N.where(N.sort(scal)[-1]==scal)[0][0]] 
        
        # Remove linear dependent
        poss=poss[N.where(N.abs(N.abs(mm(poss,a1))-\
                                    N.abs(distance(poss)*distance(a1)))>self.accuracy)]

        # Choose the one fitting the PBC best 
        vec = self.pbc[1,:]
        scal = mm(vec,N.transpose(poss))
        a2=poss[N.where(N.sort(scal)[-1]==scal)[0][0]] 

        # Remove linear dependent
        poss=poss[N.where(N.abs(N.abs(mm(poss,a2))-\
                                    N.abs(distance(poss)*distance(a2)))>self.accuracy)]

        # Choose the one fitting the PBC best 
        vec = self.pbc[2,:]
        scal = mm(vec,N.transpose(poss))
        a3=poss[N.where(N.sort(scal)[-1]==scal)[0][0]] 

        if N.linalg.det([a1,a2,a3])<-self.accuracy:
            a3=-a3

        return a1,a2,a3


    def checkLattice(self,a1,a2,a3,samexyz):
        # Check that a1..3 span lattice, return -1 for success 
        #  otherwise index to increase

        # Test of PBC
        b1, b2, b3 = N.cross(a2,a3), N.cross(a1,a3), N.cross(a1,a2)
        if distance(b1)<self.accuracy or distance(b2)<self.accuracy or\
                abs(mm(a1,N.cross(a2,a3)))<self.accuracy:
            return 2
        b1, b2, b3 = b1/mm(a1,b1), b2/mm(a2,b2), b3/mm(a3,b3)

        # Try to rule out wrong a1 and a2 quickly
        if not self.quickCheck(a1,samexyz):
            return 0
        if not self.quickCheck(a2,samexyz):
            return 1
        
        # Should be integers
        # PBC comensurate with a1...3
        Ints=mm(N.array([b1,b2,b3]),N.transpose(self.pbc))
        if N.max(N.abs(N.round(Ints)-Ints))>self.accuracy:
            return 2 

        # Do volume check
        pbcvol = abs(mm(self.pbc[0],N.cross(self.pbc[1],self.pbc[2])))
        vol = abs(mm(a1,N.cross(a2,a3)))
        numcells = pbcvol/vol
        basisatoms = self.NN/numcells
        
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

    def quickCheck(self,a,samexyz):
        # Check that first atom repeats over a
        pbc = self.pbc
        bp1, bp2, bp3 = self.bp1, self.bp2, self.bp3 

        diff = samexyz[0:]-samexyz[0].reshape((1,3))
        
        for ii in range(1,6):
            ndiff = diff-a*ii # One atom should have ndiff = m_i pbc_i
            indx = N.floor(self.accuracy+mm(ndiff,N.transpose(N.array([bp1, bp2, bp3]))))
            ndiff -= mm(indx,pbc)
            dist = N.sum(ndiff*ndiff, axis=1)
            if min(dist)>self.accuracy:
                return False

        return True

    def fullCheck(self,a):
        # Check that atoms repeats over a
        bp1, bp2, bp3 = self.bp1, self.bp2, self.bp3 
        atomtypes = N.unique(self.snr)
        passed = True
        for atomtype in atomtypes:
            sameatoms = N.argwhere(self.snr==atomtype)
            samexyz = self.xyz[sameatoms]
            samexyz = samexyz.reshape((-1, 3))
            shifted = samexyz+a.reshape((1,3))
            
            # Move inside PBC
            samexyz = moveIntoCell(samexyz,self.pbc[0,:],self.pbc[1,:],self.pbc[2,:],self.accuracy)
            shifted = moveIntoCell(shifted,self.pbc[0,:],self.pbc[1,:],self.pbc[2,:],self.accuracy)

            # Should be the same if sorted!
            ipiv = N.lexsort(N.round(N.transpose(samexyz)/self.accuracy))
            samexyz = samexyz[ipiv,:]
            ipiv = N.lexsort(N.round(N.transpose(shifted)/self.accuracy))
            shifted = shifted[ipiv,:]
            
            if not N.allclose(samexyz,shifted,atol=self.accuracy):
                passed = False
            
        return passed

    def what(self):
        # Returns directions to calculate bandstructure
        # [["text",k-direction, k-origin, Num k-points], ...]
        a1, a2, a3= self.a1, self.a2, self.a3
        b1, b2, b3= self.b1, self.b2, self.b3
        signs = [N.sign(mm(b1,[b1,b2,b3][ii])) for ii in range(3)]
        if self.latticeType == 'FCC':
            X1 = (b1-signs[1]*b2)/2
            X2 = (b1-signs[2]*b3)/2
            K = X1+X2
            L = b1
            what = [['000-100',X1,0*X1,101],
                    ['100-110',X2,X1,101],
                    ['110-000',-K,K,101],
                    ['000-111',0.5*L,0*L,51]]
        elif self.latticeType == 'BCC':
            H = abs(b1-signs[1]*b2-signs[2]*b3)/2
            NN = abs(b1-signs[1]*b2)/2
            P = abs(b1+signs[1]*b2+signs[2]*b3)/4
            G = 0*P
            what = [['G-H',H-G,G,101],
                    ['H-N',NN-H,H,101],
                    ['N-G',G-NN,NN,101],
                    ['G-P',P-G,G,101],
                    ['P-H',H-P,P,101]]
        elif self.latticeType == 'CUBIC':
            X = b1
            L = b1+b2+b3
            K = b1+b2
            what = [['000-100',X,0*X,101],
                    ['100-110',K-X,X,101],
                    ['110-000',-K,K,101],
                    ['000-111',L,0*L,101]]
        elif self.latticeType == 'POLYMER':
            ipiv = N.argsort(distance(N.array([a1,a2,a3])))
            X = [b1,b2,b3][int(ipiv[0])]
            what = [['000-100',X*0.5,0*X,101]]
        elif self.latticeType == 'GRAPHENE':
            #M = (b1+b2)/2
            #K = M+(b1-b2)/6
            M = b1/2.
            K = (b1+b2)/3.
            G = 0*M
            what = [['G-M',M,G,101],
                    ['M-K',K-M,M,101],
                    ['K-G',G-K,K,101]]
        elif self.latticeType == 'TETRAGONAL':
            ipiv = N.argsort(distance(N.array([a1,a2,a3])))
            bb1, bb2, bb3 = [b1,b2,b3][int(ipiv[0])], [b1,b2,b3][int(ipiv[1])], [b1,b2,b3][int(ipiv[2])]
            X = bb1/2
            M = (bb1+bb2)/2
            R = bb3/2
            G = 0*bb1
            what = [['X-G',G-X,X,101],
                    ['G-M',M-G,G,101],
                    ['M-X',X-M,M,101],
                    ['X-G',G-X,X,101],
                    ['G-R',R-G,G,101]]
        elif self.latticeType == 'FCT':
            ipiv = N.argsort(distance(N.array([a1,a2,a3])))
            bb1, bb2, bb3 = [b1,b2,b3][int(ipiv[0])], [b1,b2,b3][int(ipiv[1])], [b1,b2,b3][int(ipiv[2])]
            X = (bb1+bb3)/2
            Z = bb1+(bb2+bb3)/2
            L = (bb1+bb2+bb3)/2
            XP = (bb2+bb3)/2
            ZP = bb2+(bb1+bb3)/2
            G = 0*bb1
            what = [['G-X',X-G,G,101],
                    ['X-Z',Z-X,X,101],
                    ['Z-G',G-Z,Z,101],
                    ['G-Z\'',ZP-G,G,101],
                    ['Z\'-X\'',XP-ZP,ZP,101],
                    ['X\'-G',G-XP,XP,101],
                    ['G-L',L-G,G,101]]
        else:
            print "Symmetry: ERROR. Do not know what directions to calculate the phonon bandstructure for lattice %s."%self.latticeType
            kuk
        print "Symmetry: calculate along k-directions"
        for elem in what:
            print elem[0]
            print "From ",N.round(elem[2]*1e4)/1e4," to ",N.round(1e4*(elem[2]+elem[1]))/1e4 

        return what

###########################################################
# Mathematical helpers

def myUnique(set,accuracy):
    # Union, sort, remove duplicates and change signs to get positive x
    
    # Change sign to get possitive x-coordinate or if zero y>0
    changex = set[:,0]<-accuracy
    zerox = N.abs(set[:,0])<=accuracy
    changey = zerox*set[:,1]<-accuracy
    zeroy = N.abs(set[:,1])<=accuracy
    changez = zerox*zeroy*set[:,2]<-accuracy
    change = -2*((changex+changey+changez)>0)+1
    set = (N.transpose(N.array([[1],[1],[1]])*change))*set

    # Sort on z,y,x
    ipiv = N.lexsort(N.round(N.transpose(set)/accuracy))
    set = set[ipiv,:]

    # Remove duplicates
    newset=set[N.where(N.sum(N.abs(set[1:, :]-set[:-1, :]),axis=1)>accuracy)]
    newset = N.concatenate((newset, [set[-1, :]]))
    set = newset
    
    # Sort by length
    dist = N.sum(set*set,1) # Sort by length
    ipiv = dist.argsort()
    set = set[ipiv, :]

    # remove zero length
    if abs(distance(set[0]))<accuracy:
        set=set[1:]

    return set


def myUnique2(set,accuracy):
    # Union, sort, remove duplicates
    
    # Sort on z,y,x
    ipiv = N.lexsort(N.round(N.transpose(set)/accuracy))
    set = set[ipiv,:]

    # Remove duplicates
    newset=set[N.where(N.sum(N.abs(set[1:, :]-set[:-1, :]),axis=1)>accuracy)]
    newset = N.concatenate((newset, [set[-1, :]]))
    set = newset
    
    # Sort by length
    dist = N.sum(set*set,axis=1) # Sort by length
    ipiv = dist.argsort()
    set = set[ipiv, :]

    return set

def distance(x):
    return N.sqrt(N.sum(x*x, axis=-1))


def moveIntoCell(xyz,a1,a2,a3,accuracy):
    b1, b2, b3 = N.cross(a2,a3), N.cross(a1,a3), N.cross(a1,a2)
    b1, b2, b3 = b1/mm(a1,b1), b2/mm(a2,b2), b3/mm(a3,b3)
    copy = xyz.copy()
    indx = N.floor(accuracy+mm(copy,N.transpose(N.array([b1, b2, b3]))))
    copy -= mm(indx,N.array([a1, a2, a3]))
    return copy

def moveIntoClosest(xyz,a1,a2,a3):
    done = False
    while not done:
        done = True
        for ix in range(-2,3):
            for iy in range(-2,3):
                for iz in range(-2,3):
                    indx = N.where(distance(xyz)>distance(xyz+ix*a1+iy*a2+iz*a3))
                    for ii in indx[0]:
                        xyz[ii] = xyz[ii]+ix*a1+iy*a2+iz*a3
                        done = False
    return xyz

def unique(list):
    list.sort()
    uni=[list[0]]
    for ii in list[1:]:
        if not N.allclose(uni[-1],ii):
            uni+=[ii]
            print ii
    return N.array(uni)

def myPermute(list):
    if len(list)==0: return[]
    nlist=[]
    low = myPermute(list[1:])
    for ii in range(len(list[0])):
        if len(low)>0:
            for jj in low:
                nlist += [[list[0][ii]]+jj]
        else:
            nlist += [[list[0][ii]]]
    return nlist

def findRadi(a1, a2, a3):
    # Find sphere that fits in a1..3
    poss = []
    for i1 in range(-1,2):
        for i2 in range(-1,2):
            for i3 in range(-1,2):
                poss += [a1*i1+a2*i2+a3*i3]
    dist = N.sort(distance(N.array(poss)))
    return dist[1]/2.0

if __name__ == '__main__':
    b=Symmetry('/home/mpn/symmetry/test/test.XV')
    a=Symmetry('/home/mpn/symmetry/NaCl/NaCl.XV')
    c=Symmetry('/home/mpn/symmetry/FCC/FCC.XV')
