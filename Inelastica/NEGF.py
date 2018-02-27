"""

NEGF (:mod:`Inelastica.NEGF`)
=============================

Non-Equilibrium Green's function

.. currentmodule:: Inelastica.NEGF

classes
-------

.. autosummary::
   :toctree:

   SigDir
   SavedSigClass
   ElectrodeSelfEnergy
   GF

"""

import Inelastica.io.siesta as SIO
import Inelastica.MiscMath as MM
import numpy as N
import numpy.linalg as LA
import string
import pickle
import hashlib
import glob
import time
import os
import netCDF4 as NC4
import Inelastica.ValueCheck as VC

# For speed some routines can be linked as F90 code
try:
    import Inelastica.fortran.F90_lapack as F90
    F90_lapack_imp = True
except:
    F90_lapack_imp = False
    print "########################################################"
    print "Problems encountered with F90_lapack.so"
    print "The linking/finding of LAPACK routines does not work"
    print "Ensure the placement of LAPACK in your LD_LIBRARY_PATH"
    print "Try compiling manually following these steps:"
    print " $ cd Inelastica/fortran"
    print " $ source compile.bat (or compile_alternative.bat)"
    print " $ cp F90_lapack.so <python>/site-packages/Inelastica/fortran/"
    print "########################################################"

#try:
#    import scipy.linalg as SLA
#    hasSciPy = True
#except:
#    hasSciPy = False

########################################################


def AssertReal(x, label):
    VC.Check("zero-imaginary-part", abs(x.imag),
             "Imaginary part too large in quantity %s"%label)
    return x.real


def myHash(data):
    return hashlib.md5(pickle.dumps(data)).hexdigest()


def hash2dec(hash):
    return [int(ii, 16) for ii in hash]


def dec2hash(dec):
    s=''
    for ii in range(len(dec)):
        s=s+hex(dec[ii])[2:]
    return s


class SigDir:
    def __init__(self, path):
        self.path, self.data = path, {}
        self.files, self.newFile = [], None
        for ii in glob.glob(path+'/Sig*.nc'):
            print ii
            self.add(ii)

    def add(self, fn):
        ncfile = NC4.Dataset(fn, 'r')
        if 'Done' in ncfile.variables:
            print "Read ", fn
            ncv = ncfile.variables
            hash, hash2l, reE, imE, kp, LR, ispin, etaLead = ncv['hash'][:], ncv['hash2'][:], ncv['reE'][:], ncv['imE'][:], ncv['kp'][:], ncv['left'][:], ncv['ispin'][:], ncv['etaLead'][:]
            hash = [dec2hash(ii) for ii in hash]
            for ii in range(len(reE)):
                hash2 = dec2hash(hash2l[ii])
                print "Found ", hash2
                self.data[hash2] = [ncfile, ii]
            self.files+=[ncfile]
        else:
            ncfile.close()

    def getSig(self, hash, ee, kp, left, ispin, etaLead):
        if left:
            left=1
        else:
            left=0
        hash2 = myHash([N.array(hash), N.array(ee.real), N.array(ee.imag), N.array(kp), N.array(left), N.array(ispin), N.array(etaLead)])
        print "Get ", hash2
        if hash2 in self.data:
            print "Found"
            ncf, ii = self.data[hash2]
            return True, ncf.variables['reSig'][ii, :, :]+1j*ncf.variables['imSig'][ii, :, :]
        else:
            return False, None

    def addSig(self, hash, ee, kp, left, ispin, etaLead, Sig):
        #print "Add ",hash,ee,kp,left,ispin,etaLead,Sig
        if left:
            left=1
        else:
            left=0
        if self.newFile == None:
            self.newFileIndx = 0
            nf = NC4.Dataset(self.path+'/Sig_'+str(N.floor(N.random.rand()*1e9))+'.nc', 'w')
            nf.createDimension('One', 1)
            nf.createDimension('Two', 2)
            nf.createDimension('32', 32)
            nf.createDimension('Dim', len(Sig))
            nf.createDimension('List', None)
            nf.createVariable('kp', 'd', ('List', 'Two'))
            nf.createVariable('reE', 'd', ('List',))
            nf.createVariable('imE', 'd', ('List',))
            nf.createVariable('left', 'i', ('List',))
            nf.createVariable('hash', 'i', ('List', '32'))
            nf.createVariable('hash2', 'i', ('List', '32'))
            nf.createVariable('ispin', 'i', ('List',))
            nf.createVariable('etaLead', 'd', ('List',))
            nf.createVariable('reSig', 'd', ('List', 'Dim', 'Dim'))
            nf.createVariable('imSig', 'd', ('List', 'Dim', 'Dim'))
            self.newFile = nf
            self.files+=[nf]
        NN=self.newFileIndx
        hash2=myHash([N.array(hash), N.array(ee.real), N.array(ee.imag), N.array(kp), N.array(left), N.array(ispin), N.array(etaLead)])
        nfv=self.newFile.variables
        nfv['hash'][NN, :], nfv['reE'][NN], nfv['imE'][NN] = hash2dec(hash), ee.real, ee.imag
        nfv['hash2'][NN, :]= hash2dec(hash2)
        nfv['kp'][NN, :], nfv['left'][NN], nfv['ispin'][NN], nfv['etaLead'][NN] = kp, left, ispin, etaLead
        nfv['reSig'][NN, :, :], nfv['imSig'][NN, :, :] = Sig.real, Sig.imag
        print "Put ", hash2
        self.data[hash2] = [self.newFile, NN]
        self.newFileIndx+=1

    def close(self):
        if not self.newFile==None:
            var=self.newFile.createVariable('Done', 'i', ('One',))
        for f in self.files:
            f.close()


class SavedSigClass:
    """
    Saves calculated Sig in files in the directory of the TSHS file for the electrode.
    1: Each process opens a new file if it needs to write Sigma
    2: The file cannot be read before the finished flag is set
    3: The integrity of the data is maintained by a hash of HS, NA1, NA2, voltage
    """

    def __init__(self):
        self.sigs={}

    def add_hsfile(self, path):
        if not path in self.sigs:
            self.sigs[path]=SigDir(path)

    def getSig(self, path, hash, ee, kp, left, ispin, etaLead):
        return self.sigs[path].getSig(hash, ee, kp, left, ispin, etaLead)

    def addSig(self, path, hash, ee, kp, left, ispin, etaLead, Sig):
        self.sigs[path].addSig(hash, ee, kp, left, ispin, etaLead, Sig)

    def close(self):
        for ii in self.sigs: self.sigs[ii].close()
        self.sigs={}

global SavedSig
SavedSig = SavedSigClass()


class ElectrodeSelfEnergy:
    """ 
    Calculate surface Greens function and self energy
    (should probably be renamed selfEnergy ...)
    For spinpolarized use the ispin given, for nonpolarized use 
    the same self-energy for both spin
    """
    global SavedSig

    def __init__(self, fn, NA1, NA2, voltage=0.0, UseF90helpers=True):
        self.path = os.path.split(os.path.abspath(fn))[0]
        self.HS=SIO.HS(fn, UseF90helpers=UseF90helpers) # An electrode HS
        self.hash=myHash([self.HS, NA1, NA2, voltage])
        SavedSig.add_hsfile(self.path)
        if self.HS.gamma:
            raise IOError("Are you trying to sneak a Gamma point electrode calculation past me?")
        self.NA1=NA1
        self.NA2=NA2
        self.kpoint = N.array([1e10, 1e10], N.float)
        self.voltage = voltage
        self.scaling = 1.0 # Default scale factor for coupling to device

    def getSig(self, ee, qp=N.array([0, 0], N.float), left=True, Bulk=False, ispin=0, UseF90helpers=True, etaLead=0.0, useSigNCfiles=False):
        """
        Get self-energy for specified 2-D surface k-point 
        Copy out g0 (surface greens function for smaller electrode calculation) 
        onto NA1*NA2*nuo matrix with the idiotic (TS) orbital order 
        a1(0,0) a1(1,0) .. a1(0,1) a1(1,1) ...... a2(0,0)
        Where a1, a2 ... are the atoms in the electrode calculation
        and (0,0) (1,0) indicate the replicating position.

        This gives self-energy from the solution of:
        Gs = (E S - H -Sigma)^-1 for Sigma

        For Bulk = True: Return E S - H - Sig (which we substitute into the bigger H)
        The voltage is assumed to be applied symmetrically and just shifts the energies of the self-energy
        """

        eeshifted = ee-self.voltage # Shift of self energies due to voltage

        if useSigNCfiles:
            Found, Sig = SavedSig.getSig(self.path, self.hash, eeshifted, qp, left*1, ispin, etaLead)
            if Found:
                if self.scaling!=1.0:
                    print 'NEGF.getSig: Scaling self-energy with a factor', self.scaling
                return Sig*self.scaling

        if ispin>=self.HS.nspin:
            ispin=0
            print "Warning: Non-spinpolarized electrode calculation used for both spin up and down"

        NA1, NA2 = self.NA1, self.NA2
        nuo = self.HS.nuo

        # First make it easy
        # If NA1*NA2 == 1 then no repetition is needed
        if NA1 * NA2 == 1:
            SGF = self.getg0(eeshifted+1j*etaLead, qp,
                             left=left, ispin=ispin, UseF90=UseF90helpers)
            if not Bulk:
                # We only need to calculate this quantity
                # if not bulk
                ESH = eeshifted* self.S-self.H[ispin, :, :]

        elif SIO.F90imported and UseF90helpers:

            g0   = N.empty((nuo, nuo, NA1*NA2), N.complex, order='F')
            mESH = N.empty((nuo, nuo, NA1*NA2), N.complex, order='F')

            iq = -1
            for ik2 in range(NA2):
                for ik1 in range(NA1):
                    iq += 1
                    kpoint = qp.copy() # Checked against 1x1 and 3x3 electrode calculation
                    kpoint[0] = (kpoint[0]+float(ik1))/float(NA1)
                    kpoint[1] = (kpoint[1]+float(ik2))/float(NA2)
                    # Surface GF with possible extra imaginary part (etaLead):
                    g0[:, :, iq] = self.getg0(eeshifted+1j*etaLead, kpoint, left=left, ispin=ispin, UseF90=UseF90helpers)
                    mESH[:, :, iq] = eeshifted*self.S-self.H[ispin, :, :]

            ESH, SGF = SIO.F90.expansion_se(no_u=nuo, no_s=nuo*NA1*NA2, \
                                                na1=NA1, na2=NA2,\
                                                kpt=qp, \
                                                na_u=self.HS.nua, \
                                                lasto=self.HS.lasto, \
                                                esh=mESH, g0=g0)
            # Clean up before calculating the self-energy
            del g0, mESH

        else:

            # Make list for the loop containing
            #       [iatom,i1,i2,SGFstart,SGFend,g0start,g0end]
            # where iatom is the atom number in g0 corresponding to g0start and g0end orbital
            # i1 and i2 posision in the repeated lattice which together with iatom gives SGFstart/end
            nua, lasto = self.HS.nua, self.HS.lasto
            SGFstart, loop = 0, []
            for ia in range(nua):         # Atoms in electrode
                g0start, g0end = lasto[ia], lasto[ia+1]-1
                for i2 in range(NA2):     # Repetition NA2
                    for i1 in range(NA1): # Repetition NA1
                        SGFend   = SGFstart+(g0end-g0start+1)-1 # add the number of orbitals in atom ia
                        loop.append([ia, i1, i2, SGFstart, SGFend, g0start, g0end])
                        SGFstart = SGFend+1
            if SGFstart!=NA1*NA2*nuo:
                raise ValueError("Error: Check of orbitals in making Sigma not correct")
            # Complete the full Gs with atoms copied out NA1*NA2
            # To obtain Sigma we also need H expanded, i.e.,
            # Gs = (E S - H - Sig)^-1 -> Sig = E S - H-SGF^-1
            # ESmH = E S - H
            SGF = N.zeros((NA1*NA2*nuo, NA1*NA2*nuo), N.complex)
            ESH = N.zeros((NA1*NA2*nuo, NA1*NA2*nuo), N.complex) # Temporary E S00 - H00
            for ik2 in range(NA2):
                for ik1 in range(NA1):
                    kpoint=qp.copy() # Checked against 1x1 and 3x3 electrode calculation
                    kpoint[0]=kpoint[0]/NA1
                    kpoint[1]=kpoint[1]/NA2
                    kpoint[0]+=ik1*1.0/NA1
                    kpoint[1]+=ik2*1.0/NA2
                    # Surface GF with possible extra imaginary part (etaLead):
                    g0 = self.getg0(eeshifted+1j*etaLead, kpoint, left=left, ispin=ispin, UseF90=UseF90helpers)

                    matESmH = eeshifted*self.S-self.H[ispin, :, :]
                    for ia, i1, i2, iSGFs, iSGFe, ig0s, ig0e in loop:
                        for ja, j1, j2, jSGFs, jSGFe, jg0s, jg0e in loop:
                            # Same convention as for H_ij above: exp(2 pi i k * (jatom-iatom))
                            # The phases etc have been checked by comparing the self-energy from
                            # 1x1 and 3x3 electrode calculations
                            phase = 1.0/(NA1*NA2)*N.exp(-2.0j*N.pi*((j1-i1)*kpoint[0]+(j2-i2)*kpoint[1]))
                            ESH[iSGFs:iSGFe+1, jSGFs:jSGFe+1] += phase*matESmH[ig0s:ig0e+1, jg0s:jg0e+1]
                            SGF[iSGFs:iSGFe+1, jSGFs:jSGFe+1] += phase*g0[ig0s:ig0e+1, jg0s:jg0e+1]

        # Calculate self-energy or inverse of SGF for Bulk: SGF^-1 = E S - H - Sig
        if Bulk:
            Sig = LA.inv(SGF) # SGF^1
        else:
            Sig = ESH - LA.inv(SGF)

        if useSigNCfiles:
            SavedSig.addSig(self.path, self.hash, eeshifted, qp, left, ispin, etaLead, Sig)
        if self.scaling!=1.0:
            print 'NEGF.getSig: Scaling self-energy with a factor', self.scaling
        return Sig*self.scaling

    def getg0(self, ee, kpoint, left=True, ispin=0, UseF90=True):
        # Calculate surface Green's function for small electrode calculation
        self.setupHS(kpoint)
        #print "NEGF.getg0: Constructing surface GF at (ReE,ImE) = (%.6e,%6e)"%(ee.real,ee.imag)

        if F90_lapack_imp and UseF90:
            return self.F90calcg0(ee, left=left, ispin=ispin)

        return self.calcg0_old(ee, left=left, ispin=ispin)
        #Potentially faster method but seems to have numerical instability
        #if hasSciPy and :
        #    return self.calcg0(ee,left=left,ispin=ispin)
        #else:
        #    return self.calcg0_old(ee,left=left,ispin=ispin)

    def calcg0(self, ee, ispin=0, left=True):
        # Calculate surface Green's function
        # Euro Phys J B 62, 381 (2008)
        # Inverse of : NOTE, setup for "right" lead.
        # e-h00 -h01  ...
        # -h10  e-h00 ...
        h00, s00, h01, s01 = self.H[ispin, :, :], self.S, self.H01[ispin, :, :], self.S01
        NN, ee = len(h00), N.real(ee)+N.max([N.imag(ee), 1e-8])*1.0j
        if left:
            h01, s01 = MM.dagger(h01), MM.dagger(s01)

        # Solve generalized eigen-problem
        # ( e I - h00 , -I) (eps)          (h01 , 0) (eps)
        # ( h10       ,  0) (xi ) = lambda (0   , I) (xi )
        a, b = N.zeros((2*NN, 2*NN), N.complex), N.zeros((2*NN, 2*NN), N.complex)
        a[0:NN, 0:NN] = ee*s00-h00
        a[0:NN, NN:2*NN] = -N.eye(NN)
        a[NN:2*NN, 0:NN] = MM.dagger(h01)-ee*MM.dagger(s01)
        b[0:NN, 0:NN] = h01-ee*s01
        b[NN:2*NN, NN:2*NN] = N.eye(NN)
        ev, evec = SLA.eig(a, b)

        # Select lambda <0 and the eps part of the evec
        ipiv = N.where(N.abs(ev)<1.0)[0]
        ev, evec = ev[ipiv], N.transpose(evec[:NN, ipiv])
        # Normalize evec
        norm = N.sqrt(N.diag(MM.mm(evec, MM.dagger(evec))))
        evec = MM.mm(N.diag(1.0/norm), evec)

        # E^+ Lambda_+ (E^+)^-1 --->>> g00
        EP = N.transpose(evec)
        FP = MM.mm(EP, N.diag(ev), LA.inv(MM.mm(MM.dagger(EP), EP)), MM.dagger(EP))
        g00 = LA.inv(ee*s00-h00-MM.mm(h01-ee*s01, FP))

        # Check!
        err=N.max(N.abs(g00-LA.inv(ee*s00-h00-\
                         MM.mm(h01-ee*s01, g00, MM.dagger(h01)-ee*MM.dagger(s01)))))
        if err>1.0e-8 and left:
            print "WARNING: Lopez-scheme not-so-well converged for LEFT electrode at E = %.4f eV:"%ee, err
        if err>1.0e-8 and not left:
            print "WARNING: Lopez-scheme not-so-well converged for RIGHT electrode at E = %.4f eV:"%ee, err
        return g00

    def F90calcg0(self, ee, ispin=0, left=True):
        """
        Call the fortran equivalent routine of the Lopez-Sancho algorithm also 
        utilised in tbtrans.
        Coded by Nick Papior Andersen
        """
        no = len(self.H[ispin, :, :])
        oHs = self.H.shape
        oSs = self.S.shape
        # due to an implementation not fully complete we need to transfer the shapes (luckily
        # the shape is just an internal parameter, and any data transfer is not needed...)
        self.H.shape = (len(self.H), -1)
        self.S.shape = (self.S.size,)
        self.H01.shape = self.H.shape
        self.S01.shape = self.S.shape
        tmp = F90.surfacegreen(no=no, ze=ee, h00=self.H[ispin], s00=self.S,
                               h01=self.H01[ispin], s01=self.S01,
                               accur=1.e-15, is_left=left)
        tmp.shape = (no, no)
        tmp = N.require(tmp, requirements=['A', 'C'])
        self.H.shape = oHs
        self.S.shape = oSs
        self.H01.shape = oHs
        self.S01.shape = oSs
        return tmp

    def calcg0_old(self, ee, ispin=0, left=True):
        """
        Only used if SciPy is not installed!
        For the left surface Green's function  (1 is surface layer, 0 is all the other atoms):
        (E S00-H00  E S01-H01)   (g00 g01)    ( I 0 )
        (E S10-H10  E S11-H11) * (g01 g11)  = ( 0 I ) ->
        call E S - H for t ...

        t00 g01 + t01 g11 = 0  -> g01 = - t00^-1 t01 g11
        t10 g01 + t11 g11 = I -> - t10 t00^-1 t01 g11 + t11 g11 = I -> 

        And we get the surface Green's function:

        g11 = (t11 - t10 t00^-1 t01)^-1 with the right size of unitcell t00^-1 = g11!
        g11 = (E S11 - H11 - (E S10 - H10) g11 (E S01 - H01))^-1

        In the calculations H01^+ and S01^+ are used instead of S10 and H10.
        (For complex energies (E S01 -H01)^+ is not (E S10 -H10) because the conjugate of the energy!!!!)

        For the right surface greens function same but different order on the MM.daggers!
        i.e., (E S - H - (E S01 - H01) gs (E S01^+ -H01^+)

        Algorith: Lopez Sancho*2 J Phys F:Met Phys 15 (1985) 851

        I'm still very suspicios of this algorithm ... but it works and is really quick! 
        The convergence is always checked against gs (E S - H - (E S01^+ - H01^+) gs (E S01 -H01) ) = I!
        """
        H, S, H01, S01 = self.H[ispin, :, :], self.S, self.H01[ispin, :, :], self.S01

        alpha, beta = MM.dagger(H01)-ee*MM.dagger(S01), H01-ee*S01
        eps, epss = H.copy(), H.copy()

        converged=False
        iteration=0
        while not converged:
            iteration+=1
            oldeps, oldepss = eps.copy(), epss.copy()
            oldalpha, oldbeta = alpha.copy(), beta.copy()
            tmpa=LA.solve(ee*S - oldeps, oldalpha)
            tmpb=LA.solve(ee*S - oldeps, oldbeta)
            alpha, beta = MM.mm(oldalpha, tmpa), MM.mm(oldbeta, tmpb)
            eps = oldeps + MM.mm(oldalpha, tmpb)+MM.mm(oldbeta, tmpa)
            if left:
                epss = oldepss + MM.mm(oldalpha, tmpb)
            else:
                epss = oldepss + MM.mm(oldbeta, tmpa)
            LopezConvTest=N.max(abs(alpha)+abs(beta))
            if LopezConvTest<1.0e-40:
                gs=LA.inv(ee*S-epss)
                if left:
                    test=ee*S-H-MM.mm(ee*MM.dagger(S01)-MM.dagger(H01), gs, ee*S01-H01)
                else:
                    test=ee*S-H-MM.mm(ee*S01-H01, gs, ee*MM.dagger(S01)-MM.dagger(H01))
                myConvTest=N.max(abs(MM.mm(test, gs)-N.identity((self.HS.nuo), N.complex)))
                if myConvTest<VC.GetCheck("Lopez-Sancho"):
                    converged=True
                    if myConvTest > VC.GetCheck("Lopez-Sancho-warning"):
                        v = "RIGHT"
                        if left: v = "LEFT"
                        print "WARNING: Lopez-scheme not-so-well converged for "+v+" electrode at E = %.4f eV:"%ee, myConvTest
                else:
                    VC.Check("Lopez-Sancho", myConvTest,
                             "Error: gs iteration {0}".format(iteration))
        return gs

    def setupHS(self, kpoint):
        """
        Setup H, S, H01 and S01 where H01 has large elements in the lower left corner, i.e., H01 = Hi,i+1
        (... 0     H01^+ H     H01   0      ...  )
        (... 0     0     H01^+ H     H01    0      ...  )
        """
        # Save time by not repeating too often
        if N.max(abs(kpoint-self.kpoint))>1e-10:
            self.kpoint = kpoint.copy()
            # Do the trick:
            # H(k=0)+H(kz=0.5) = H + H01 + H10 + H - H01 - H10 = 2 H
            kp = N.zeros((3), N.float)

            kp[0:2] = kpoint
            self.HS.setkpoint(kp, verbose=False)
            tmpH, tmpS = self.HS.H.copy(), self.HS.S.copy()

            kp[self.semiinf] = 0.5
            self.HS.setkpoint(kp, verbose=False)
            self.H = 0.5 * (tmpH + self.HS.H)
            self.S = 0.5 * (tmpS + self.HS.S)

            # Additional trick:
            # 1: -i*(H(kz=0.25)-H) = -i*(H + i*H01 - i*H10-H) = H01-H10
            # 2: H(kz=0)-H  = H + H01 + H10 - H =  H01+H10
            # -> H10 = (-i*(H(kz=0.25)-H) + H(kz=0)-H)/2
            kp[self.semiinf] = 0.25
            self.HS.setkpoint(kp, verbose=False)
            self.H01 = 0.5*(-1j*(self.HS.H - self.H) + tmpH - self.H)
            self.S01 = 0.5*(-1j*(self.HS.S - self.S) + tmpS - self.S)

            self.HS.resetkpoint()
            del tmpH, tmpS

#############################################################################


class GF:
    def __init__(self, TSHSfile, elecL, elecR, Bulk=True, DeviceAtoms=[0, 0], BufferAtoms=N.empty((0,))):
        """
        Calculate Green's functions etc for TSHSfile connected to left/right 
        electrode (class ElectrodeSelfEnergy). 
        To speed up calculations folding to smaller device region suggested
        For spin-polarized calcGF has to be called for each ispin
        Variables:
        Gr     : Retarded Green's function, folded
        H,S    : Hamiltonian, overlap, folded
        H0,S0  : Hamiltonian, overlap, not folded
        nuo    : Size of Gr
        SigL, SigR, GamL, GamR : Self energy, Gamma NOTE, not same size as Gr! 
        nuoL, nuoR : Size of Sig, Gam
        nuo0, nuoL0, nuoR0 : Non-folded sizes
        FoldedL, FoldedR : True/False
        DeviceAtoms : start/end Siesta numbering of atoms included in device
        DeviceOrbs : Start/end of orbitals. Siesta ordering.
        BufferAtoms: A list of buffer atoms
        """
        self.elecL, self.elecR, self.Bulk = elecL, elecR, Bulk
        self.HS = SIO.HS(TSHSfile, BufferAtoms=BufferAtoms)
        print 'GF: UseBulk=', Bulk
        self.DeviceAtoms=DeviceAtoms
        if DeviceAtoms[0]<=1:
            self.DeviceAtoms[0]=1
            self.FoldedL = False
        else:
            self.FoldedL = True
        if DeviceAtoms[1]==0 or DeviceAtoms[1]>=self.HS.nua:
            self.DeviceAtoms[1]=self.HS.nua
            self.FoldedR = False
        else:
            self.FoldedR = True
        self.DeviceOrbs = [self.HS.lasto[DeviceAtoms[0]-1]+1, self.HS.lasto[DeviceAtoms[1]]]

        self.nuo0, self.nuoL0, self.nuoR0 = self.HS.nuo, elecL.NA1*elecL.NA2*elecL.HS.nuo, elecR.NA1*elecR.NA2*elecR.HS.nuo
        self.nuo = self.DeviceOrbs[1]-self.DeviceOrbs[0]+1
        self.nuoL, self.nuoR = self.nuoL0, self.nuoR0 # Not folded, for folded case changed below

        print "GF:", TSHSfile
        print "Device atoms %i-%i, orbitals %i-%i"%(tuple(self.DeviceAtoms+self.DeviceOrbs))
        if not self.FoldedL:
            print "Suggest left folding to atom : ", self.elecL.HS.nua*self.elecL.NA1*self.elecL.NA2+1
        if not self.FoldedR:
            print "Suggest right folding to atom : ", self.HS.nua-self.elecR.HS.nua*self.elecR.NA1*self.elecR.NA2

        if self.FoldedL or self.FoldedR:
            # Check that device region is large enough!
            kpoint=N.zeros((2,), N.float)
            self.setkpoint(kpoint, ispin=0) # At least for one spin

            devSt, devEnd = self.DeviceOrbs[0], self.DeviceOrbs[1]
            VC.Check("Device-Elec-overlap", N.abs(self.S0[0:devSt, devEnd:self.nuo0]),
                     "Too much overlap directly from left-top right",
                     "Make device region larger")
            VC.Check("Device-Elec-overlap", N.abs(self.H0[0:devSt, devEnd:self.nuo0]),
                     "Too large Hamiltonian directly from left-top right.",
                     "Make device region larger")
        if self.FoldedL:
            # Find orbitals in device region coupling to left and right.
            tau  = abs(self.S0[0:devSt-1, 0:devEnd])
            coupling = N.sum(tau, axis=0)
            ii=devEnd-1
            while coupling[ii]<1e-10: ii=ii-1
            self.devEndL = max(ii+1, self.nuoL0)
            self.nuoL = self.devEndL-devSt+1
            print "Left self energy on orbitals %i-%i"%(devSt, self.devEndL)
        if self.FoldedR:
            tau  = abs(self.S0[devEnd-1:self.nuo0, 0:self.nuo0])
            coupling = N.sum(tau, axis=0)
            ii=devSt-1
            while coupling[ii]<1e-10: ii=ii+1
            self.devStR = min(ii+1, self.nuo0-self.nuoR0+1)
            self.nuoR = devEnd-self.devStR+1
            print "Right self energy on orbitals %i-%i"%(self.devStR, devEnd)
        # Quantities expressed in nonorthogonal basis:
        self.OrthogonalDeviceRegion = False

    def calcSigLR(self, ee, kpoint, ispin=0, etaLead=0.0, useSigNCfiles=False, SpectralCutoff=0.0):
        """
        Calculate (folded) self-energy at energy ee and 2d k-point
        Uses SpectralMatrix format for the spectralfunction matrices, see MiscMath, if cutoff>0.0
        """

        nuo, nuoL, nuoR = self.nuo, self.nuoL, self.nuoR
        nuo0, nuoL0, nuoR0 = self.nuo0, self.nuoL0, self.nuoR0
        FoldedL, FoldedR = self.FoldedL, self.FoldedR
        devSt, devEnd = self.DeviceOrbs[0], self.DeviceOrbs[1]
        # Calculate Sigma without folding
        self.setkpoint(kpoint, ispin)
        SigL0 = self.elecL.getSig(ee, kpoint, left=True, Bulk=self.Bulk, ispin=ispin, etaLead=etaLead, useSigNCfiles=useSigNCfiles)
        SigR0 = self.elecR.getSig(ee, kpoint, left=False, Bulk=self.Bulk, ispin=ispin, etaLead=etaLead, useSigNCfiles=useSigNCfiles)

        if FoldedL:
            # Fold down from nuoL0 to the device region
            # A11 A12     g11 g12    I 0
            # A21 A22  *  g21 g22  = 0 I ->
            # g22 = (A22-A21.A11^-1.A12)^-1 ->
            # Sigma = A21.A11^-1.A12          (tau=A12)
            devEndL = self.devEndL
            # Do folding
            eSmH = ee*self.S0-self.H0
            eSmHmS = eSmH[0:devEndL, 0:devEndL].copy()
            if self.Bulk:
                eSmHmS[0:nuoL0, 0:nuoL0] = SigL0 # SGF^1
            else:
                eSmHmS[0:nuoL0, 0:nuoL0] = eSmHmS[0:nuoL0, 0:nuoL0]-SigL0
            tau  = eSmHmS[0:devSt-1, devSt-1:devEndL].copy()
            taud = eSmHmS[devSt-1:devEndL, 0:devSt-1].copy()
            inv = LA.inv(eSmHmS[0:devSt-1, 0:devSt-1])
            eSmHmS[devSt-1:devEndL, devSt-1:devEndL]=eSmHmS[devSt-1:devEndL, devSt-1:devEndL]-\
                MM.mm(taud, inv, tau)
            self.SigL = eSmH[devSt-1:devEndL, devSt-1:devEndL]-eSmHmS[devSt-1:devEndL, devSt-1:devEndL]
        else:
            self.SigL = SigL0
        self.GamL = 1.0j*(self.SigL-MM.dagger(self.SigL))
        if self.Bulk and not FoldedL:
            # Reverse sign since SigL is really SGF^-1
            self.GamL = -1.0*self.GamL
        AssertReal(N.diag(self.GamL), 'GamL')

        if FoldedR:
            # Fold down from nuoR0 to the device region
            devStR = self.devStR
            eSmH = ee*self.S0-self.H0
            eSmHmS = eSmH[devStR-1:nuo0, devStR-1:nuo0].copy()
            tmpnuo=len(eSmHmS)
            if self.Bulk:
                eSmHmS[tmpnuo-nuoR0:tmpnuo, tmpnuo-nuoR0:tmpnuo] = SigR0 # SGF^1
            else:
                eSmHmS[tmpnuo-nuoR0:tmpnuo, tmpnuo-nuoR0:tmpnuo] = eSmHmS[tmpnuo-nuoR0:tmpnuo, tmpnuo-nuoR0:tmpnuo]-SigR0
            tau  = eSmHmS[0:nuoR, nuoR:tmpnuo].copy()
            taud = eSmHmS[nuoR:tmpnuo, 0:nuoR].copy()
            inv = LA.inv(eSmHmS[nuoR:tmpnuo, nuoR:tmpnuo])
            eSmHmS[0:nuoR, 0:nuoR]=eSmHmS[0:nuoR, 0:nuoR]-MM.mm(tau, inv, taud)
            self.SigR = eSmH[devStR-1:devEnd, devStR-1:devEnd]-eSmHmS[0:nuoR, 0:nuoR]
        else:
            self.SigR = SigR0
        self.GamR = 1.0j*(self.SigR-MM.dagger(self.SigR))
        if self.Bulk and not FoldedR:
            # Reverse sign since SigR is really SGF^-1
            self.GamR = -1.0*self.GamR
        AssertReal(N.diag(self.GamR), 'GamR')

    def calcGF(self, ee, kpoint, ispin=0, etaLead=0.0, useSigNCfiles=False, SpectralCutoff=0.0):
        "Calculate GF etc at energy ee and 2d k-point"
        nuo, nuoL, nuoR = self.nuo, self.nuoL, self.nuoR
        nuo0, nuoL0, nuoR0 = self.nuo0, self.nuoL0, self.nuoR0
        FoldedL, FoldedR = self.FoldedL, self.FoldedR
        devSt, devEnd = self.DeviceOrbs[0], self.DeviceOrbs[1]

        # Determine whether electrode self-energies should be k-sampled or not
        try:    mesh = self.elecL.mesh # a mesh was attached
        except: mesh = False
        # Calculate electrode self-energies
        if mesh:
            try:    self.SigAvg # Averaged self-energies exist
            except: self.SigAvg = [False, -1]
            if self.SigAvg[0] == ee and self.SigAvg[1] == ispin:
                # We have already the averaged self-energies
                print 'NEGF: Reusing sampled electrode self-energies', mesh.Nk, mesh.type, 'for ispin= %i e= %f'%(ispin, ee)
            else:
                # k-sampling performed over folded electrode self-energies
                print 'NEGF: Sampling electrode self-energies', mesh.Nk, mesh.type, 'for ispin= %i e= %f'%(ispin, ee)
                self.calcSigLR(ee, mesh.k[0, :2], ispin, etaLead, useSigNCfiles, SpectralCutoff)
                AvgSigL = mesh.w[0, 0]*self.SigL
                AvgSigR = mesh.w[0, 0]*self.SigR
                for i in range(1, len(mesh.k)):
                    self.calcSigLR(ee, mesh.k[i, :2], ispin, etaLead, useSigNCfiles, SpectralCutoff)
                    AvgSigL += mesh.w[0, i]*self.SigL
                    AvgSigR += mesh.w[0, i]*self.SigR
                # We now simply continue with the averaged self-energies
                self.SigL = AvgSigL
                self.SigR = AvgSigR
                self.SigAvg = [ee, ispin]
        else:
            # We sample k-points the usual way
            self.calcSigLR(ee, kpoint, ispin, etaLead, useSigNCfiles)

        # Ready to calculate Gr
        self.setkpoint(kpoint, ispin)
        eSmH=ee*self.S-self.H
        if FoldedL:
            eSmH[0:nuoL, 0:nuoL]=eSmH[0:nuoL, 0:nuoL]-self.SigL
        else:
            if self.Bulk:
                eSmH[0:nuoL, 0:nuoL] = self.SigL # SGF^1
            else:
                eSmH[0:nuoL, 0:nuoL] = eSmH[0:nuoL, 0:nuoL]-self.SigL
        if FoldedR:
            eSmH[nuo-nuoR:nuo, nuo-nuoR:nuo]=eSmH[nuo-nuoR:nuo, nuo-nuoR:nuo]-self.SigR
        else:
            if self.Bulk:
                eSmH[nuo-nuoR:nuo, nuo-nuoR:nuo] = self.SigR # SGF^1
            else:
                eSmH[nuo-nuoR:nuo, nuo-nuoR:nuo] = eSmH[nuo-nuoR:nuo, nuo-nuoR:nuo]-self.SigR
        self.Gr = LA.inv(eSmH)
        self.Ga = MM.dagger(self.Gr)
        # Calculate spectral functions
        if SpectralCutoff>0.0:
            self.AL = MM.SpectralMatrix(MM.mm(self.Gr[:, 0:nuoL], self.GamL, self.Ga[0:nuoL, :]), cutoff=SpectralCutoff)
            tmp = MM.mm(self.GamL, self.Gr[0:nuoL, :])
            self.ALT = MM.SpectralMatrix(MM.mm(self.Ga[:, 0:nuoL], tmp), cutoff=SpectralCutoff)
            self.AR = MM.SpectralMatrix(MM.mm(self.Gr[:, nuo-nuoR:nuo], self.GamR, self.Ga[nuo-nuoR:nuo, :]), cutoff=SpectralCutoff)
            self.ARGLG = MM.mm(self.AR.L, self.AR.R[:, 0:nuoL], tmp)
            self.A = self.AL+self.AR
            # transmission matrix AL.GamR
            self.TT = MM.mm(self.AL.R[:, nuo-nuoR:nuo], self.GamR, self.AL.L[nuo-nuoR:nuo, :])
        else:
            self.AL = MM.mm(self.Gr[:, 0:nuoL], self.GamL, self.Ga[0:nuoL, :])
            tmp = MM.mm(self.GamL, self.Gr[0:nuoL, :])
            self.ALT = MM.mm(self.Ga[:, 0:nuoL], tmp)
            self.AR = MM.mm(self.Gr[:, nuo-nuoR:nuo], self.GamR, self.Ga[nuo-nuoR:nuo, :])
            self.ARGLG = MM.mm(self.AR[:, 0:nuoL], tmp)
            self.A = self.AL+self.AR
            # transmission matrix AL.GamR
            self.TT = MM.mm(self.AL[nuo-nuoR:nuo, nuo-nuoR:nuo], self.GamR)

        print 'NEGF.calcGF: Shape of transmission matrix (TT):', self.TT.shape
        print 'NEGF.calcGF: Energy and total transmission Tr[TT].real:', ee, N.trace(self.TT).real
        # Write also the Gammas in the full space of Gr/Ga/A
        # (needed for the inelastic shot noise)
        self.GammaL = N.zeros(self.Gr.shape, N.complex)
        self.GammaL[0:nuoL, 0:nuoL] = self.GamL
        self.GammaR = N.zeros(self.Gr.shape, N.complex)
        self.GammaR[nuo-nuoR:nuo, nuo-nuoR:nuo] = self.GamR

    def setkpoint(self, kpoint, ispin=0):
        # Initiate H, S to correct kpoint
        nuo, nuoL, nuoR = self.nuo0, self.nuoL0, self.nuoR0

        kpoint3 = N.zeros((3), N.float)
        kpoint3[0:2] = kpoint[:]
        self.HS.setkpoint(kpoint3, verbose=False)
        # Remove PBC in z-direction
        if self.HS.gamma:
            self.H0 = self.HS.H[ispin, :, :].copy()
            self.S0 = self.HS.S.copy()
            # Remove direct left/right coupling
            self.H0[0:nuoL, nuo-nuoR:nuo] = 0.
            self.H0[nuo-nuoR:nuo, 0:nuoL] = 0.
            self.S0[0:nuoL, nuo-nuoR:nuo] = 0.
            self.S0[nuo-nuoR:nuo, 0:nuoL] = 0.
        else:
            # Do trick with kz
            tmpH, tmpS = self.HS.H[ispin, :, :].copy(), self.HS.S.copy()
            if self.elecL.semiinf==0 and self.elecR.semiinf==0:
                # Periodicity along A1
                if kpoint[0] == 0.0:
                    kpoint3[0] = 0.5
                else:
                    print 'Specified 2D k-point=', kpoint
                    raise IOError('Incompatible 2D k-point - use z as transport direction')
            elif self.elecL.semiinf==1 and self.elecR.semiinf==1:
                # Periodicity along A1
                if kpoint[1] == 0.0:
                    kpoint3[1] = 0.5
                else:
                    print 'Specified 2D k-point=', kpoint
                    raise IOError('Incompatible 2D k-point - use z as transport direction')
            else:
                # Default is along A3
                kpoint3[2] = 0.5
            self.HS.setkpoint(kpoint3, verbose=False)
            self.H0 = 0.5 * (tmpH + self.HS.H[ispin, :, :])
            self.S0 = 0.5 * (tmpS + self.HS.S)

        if self.FoldedL or self.FoldedR:
            devSt, devEnd = self.DeviceOrbs[0], self.DeviceOrbs[1]
            self.H = self.H0[devSt-1:devEnd, devSt-1:devEnd]
            self.S = self.S0[devSt-1:devEnd, devSt-1:devEnd]
        else:
            self.H, self.S = self.H0, self.S0

        self.OrthogonalDeviceRegion = False

    def calcTEIG(self, channels=10):
        # Transmission matrix (complex array)
        TT = self.TT
        Trans = N.trace(TT)
        VC.Check("trans-imaginary-part", Trans.imag,
                 "Transmission has large imaginary part")
        # Calculate eigenchannel transmissions too
        tval, tvec = LA.eig(TT)
        idx = (tval.real).argsort()[::-1] # sort from largest to smallest
        tval = tval[idx]
        tvec = tvec[:, idx]
        # Compute shot noise
        Smat = MM.mm(TT, N.identity(len(TT))-TT)
        sval = N.diag(MM.mm(MM.dagger(tvec), Smat, tvec))
        # set up arrays
        T = N.zeros(channels+1)
        SN = N.zeros(channels+1)
        T[0] = Trans.real
        SN[0] = N.trace(Smat).real
        for i in range(min(channels, len(TT))):
            T[i+1] = tval[i].real
            SN[i+1] = sval[i].real
        return T, SN

    def orthogonalize(self):
        print 'NEGF.GF.orthogonalize: Orthogonalizing device region quantities'
        self.OrthogonalDeviceRegion = True
        self.HNO = self.H.copy() # nonorthogonal device Hamiltonian (needed)

        # Device part
        Usi = MM.mysqrt(self.S) # Folded S
        Us = LA.inv(Usi)
        # Store transformation matrices
        self.Usi, self.Us = Usi, Us

        # Transform S and H
        self.S, self.H = MM.mm(Us, self.S, Us), MM.mm(Us, self.H, Us)

        # Sigmas/Gammas in pyTBT GF can be smaller than device region
        # First give them the shape of the device region
        nnL, nnR = len(self.SigL),  len(self.SigR)
        S1, S2 = N.zeros(self.H.shape, N.complex), N.zeros(self.H.shape, N.complex)
        S1[0:nnL, 0:nnL],  S2[-nnR:, -nnR:] = self.SigL, self.SigR
        # Resetting Sigmas to orthogonalized quantities
        self.SigL, self.SigR = MM.mm(Us, S1, Us), MM.mm(Us, S2, Us)
        # ... now the same for the Gammas
        G1, G2 = N.zeros(self.H.shape, N.complex), N.zeros(self.H.shape, N.complex)
        G1[0:nnL, 0:nnL],  G2[-nnR:, -nnR:] = self.GamL, self.GamR
        # Resetting Gammas to orthogonalized quantities
        self.GamL, self.GamR = MM.mm(Us, G1, Us), MM.mm(Us, G2, Us)

        # Orthogonalize Greens functions
        self.Gr = MM.mm(Usi, self.Gr, Usi)
        self.Ga = MM.dagger(self.Gr)

    def __calcEigChan(self, A1, G2, Left, channels=10):
        # Calculate Eigenchannels using recipe from PRB
        # For right eigenchannels, A1=A2, G2=G1 !!!
        if isinstance(A1, MM.SpectralMatrix):
            ev, U = LA.eigh(MM.mm(A1.L, A1.R))
        else:
            ev, U = LA.eigh(A1)

        # This small trick will remove all zero contribution vectors
        # and will diagonalize the tt matrix in the subspace where there
        # are values.
        idx = (ev > 0).nonzero()[0]
        ev = N.sqrt(ev[idx] / (2 * N.pi))
        ev.shape = (1, -1)
        Utilde = ev * U[:, idx]

        nuo, nuoL, nuoR = self.nuo, self.nuoL, self.nuoR
        if Left:
            tt=MM.mm(MM.dagger(Utilde[nuo-nuoR:nuo, :]), 2*N.pi*G2, Utilde[nuo-nuoR:nuo, :])
        else:
            tt=MM.mm(MM.dagger(Utilde[:nuoL, :]), 2*N.pi*G2, Utilde[:nuoL, :])

        # Diagonalize (note that this is on a reduced tt matrix (no 0 contributing columns)
        evF, UF = LA.eigh(tt)
        EC = MM.mm(Utilde, UF[:, -channels:]).T
        return EC[::-1, :], evF[::-1] # reverse eigenvalues

    def calcEigChan(self, channels=10):
        # Calculate Eigenchannels from left
        self.ECleft, self.EigTleft = self.__calcEigChan(self.AL, self.GamR, True, channels)
        print 'NEGF.calcEigChan: Left eigenchannel transmissions [T1, ..., Tn]:\n', self.EigTleft[:channels]
        # Calculate Eigenchannels from right
        self.ECright, self.EigTright = self.__calcEigChan(self.AR, self.GamL, False, channels)
        print 'NEGF.calcEigChan: Right eigenchannel transmissions [T1, ..., Tn]:\n', self.EigTright[:channels]
