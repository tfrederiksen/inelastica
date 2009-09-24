#################################################################
#
# python TBTrans 
# Magnus Paulsson magnus.paulsson@hik.se
#
# Requires: numpy (compile it linked with mkl, acml or atlas!)
#           ScientificPython (vers. >= 2.8)
#           For speed compile the fortran subroutines in F90 
#           (cd F90;source compile.bat)
#
# Contains:
#  class HS moved to SiestaIO.py
#  ( class HS : Reads .TSHS.nc file and folds spcific k-point into 
#              normal matrix from sparse format. )
#   class surfaceGF : Reads electrode TSHS.nc and returns surface
#                     Green's function of supercell surface 
#   class GF : Calculates Green's function etc 
#
#  UNITS! Always eV and Angstrom!
#         k-values always given in range [0,1.0] (or [-0.5,0.5])
#         They are not in reciprocal space. Instead they corresponds
#         to the mathematical orthogonal space that is fourier 
#         transformed.
#
#################################################################



import SiestaIO as SIO
import numpy as N
import numpy.linalg as LA
import profile, sys, string
from optparse import OptionParser, OptionGroup

################### Help functions ############################
try:
    import F90helpers as F90
    F90imported = True
except:
    F90imported = False
    print "########################################################"
    print "Perhaps time to compile F90/setkpointhelper"
    print "Try:" 
    print "        cd F90;source compile.bat"
    print "########################################################"

def mm(* args):
    # Matrix multiplication with arbitrary number of arguments
    tmp=N.dot(args[0],args[1])
    for ii in range(len(args)-2):
        tmp=N.dot(tmp,args[ii+2])
    return tmp

def dag(x):
    return N.transpose(N.conjugate(x))


################### Main program ############################
def main(pyTBT=True,deviceRegion=[0,0],fn=None):
    """
    Running standalone to calculate transmission 
    *OR* 
      called from Eigenchannels or Inelastica 
      returning elecL, elecR, GF, deviceStart, deviceEnd
    """
    
    if pyTBT: 
        usage = "usage: %prog RUN.fdf"
        descr = "pyTBT is the Python version of TBtrans originally developed by Mads Brandbyge."
        intro = """
pyTBT is the Python version of TBtrans originally developed by Mads Brandbyge.

pyTBT reads some of the TBT and TS keywords from the fdf file:

Electrodes:
TS.HSFileLeft         filename.TSHS
TS.ReplicateA1Left    1
TS.ReplicateA2Left    1
TS.HSFileRight        filename.TSHS
TS.ReplicateA1Right   1
TS.ReplicateA2Right   1
Note: Fredericos TBtrans and the Transiesta version planned to be release in 2009 cannot use ReplicateA1,2 but pyTBT can.

Device region:
TS.TBT.PDOSFrom       10       [default=1]
TS.TBT.PDOSTo         20       [default=last atom]
Note: If you just want transmission pyTBT is quickest if the device region
      is the middle 1/3 of the orbitals.

Transmission energies [default]:
TS.TBT.NPoints        21              
TS.TBT.Emin          -1.000000 eV  
TS.TBT.Emax           1.000000 eV 

How self-energies are applied:
TS.UseBulkInElectrodes .True.
Note, False for this option does not seem to be a good choice.

NEW KEYWORDS:
pyTBT.eta             0.000001 eV [default, imaginary part of energy]

Kpoint sampling of transmission:
pyTBT.K_A1            1           [default=1]
pyTBT.K_A2            1


Ouputfiles:
SystemLabel[.UP/.DOWN].TRANS     Transmission k-point dependent.
SystemLabel[.UP/.DOWN].AVTRANS   Averaged over k-points.

"""
        parser = OptionParser(usage,description=descr)
        print intro
        parser.parse_args()

    # Read options
    ##############################################################################
    if fn==None:
        try: 
            fn = sys.argv[1]
            print "pyTBT reading keywords from ",fn
        except:
            fn = 'RUN.fdf'
            print "pyTBT::WARNING reading keywords from default file : ",fn

    # Electrodes
    fnL  =SIO.GetFDFlineWithDefault(fn,'TS.HSFileLeft', str, None, 'pyTBT')
    NA1L =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA1Left', int, 1, 'pyTBT')
    NA2L =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA2Left', int, 1, 'pyTBT')
    fnR  =SIO.GetFDFlineWithDefault(fn,'TS.HSFileRight', str, None, 'pyTBT')
    NA1R =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA1Right', int, 1, 'pyTBT')
    NA2R =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA2Right', int, 1, 'pyTBT')

    # Device region
    if deviceRegion[0]==0:
        devSt =SIO.GetFDFlineWithDefault(fn,'TS.TBT.PDOSFrom', int, 0, 'pyTBT')
    else:
        devSt=deviceRegion[0]
    if deviceRegion[1]==0:
        devEnd=SIO.GetFDFlineWithDefault(fn,'TS.TBT.PDOSTo', int, 0, 'pyTBT')
    else:
        devEnd=deviceRegion[1]

    # Energy range
    nE  =SIO.GetFDFlineWithDefault(fn,'TS.TBT.NPoints', int, 21, 'pyTBT')
    minE=SIO.GetFDFlineWithDefault(fn,'TS.TBT.Emin', float, -1.0, 'pyTBT')
    maxE=SIO.GetFDFlineWithDefault(fn,'TS.TBT.Emax', float, 1.0, 'pyTBT')
    if nE>1:
        dE = (maxE-minE)/float(nE-1)
        Elist = N.array(range(int((maxE-minE+1e-9)/dE)+1),N.float)*dE+minE
    else:
        dE=0.0
        Elist=N.array((minE,),N.float)

    UseBulk=SIO.GetFDFlineWithDefault(fn,'TS.UseBulkInElectrodes', bool, True, 'pyTBT')

    eta=SIO.GetFDFlineWithDefault(fn,'pyTBT.eta', float, 0.000001, 'pyTBT')
    Nk1=SIO.GetFDFlineWithDefault(fn,'pyTBT.K_A1', int, 1, 'pyTBT')
    Nk2=SIO.GetFDFlineWithDefault(fn,'pyTBT.K_A2', int, 1, 'pyTBT')

    outFile=SIO.GetFDFlineWithDefault(fn,'SystemLabel', str, None, 'pyTBT')

    #=SIO.GetFDFlineWithDefault(fn,'', int, , 'pyTBT')

    ##############################################################################
    # Define electrodes and device

    elecL=surfaceGF(fnL,NA1L,NA2L)
    elecR=surfaceGF(fnR,NA2R,NA2R)
    myGF = GF(outFile+'.TSHS',elecL,elecR,Bulk=UseBulk,DeviceAtoms=[devSt, devEnd])
    nspin = myGF.HS.nspin
    if devSt==0:
        devSt=GF.DeviceAtoms[0]
    if devEnd==0:
        devEnd=GF.DeviceAtoms[1]
        
    print """
##############################################################
pyTBT

Energy [eV]                     : %f:%f:%f
kpoints                         : %i, %i 
eta [eV]                        : %f 
Device [Atoms Siesta numbering] : %i:%i 
Bulk                            : %s
SpinPolarization                : %i
##############################################################

"""%(minE,dE,maxE,Nk1,Nk2,eta,devSt,devEnd,UseBulk,nspin)

    if not pyTBT:
        # For Eigenchannels and inelastica!
        return elecL, elecR, myGF, devSt, devEnd, Elist, eta, outFile
    else:
        Tkpt=N.zeros((len(Elist),Nk1,Nk2),N.float)
        for iSpin in range(nspin):
            if nspin<2:
                fo=open(outFile+'.AVTRANS','write')
            else:
                fo=open(outFile+['.UP','.DOWN'][iSpin]+'.AVTRANS','write')
            for ie, ee in enumerate(Elist):
                Tavg = 0.0
                for ik1 in range(Nk1):
                    for ik2 in range(Nk2):
                        kpt=N.array([ik1/float(Nk1),ik2/float(Nk2)],N.float)
                        myGF.calcGF(ee+eta*1.0j,kpt,ispin=iSpin)
                        T = myGF.calcT()
                        Tavg = Tavg + T/Nk1/Nk2
                        Tkpt[ie,ik1,ik2] = T
                print ee," ",Tavg
                fo.write('%f %f \n'%(ee,Tavg))
            fo.close()
        
            # Write k-point transmission
            if nspin<2:
                fo=open(outFile+'.TRANS','write')
            else:
                fo=open(outFile+['.UP','.DOWN'][iSpin]+'.TRANS','write')
            for ik1 in range(Nk1):
                for ik2 in range(Nk2):
                    fo.write('\n# k = %f, %f \n'%(ik1/float(Nk1),ik2/float(Nk2)))
                    for ie, ee in enumerate(Elist):
                        fo.write('%f %f \n'%(ee,Tkpt[ie,ik1,ik2]))
            fo.close()


class surfaceGF:
    """ 
    Calculate surface Greensfunction and self energy
    (should probably be renamed selfEnergy ...)
    For spinpolarized use the ispin given, for nonpolarized use 
    the same self-energy for both spin 
    """
    def __init__(self,fn,NA1,NA2,UseF90helpers=True):
        self.HS=SIO.HS(fn,UseF90helpers=UseF90helpers)
        if self.HS.gamma:
            print "Are you trying to sneak a Gamma point electrode calculation past me?"
            kuk
        self.NA1=NA1
        self.NA2=NA2
        self.kpoint = N.array([1e10,1e10],N.float)

    def getSig(self,ee,qp=N.array([0,0],N.float),left=True,Bulk=False,ispin=0,UseF90helpers=True):
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
        """
        if ispin>self.HS.nspin:
            ispin=0
            print "Warning: Non-spinpolarized electrode calculation used for both spin up and down"

        # Make list for the loop containing 
        #       [iatom,i1,i2,SGFstart,SGFend,g0start,g0end] 
        # where iatom is the atom number in g0 corresponding to g0start and g0end orbital
        # i1 and i2 posision in the repeated lattice which together with iatom gives SGFstart/end
        NA1, NA2 = self.NA1, self.NA2
        nua, lasto, nuo = self.HS.nua, self.HS.lasto, self.HS.nuo 
        SGFstart, loop = 0, []
        for ia in range(nua):         # Atoms in electrode
            for i2 in range(NA2):     # Repetition NA2
                for i1 in range(NA1): # Repetition NA1
                    g0start, g0end = lasto[ia],lasto[ia+1]-1          
                    SGFend=SGFstart+(g0end-g0start+1)-1               # add the number of orbitals in atom ia
                    tmp=[ia,i1,i2,SGFstart,SGFend,g0start,g0end]
                    loop.append(tmp)
                    SGFstart=SGFend+1

        if SGFstart!=NA1*NA2*nuo:
            print "Error: Check of orbitals in making Sigma not correct"
            kuk

        # Complete the full Gs with atoms copied out NA1*NA2
        # To obtain Sigma we also need H expanded, i.e., 
        # Gs = (E S - H - Sig)^-1 -> Sig = E S - H-SGF^-1 
        # ESmH = E S - H
        SGF = N.zeros((NA1*NA2*nuo,NA1*NA2*nuo),N.complex)
        ESmH =N.zeros((NA1*NA2*nuo,NA1*NA2*nuo),N.complex) # Temporary E S00 - H00
        for ik1 in range(NA1):
            for ik2 in range(NA2):
                kpoint=qp.copy()                     # Checked against 1x1 and 3x3 electrode calculation
                kpoint[0]=kpoint[0]/NA1
                kpoint[1]=kpoint[1]/NA2
                kpoint[0]+=ik1*1.0/NA1
                kpoint[1]+=ik2*1.0/NA2
                g0=self.getg0(ee,kpoint,left=left,ispin=ispin)             
                matESmH = ee*self.S-self.H[ispin,:,:]
                if F90imported and UseF90helpers:
                    ESmH, SGF = F90.f90distributegs(loop=N.array(loop,N.int), nuo=nuo,\
                                  nua=nua, na1=NA1, na2=NA2, kpoint=kpoint,\
                                  matesmh=matESmH, g0=g0, esmh=ESmH, sgf=SGF)
                else:
                    for ia, i1, i2, iSGFs, iSGFe, ig0s, ig0e in loop:
                        for ja, j1, j2, jSGFs, jSGFe, jg0s, jg0e in loop:
                            # Same convention as for H_ij above: exp(2 pi i k * (jatom-iatom))
                            # The phases etc have been checked by comparing the self-energy from 
                            # 1x1 and 3x3 electrode calculations
                            phase = 1.0/(NA1*NA2)*N.exp(2.0j*N.pi*((j1-i1)*kpoint[0]+(j2-i2)*kpoint[1])) 
                            ESmH[iSGFs:iSGFe+1,jSGFs:jSGFe+1]=ESmH[iSGFs:iSGFe+1,jSGFs:jSGFe+1]+\
                                N.conjugate(phase)*matESmH[ig0s:ig0e+1,jg0s:jg0e+1]       
                            SGF[iSGFs:iSGFe+1,jSGFs:jSGFe+1]=SGF[iSGFs:iSGFe+1,jSGFs:jSGFe+1]+\
                                N.conjugate(phase)*g0[ig0s:ig0e+1,jg0s:jg0e+1]

        # Calculate self-energy or inverse of SGF for Bulk: SGF^-1 = E S - H - Sig
        if not Bulk:
            Sig=ESmH-LA.inv(SGF)
        else:
            Sig=LA.inv(SGF)
        return Sig

    def getg0(self,ee,kpoint,left=True,ispin=0):
        # Calculate surface Green's function for small electrode calculation
        self.setupHS(kpoint)
        return self.calcg0(ee,left=left,ispin=ispin)
        
    def calcg0(self,ee,ispin=0,left=True):
        """
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

        For the right surface greens function same but different order on the daggers!
        i.e., (E S - H - (E S01 - H01) gs (E S01^+ -H01^+)

        Algorith: Lopez Sancho*2 J Phys F:Met Phys 15 (1985) 851
        
        I'm still very suspicios of this algorithm ... but it works and is really quick! 
        The convergence is always checked against gs (E S - H - (E S01^+ - H01^+) gs (E S01 -H01) ) = I!
        """
        H, S, H01, S01 = self.H[ispin,:,:] ,self.S ,self.H01[ispin,:,:], self.S01

        alpha, beta = dag(H01)-ee*dag(S01), H01-ee*S01
        eps, epss = H.copy(), H.copy()
        
        converged=False
        iteration=0
        while not converged:
            iteration+=1
            oldeps, oldepss = eps.copy(), epss.copy()
            oldalpha, oldbeta = alpha.copy(), beta.copy()
            tmp=LA.inv(ee*S - oldeps)
            alpha, beta = mm(oldalpha,tmp,oldalpha), mm(oldbeta,tmp,oldbeta)
            eps = oldeps + mm(oldalpha,tmp,oldbeta)+mm(oldbeta,tmp,oldalpha)
            if left:
                epss = oldepss + mm(oldalpha,tmp,oldbeta)
            else:
                epss = oldepss + mm(oldbeta,tmp,oldalpha)
            LopezConvTest=N.max(abs(alpha)+abs(beta))
            if LopezConvTest<1.0e-40:
                gs=LA.inv(ee*S-epss)
                if left:
                    test=ee*S-H-mm(ee*dag(S01)-dag(H01),gs,ee*S01-H01)
                else:
                    test=ee*S-H-mm(ee*S01-H01,gs,ee*dag(S01)-dag(H01))
                myConvTest=N.max(abs(mm(test,gs)-N.identity((self.HS.nuo),N.complex)))
                if myConvTest<1.0e-5: # THF: tolerance slightly raised from originally 2.0e-7
                    converged=True
                    if myConvTest>1.0e-8 and left:
                        print "WARNING: Lopez-scheme not-so-well converged for LEFT electrode at E = %.4f eV:"%ee, myConvTest
                    if myConvTest>1.0e-8 and not left:
                        print "WARNING: Lopez-scheme not-so-well converged for RIGHT electrode at E = %.4f eV:"%ee, myConvTest
                else:
                    print "Error: gs iteration: ", iteration
                    print "Lopez report conv : ",LopezConvTest," but not converged :",myConvTest
                    kuk
        return gs        
        
    def setupHS(self,kpoint):
        """
        Setup H, S, H01 and S01 where H01 has large elements in the lower left corner, i.e., H01 = Hi,i+1
        (... 0     H01^+ H     H01   0      ...  )
        (... 0     0     H01^+ H     H01    0      ...  )
        """
        # Save time by not repeating too often
        if N.max(abs(kpoint-self.kpoint))>1e-10:
            self.kpoint=kpoint.copy()
            # Do the trick:
            # H(k=0)+H(kz=0.5) = H + H01 + H10 + H - H01 - H10 = 2 H 
            kp=N.zeros((3),N.float)
            kp[0:2]=kpoint
            self.HS.setkpoint(kp)
            tmpH1, tmpS1 = self.HS.H, self.HS.S
            kp[2]=0.5
            self.HS.setkpoint(kp)
            tmpH2, tmpS2 = self.HS.H, self.HS.S
            self.H, self.S = 0.5*(tmpH1+tmpH2), 0.5*(tmpS1+tmpS2)

            # Additional trick:
            # 1: -i*(H(kz=0.25)-H) = -i*(H + i*H01 - i*H10-H) = H01-H10 
            # 2: H(kz=0)-H  = H + H01 + H10 - H =  H01+H10
            # -> H10 = (-i*(H(kz=0.25)-H) + H(kz=0)-H)/2
            kp[2]=0.25
            self.HS.setkpoint(kp)
            tmpH3, tmpS3 = self.HS.H, self.HS.S
            self.H01, self.S01 = 0.5*(-1j*(tmpH3-self.H)+tmpH1-self.H),\
                0.5*(-1j*(tmpS3-self.S)+tmpS1-self.S)

#############################################################################            
            
class GF:
    def __init__(self,TSHSfile, elecL, elecR, Bulk=False, DeviceAtoms=[0,0]):
        """
        Calculate Green's functions etc for TSHSfile connected to left/right 
        electrode (class surfaceGF). 
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
        DeviceOrbs : Start / end of orbitals. Siesta ordering.
        """
        self.elecL, self.elecR, self.Bulk = elecL, elecR, Bulk
        self.HS = SIO.HS(TSHSfile)

        self.DeviceAtoms=DeviceAtoms
        if DeviceAtoms[0]==0:
            self.DeviceAtoms[0]=1
            self.FoldedL = False
        else:
            self.FoldedL = True
        if DeviceAtoms[1]==0:
            self.DeviceAtoms[1]=self.HS.nua
            self.FoldedR = False
        else:
            self.FoldedR = True
        self.DeviceOrbs = [self.HS.lasto[DeviceAtoms[0]-1]+1, self.HS.lasto[DeviceAtoms[1]]]

        self.nuo0, self.nuoL0, self.nuoR0 = self.HS.nuo, elecL.NA1*elecL.NA2*elecL.HS.nuo, elecR.NA1*elecR.NA2*elecR.HS.nuo 
        self.nuo = self.DeviceOrbs[1]-self.DeviceOrbs[0]+1
        self.nuoL, self.nuoR = self.nuoL0, self.nuoR0 # Not folded, for folded case changed below
        
        print "GF : ",TSHSfile
        print "Device atoms %i-%i, orbitals %i-%i"%(tuple(self.DeviceAtoms+self.DeviceOrbs))
        if not self.FoldedL:
            print "Suggest left folding to atom : ",self.elecL.HS.nua*self.elecL.NA1*self.elecL.NA2+1
        if not self.FoldedR:
            print "Suggest right folding to atom : ",self.HS.nua-self.elecR.HS.nua*self.elecR.NA1*self.elecR.NA2

        if self.FoldedL and self.FoldedR:
            # Check that device region is large enough!
            kpoint=N.zeros((2,),N.float)
            self.setkpoint(kpoint,ispin=0) # At least for one spin
            
            devSt, devEnd = self.DeviceOrbs[0], self.DeviceOrbs[1]
            Soverlap, Hoverlap = N.max(abs(self.S0[0:devSt,devEnd:self.nuo0])), N.max(abs(self.H0[0:devSt,devEnd:self.nuo0]))
            if max(Soverlap,Hoverlap) > 1e-10 :
                print "ERROR! Too much overlap directly from left-top right"
                print "Make Device Region larger!"
                sys.exit(1)
            
            # Find orbitals in device region coupling to left and right.
            tau  = abs(self.S0[0:devSt-1,0:devEnd])
            coupling = N.sum(tau,axis=0)
            ii=devEnd-1
            while coupling[ii]<1e-10: ii=ii-1
            self.devEndL = max(ii+1,self.nuoL0)
            self.nuoL = self.devEndL-devSt+1
            print "Left self energy on orbitals %i-%i"%(devSt,self.devEndL)

            tau  = abs(self.S0[devEnd-1:self.nuo0,0:self.nuo0])
            coupling = N.sum(tau,axis=0)
            ii=devSt-1
            while coupling[ii]<1e-10: ii=ii+1
            self.devStR = min(ii+1,self.nuo0-self.nuoR0+1)
            self.nuoR = devEnd-self.devStR+1
            print "Right self energy on orbitals %i-%i"%(self.devStR,devEnd)


    def calcGF(self,ee,kpoint,ispin=0):
        "Calculate GF etc at energy ee and 2d k-point"

        nuo, nuoL, nuoR = self.nuo, self.nuoL, self.nuoR
        nuo0, nuoL0, nuoR0 = self.nuo0, self.nuoL0, self.nuoR0
        FoldedL, FoldedR = self.FoldedL, self.FoldedR
        devSt, devEnd = self.DeviceOrbs[0], self.DeviceOrbs[1]
        self.setkpoint(kpoint,ispin=ispin)

        # Calculate Sigma without folding
        SigL0 = self.elecL.getSig(ee,kpoint,left=True,Bulk=self.Bulk,ispin=ispin)
        SigR0 = self.elecR.getSig(ee,kpoint,left=False,Bulk=self.Bulk,ispin=ispin)
        
        if FoldedL:
            # Fold down from nuoL0 to the device region
            # A11 A12     g11 g12    I 0
            # A21 A22  *  g21 g22  = 0 I ->
            # g22 = (A22-A21.A11^-1.A12)^-1 ->
            # Sigma = A21.A11^-1.A12          (tau=A12)
            
            devEndL = self.devEndL

            # Do folding
            eSmH = ee*self.S0-self.H0                                        
            eSmHmS = eSmH[0:devEndL,0:devEndL].copy()                             
            if self.Bulk:
                eSmHmS[0:nuoL0,0:nuoL0] = SigL0     
            else:
                eSmHmS[0:nuoL0,0:nuoL0] = eSmHmS[0:nuoL0,0:nuoL0]-SigL0     
            tau  = eSmHmS[0:devSt-1,devSt-1:devEndL].copy()
            taud = eSmHmS[devSt-1:devEndL,0:devSt-1].copy()
            inv = LA.inv(eSmHmS[0:devSt-1,0:devSt-1])
            eSmHmS[devSt-1:devEndL,devSt-1:devEndL]=eSmHmS[devSt-1:devEndL,devSt-1:devEndL]-\
                mm(taud,inv,tau)
            self.SigL = eSmH[devSt-1:devEndL,devSt-1:devEndL]-eSmHmS[devSt-1:devEndL,devSt-1:devEndL]
        else:
            self.SigL=SigL0

        if FoldedR:
            # Fold down from nuoR0 to the device region
            
            devStR = self.devStR

            eSmH = ee*self.S0-self.H0                      
            eSmHmS = eSmH[devStR-1:nuo0,devStR-1:nuo0].copy()
            tmpnuo=len(eSmHmS)                             
            if self.Bulk:
                eSmHmS[tmpnuo-nuoR0:tmpnuo,tmpnuo-nuoR0:tmpnuo] = SigR0     
            else:
                eSmHmS[tmpnuo-nuoR0:tmpnuo,tmpnuo-nuoR0:tmpnuo] = eSmHmS[tmpnuo-nuoR0:tmpnuo,tmpnuo-nuoR0:tmpnuo]-SigR0     
            tau  = eSmHmS[0:nuoR,nuoR:tmpnuo].copy()
            taud = eSmHmS[nuoR:tmpnuo,0:nuoR].copy()
            inv = LA.inv(eSmHmS[nuoR:tmpnuo,nuoR:tmpnuo])
            eSmHmS[0:nuoR,0:nuoR]=eSmHmS[0:nuoR,0:nuoR]-mm(tau,inv,taud)
            self.SigR = eSmH[devStR-1:devEnd,devStR-1:devEnd]-eSmHmS[0:nuoR,0:nuoR]
        else:
            self.SigR=SigR0

        self.GamL, self.GamR = 1.0j*(self.SigL-dag(self.SigL)), 1.0j*(self.SigR-dag(self.SigR))

        # Finally ready to calculate Gr
        eSmH=ee*self.S-self.H
        if FoldedL:
            eSmH[0:nuoL,0:nuoL]=eSmH[0:nuoL,0:nuoL]-self.SigL
        else:
            if self.Bulk:
                eSmH[0:nuoL,0:nuoL]=self.SigL
            else:
                eSmH[0:nuoL,0:nuoL]=eSmH[0:nuoL,0:nuoL]-self.SigL
        if FoldedR:
            eSmH[nuo-nuoR:nuo,nuo-nuoR:nuo]=eSmH[nuo-nuoR:nuo,nuo-nuoR:nuo]-self.SigR
        else:
            if self.Bulk:
                eSmH[nuo-nuoR:nuo,nuo-nuoR:nuo]=self.SigR
            else:
                eSmH[nuo-nuoR:nuo,nuo-nuoR:nuo]=eSmH[nuo-nuoR:nuo,nuo-nuoR:nuo]-self.SigR
        self.Gr= LA.inv(eSmH)
        
    def setkpoint(self,kpoint,ispin=0):
        # Initiate H, S to correct kpoint
        nuo, nuoL, nuoR = self.nuo0, self.nuoL0, self.nuoR0

        kpoint3 = N.zeros((3),N.float)
        kpoint3[0:2]=kpoint[:]
        self.HS.setkpoint(kpoint3)
        # Remove PBC in z-direction
        if self.HS.gamma:
            self.H0=self.HS.H[ispin,:,:].copy()
            self.S0=self.HS.S.copy()
            # Remove direct left/right coupling 
            self.H0[0:nuoL,nuo-nuoR:nuo]=N.zeros((nuoL,nuoR),N.complex)
            self.H0[nuo-nuoR:nuo,0:nuoL]=N.zeros((nuoR,nuoL),N.complex)
            self.S0[0:nuoL,nuo-nuoR:nuo]=N.zeros((nuoL,nuoR),N.complex)
            self.S0[nuo-nuoR:nuo,0:nuoL]=N.zeros((nuoR,nuoL),N.complex)
        else:
            # Do trick with kz
            tmpH1, tmpS1 = self.HS.H[ispin,:,:].copy(), self.HS.S.copy()
            kpoint3[2]=0.5
            self.HS.setkpoint(kpoint3)
            tmpH2, tmpS2 = self.HS.H[ispin,:,:].copy(), self.HS.S.copy()
            self.H0, self.S0 = 0.5*(tmpH1+tmpH2), 0.5*(tmpS1+tmpS2)
        
        if self.FoldedL or self.FoldedR:
            devSt,devEnd=self.DeviceOrbs[0],self.DeviceOrbs[1]
            self.H = self.H0[devSt-1:devEnd,devSt-1:devEnd]
            self.S = self.S0[devSt-1:devEnd,devSt-1:devEnd]
        else:
            self.H, self.S = self.H0, self.S0

    def calcT(self):
        # Calculate transmission
        # Note that size of matrices not uniform and care is taken to minimize computation times
        GamL, GamR, Gr = self.GamL, self.GamR, self.Gr
        nuo, nuoL, nuoR = self.nuo, self.nuoL, self.nuoR

        tmp=mm(GamL,Gr[0:nuoL,nuo-nuoR:nuo])
        tmp=mm(tmp,GamR)
        tmp2=dag(Gr)
        tmp=mm(tmp,tmp2[nuo-nuoR:nuo,0:nuoL])
        Trans= N.trace(tmp)
        if Trans.imag>1e-10: 
            print "Error transmission has large imaginary value :", Trans
            kuk
        return Trans.real

#############################################################################            

if __name__ == '__main__':
    main()
    #profile.run('main()')
