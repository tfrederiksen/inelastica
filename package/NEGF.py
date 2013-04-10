print "SVN $Id$"

import SiestaIO as SIO
import MiscMath as MM
import numpy as N
import numpy.linalg as LA
import sys, string


################## Help functions ############################
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


class ElectrodeSelfEnergy:
    """ 
    Calculate surface Greensfunction and self energy
    (should probably be renamed selfEnergy ...)
    For spinpolarized use the ispin given, for nonpolarized use 
    the same self-energy for both spin 
    """
    def __init__(self,fn,NA1,NA2,voltage=0.0,UseF90helpers=True):
        self.HS=SIO.HS(fn,UseF90helpers=UseF90helpers)
        if self.HS.gamma:
            print "Are you trying to sneak a Gamma point electrode calculation past me?"
            kuk
        self.NA1=NA1
        self.NA2=NA2
        self.kpoint = N.array([1e10,1e10],N.float)
        self.voltage = voltage

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
        The voltage is assumed to be applied symmetrically and just shifts the energies of the self-energy
        """

        eeshifted = ee-self.voltage # Shift of self energies due to voltage
        
        if ispin>=self.HS.nspin:
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
                g0=self.getg0(eeshifted,kpoint,left=left,ispin=ispin)             
                matESmH = eeshifted*self.S-self.H[ispin,:,:]
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
        return self.calcg0_old(ee,left=left,ispin=ispin)
        #Potentially faster method but seems to have numerical instability
        #if hasSciPy and :
        #    return self.calcg0(ee,left=left,ispin=ispin)
        #else:
        #    return self.calcg0_old(ee,left=left,ispin=ispin)
            
    def calcg0(self,ee,ispin=0,left=True):
        # Calculate surface Green's function
        # Euro Phys J B 62, 381 (2008)
        # Inverse of : NOTE, setup for "right" lead.
        # e-h00 -h01  ...
        # -h10  e-h00 ...
        h00,s00,h01,s01 = self.H[ispin,:,:],self.S,self.H01[ispin,:,:],self.S01
        NN, ee = len(h00), N.real(ee)+N.max([N.imag(ee),1e-8])*1.0j
        if left:
            h01, s01 = MM.dagger(h01), MM.dagger(s01)

        # Solve generalized eigen-problem
        # ( e I - h00 , -I) (eps)          (h01 , 0) (eps)
        # ( h10       ,  0) (xi ) = lambda (0   , I) (xi )
        a, b = N.zeros((2*NN,2*NN),N.complex), N.zeros((2*NN,2*NN),N.complex)
        a[0:NN,0:NN] = ee*s00-h00
        a[0:NN,NN:2*NN] = -N.eye(NN)
        a[NN:2*NN,0:NN] = MM.dagger(h01)-ee*MM.dagger(s01)
        b[0:NN,0:NN] = h01-ee*s01
        b[NN:2*NN,NN:2*NN] = N.eye(NN)
        ev, evec = SLA.eig(a,b)

        # Select lambda <0 and the eps part of the evec
        ipiv = N.where(N.abs(ev)<1.0)[0]
        ev, evec = ev[ipiv], N.transpose(evec[:NN,ipiv])        
        # Normalize evec
        norm = N.sqrt(N.diag(MM.mm(evec,MM.dagger(evec))))
        evec = MM.mm(N.diag(1.0/norm),evec)

        # E^+ Lambda_+ (E^+)^-1 --->>> g00
        EP = N.transpose(evec)
        FP = MM.mm(EP,N.diag(ev),LA.inv(MM.mm(MM.dagger(EP),EP)),MM.dagger(EP))
        g00 = LA.inv(ee*s00-h00-MM.mm(h01-ee*s01,FP))

        # Check!
        err=N.max(N.abs(g00-LA.inv(ee*s00-h00-\
                         MM.mm(h01-ee*s01,g00,MM.dagger(h01)-ee*MM.dagger(s01)))))
        if err>1.0e-8 and left:
            print "WARNING: Lopez-scheme not-so-well converged for LEFT electrode at E = %.4f eV:"%ee, err
        if err>1.0e-8 and not left:
            print "WARNING: Lopez-scheme not-so-well converged for RIGHT electrode at E = %.4f eV:"%ee, err
        return g00

    def calcg0_old(self,ee,ispin=0,left=True):
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
        H, S, H01, S01 = self.H[ispin,:,:] ,self.S ,self.H01[ispin,:,:], self.S01

        alpha, beta = MM.dagger(H01)-ee*MM.dagger(S01), H01-ee*S01
        eps, epss = H.copy(), H.copy()
        
        converged=False
        iteration=0
        while not converged:
            iteration+=1
            oldeps, oldepss = eps.copy(), epss.copy()
            oldalpha, oldbeta = alpha.copy(), beta.copy()
            tmpa=LA.solve(ee*S - oldeps,oldalpha)
            tmpb=LA.solve(ee*S - oldeps,oldbeta)
            alpha, beta = MM.mm(oldalpha,tmpa), MM.mm(oldbeta,tmpb)
            eps = oldeps + MM.mm(oldalpha,tmpb)+MM.mm(oldbeta,tmpa)
            if left:
                epss = oldepss + MM.mm(oldalpha,tmpb)
            else:
                epss = oldepss + MM.mm(oldbeta,tmpa)
            LopezConvTest=N.max(abs(alpha)+abs(beta))
            if LopezConvTest<1.0e-40:
                gs=LA.inv(ee*S-epss)
                if left:
                    test=ee*S-H-MM.mm(ee*MM.dagger(S01)-MM.dagger(H01),gs,ee*S01-H01)
                else:
                    test=ee*S-H-MM.mm(ee*S01-H01,gs,ee*MM.dagger(S01)-MM.dagger(H01))
                myConvTest=N.max(abs(MM.mm(test,gs)-N.identity((self.HS.nuo),N.complex)))
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
    def __init__(self, TSHSfile, elecL, elecR, Bulk=True, DeviceAtoms=[0,0]):
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
        # Quantities expressed in nonorthogonal basis:
        self.OrthogonalDeviceRegion = False

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
                MM.mm(taud,inv,tau)
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
            eSmHmS[0:nuoR,0:nuoR]=eSmHmS[0:nuoR,0:nuoR]-MM.mm(tau,inv,taud)
            self.SigR = eSmH[devStR-1:devEnd,devStR-1:devEnd]-eSmHmS[0:nuoR,0:nuoR]
        else:
            self.SigR=SigR0

        self.GamL, self.GamR = 1.0j*(self.SigL-MM.dagger(self.SigL)), 1.0j*(self.SigR-MM.dagger(self.SigR))

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
        self.Gr = LA.inv(eSmH)
        self.Ga = MM.dagger(self.Gr)
        
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
        self.OrthogonalDeviceRegion = False


    def calcT(self,channels=10):
        # Calculate transmission
        # Note that size of matrices may not be not uniform and care is taken to minimize computation times
        GamL, GamR, Gr = self.GamL, self.GamR, self.Gr
        if GamL.shape == Gr.shape and GamR.shape == Gr.shape:
            # Quantities have same shape and the matrix product is
            # straight-forward (but perhaps unnecesarily time consuming)
            Tmat = MM.mm(GamL,Gr,GamR,MM.dagger(Gr))
        else:
            # Nonorthogonal quantities, computation time minimized
            nuo, nuoL, nuoR = self.nuo, self.nuoL, self.nuoR
            tmp = MM.mm(GamL,Gr[0:nuoL,nuo-nuoR:nuo])
            tmp = MM.mm(tmp,GamR)
            tmp2 = MM.dagger(Gr)
            Tmat = MM.mm(tmp,tmp2[nuo-nuoR:nuo,0:nuoL])
        Trans = N.trace(Tmat)
        if Trans.imag>1e-10: 
            print "Error transmission has large imaginary value :", Trans
            kuk
        # Calculate eigenchannel transmissions too
        tval,tvec = LA.eig(Tmat)
        tval = sorted(tval,reverse=True) # Sort eigenvalues descending
        T = [Trans.real]
        for i in range(channels):
            T += [tval[i].real]
        return N.array(T)

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
        self.S, self.H = MM.mm(Us,self.S,Us), MM.mm(Us,self.H,Us)

        # Sigmas/Gammas in pyTBT GF can be smaller than device region
        # First give them the shape of the device region
        nnL, nnR = len(self.SigL),  len(self.SigR)
        S1, S2 = N.zeros(self.H.shape, N.complex), N.zeros(self.H.shape, N.complex)
        S1[0:nnL, 0:nnL] ,  S2[-nnR:, -nnR:] = self.SigL, self.SigR
        # Resetting Sigmas to orthogonalized quantities
        self.SigL, self.SigR = MM.mm(Us,S1,Us), MM.mm(Us,S2,Us)
        # ... now the same for the Gammas
        G1, G2 = N.zeros(self.H.shape, N.complex), N.zeros(self.H.shape, N.complex)
        G1[0:nnL, 0:nnL] ,  G2[-nnR:, -nnR:] = self.GamL, self.GamR
        # Resetting Gammas to orthogonalized quantities
        self.GamL, self.GamR = MM.mm(Us,G1,Us), MM.mm(Us,G2,Us)
        
        # Orthogonalize Greens functions
        self.Gr = MM.mm(Usi,self.Gr,Usi)
        self.Ga = MM.dagger(self.Gr)

    def __calcEigChan(self,A1,G2,channels=10):
        # Calculate Eigenchannels using recipe from PRB
        # For right eigenchannels, A1=A2, G2=G1 !!!
        ev, U = LA.eigh(A1)
        U = N.transpose(U)

        Utilde=0.0*U

        for jj in range(len(ev)): # Problems with negative numbers
            if ev[jj]<0: ev[jj]=0
            Utilde[jj,:]=N.sqrt(ev[jj]/(2*N.pi))*U[jj,:]
        Utilde=N.transpose(Utilde)

        tt=MM.mm(MM.dagger(Utilde),2*N.pi*G2,Utilde)
        evF, UF = LA.eigh(tt)
        UF = N.transpose(UF)

        EC=[]
        for jj in range(channels):
            tmp=MM.mm(self.Us,Utilde,UF[len(evF)-jj-1,:])
            EC.append(tmp.copy())

        return EC, evF

    def calcEigChan(self,channels=10):
        if not self.OrthogonalDeviceRegion:
            self.orthogonalize()
        # Calculate Eigenchannels from left
        self.A1 = MM.mm(self.Gr,self.GamL,self.Ga)
        self.ECleft, self.EigTleft = self.__calcEigChan(self.A1,self.GamR,channels)
        print 'Left eigenchannel transmissions:',self.EigTleft
        # Calculate Eigenchannels from right
        self.A2 = MM.mm(self.Gr,self.GamR,self.Ga)
        self.ECright, self.EigTright = self.__calcEigChan(self.A2,self.GamL,channels)
        print 'Right eigenchannel transmissions:',self.EigTright



#############################################################################            

