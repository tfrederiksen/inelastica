version = "SVN $Id$"
print version

"""
A completely rewritten Phonons.py script with various improvements
and additional flexibility:

* The dynamic atoms no longer need to be a complete range(FCfirst,FClast+1).
  Any arbitrary list of atoms (options.DynamicAtoms) can now be specified.

* The displacement amplitude in the various FCrun and OSrun directories
  may correspond to different values.

* The code accepts that some atoms may have been displaced in several FCrun directories.
  Only the first instance (first FCrun directory) encountered is read/used.

* The auxiliary NetCDF file has been eliminated by simply interchanging 
  the loops over gradients and phonon modes. Thus, only one gradient
  need to be present in the memory at one time.

Thomas Frederiksen, August 2014
"""

import SiestaIO as SIO
import Symmetry
import CommonFunctions as CF
import MakeGeom as MG
import PhysicalConstants as PC
import MiscMath as MM
import WriteNetCDF as NCDF
import ValueCheck as VC

import numpy as N
import numpy.linalg as LA
import glob, os,sys,string

vinfo = [version,SIO.version,Symmetry.version,CF.version,
         MG.version,PC.version,MM.version,NCDF.version,VC.version]

def GetOptions(argv,**kwargs):
    # if text string is specified, convert to list
    if type(argv)==type(''): argv = argv.split()

    import optparse as o

    usage = "usage: %prog [options] DestinationDirectory"
    description = "Methods to calculate vibrations and e-ph couplings from SIESTA output"

    p = o.OptionParser(description=description,usage=usage)

    p.add_option("-c", "--CalcCoupl",dest="CalcCoupl",
                 help="Calculate e-ph couplings [default=%default]",
                 action="store_true",default=False)
    
    p.add_option("-F","--DeviceFirst",dest="DeviceFirst",
                 help="First device atom index (in the electronic basis) [default=%default]",
                 type="int",default=1)
    p.add_option("-L","--DeviceLast",dest="DeviceLast",
                 help="Last device atom index (in the electronic basis) [default=%default]",
                 type="int",default=1000)
    
    p.add_option("--FCfirst",dest="FCfirst",
                 help="First FC atom index [default=%default]",
                 type="int",default=1)
    p.add_option("--FClast", dest="FClast",
                 help="Last FC atom index [default=%default]" ,
                 type="int",default=1000)
    
    p.add_option("--PBCFirst", dest="PBCFirst",\
                 help="For eliminating interactions through periodic boundary conditions in z-direction [default=%default]",
                 type="int",default=1)
    p.add_option("--PBCLast", dest="PBCLast",\
                 help="For eliminating interactions through periodic boundary conditions in z-direction [default=%default]",
                 type="int",default=1000)
    
    p.add_option("--FCwildcard",dest="FCwildcard",
                 help="Wildcard for FC directories [default=%default]",
                 type="str",default="./FC*")
    
    p.add_option("--OSdir",dest="onlySdir",
                 help="Location of OnlyS directory [default=%default]",
                 type="str",default="./OSrun")
    
    p.add_option("-a", "--AbsoluteEnergyReference",dest="AbsEref",
                help="Use an absolute energy reference (Fermi energy of equilibrium structure) for displaced Hamiltonians (e.g., when eF is not well-defined) instead of the instantaneous Fermi energy for the displaced geometries, cf. Eq.(17) in PRB 75, 205413 (2007) [default=%default]",action="store_true",default=False)
    
    p.add_option("-i", "--Isotopes",dest="Isotopes",
                 help="List of substitutions [[i1, anr1],...], where atom index i1 (SIESTA numbering) is set to be of type anr1. Alternatively, the argument can be a file with the input string [default=%default]",default=[])

    p.add_option("--AtomicMass", dest='AtomicMass', default='[]',
                 help="Option to add to (or override!) existing dictionary of atomic masses. Format is a list [[anr1,mass1(,label)],...] [default=%default]")
    
    p.add_option("-x","--k1", dest='k1', default=0.0,type='float',
                 help="k-point along a1 where e-ph couplings are evaluated [%default]")
    p.add_option("-y","--k2", dest='k2', default=0.0,type='float',
                 help="k-point along a2 where e-ph couplings are evaluated [%default]")
    p.add_option("-z","--k3", dest='k3', default=0.0,type='float',
                 help="k-point along a3 where e-ph couplings are evaluated [%default]")

    (options, args) = p.parse_args(argv)

    # Get the last positional argument
    options.DestDir = VC.GetPositional(args,"You need to specify a destination directory")

    # With this one can overwrite the logging information
    if "log" in kwargs:
        options.Logfile = kwargs["log"]
    else:
        options.Logfile = 'Phonons.log'

    # k-point
    options.kpoint = N.array([options.k1,options.k2,options.k3],N.float)
    del options.k1,options.k2,options.k3

    # Dynamic atoms
    options.DynamicAtoms = range(options.FCfirst,options.FClast+1)
    del options.FCfirst, options.FClast
 
    # PBCFirst/PBCLast
    if options.PBCFirst<options.DeviceFirst:
        options.PBCFirst = options.DeviceFirst
    if options.PBCLast>options.DeviceLast:
        options.PBCLast = options.DeviceLast

    # Isotopes specified in separate file?
    if type(options.Isotopes)!=type([]): # i.e., not a list
        if os.path.isfile(options.Isotopes):
            f = open(options.Isotopes)
            s = ''
            for line in f.readlines():
                s += line.replace('\n','').replace(' ','')
            options.Isotopes = s
            f.close()

    # Check if AtomicMasses are specified
    if options.AtomicMass!='[]':
        from Inelastica import PhysicalConstants as PC
        masslist = eval(options.AtomicMass.replace('\n','').replace(' ',''))
        for elm in masslist:
            anr = int(elm[0])
            mass = float(elm[1])
            PC.AtomicMass[anr] = mass
            if len(elm)==3:
                label = elm[2]
                PC.PeriodicTable[anr] = label
                PC.PeriodicTable[label] = anr
        print 'AtomicMass =',PC.AtomicMass
        print 'PeriodicTable =',PC.PeriodicTable

    return options
  

class FCrun():

    def __init__(self,runfdf):
        self.fdf = runfdf
        self.directory,self.tail =  os.path.split(runfdf)
        self.systemlabel = SIO.GetFDFlineWithDefault(runfdf,'SystemLabel', str, 'siesta','Phonons')
        FCfirst = SIO.GetFDFlineWithDefault(runfdf,'MD.FCfirst',int,0,'Phonons')
        FClast = SIO.GetFDFlineWithDefault(runfdf,'MD.FClast',int,0,'Phonons')
        # Finite-displacement amplitude
        ampl,unit = SIO.GetFDFline(runfdf,KeyWord='MD.FCDispl')
        if unit.upper()=='ANG':
            self.Displ = float(ampl)
        elif unit.upper()=='BOHR':
            self.Displ = float(ampl)*PC.Bohr2Ang
        print 'Displacement = %.6f Ang'%self.Displ
        # Read geometry
        self.geom = MG.Geom(runfdf)
        # Compare with XV file corrected for last displacement
        XV = self.directory+'/%s.XV'%self.systemlabel
        geomXV = MG.Geom(XV)
        geomXV.xyz[FClast-1,2] -= self.Displ
        if not N.allclose(geomXV.xyz,self.geom.xyz):
            sys.exit('Error: Geometries %s and %s should differ ONLY by displacement of atom %s in z'\
                     %(runfdf,XV,FClast))
        # Set up FC[i,a,j,b]: Force constant (eV/A^2) from moved atom i, axis a to atom j, axis b
        natoms = self.geom.natoms
        self.m = N.zeros((FClast-FCfirst+1,3,natoms,3), N.float)
        self.p = N.zeros((FClast-FCfirst+1,3,natoms,3), N.float)
        fc = N.array(SIO.ReadFCFile(self.directory+'/%s.FC'%self.systemlabel))
        for i in range(FClast-FCfirst+1):
            for j in range(3):
                self.m[i,j] = fc[2*(3*i+j)*natoms:(2*(3*i+j)+1)*natoms]
                self.p[i,j] = fc[(2*(3*i+j)+1)*natoms:(2*(3*i+j)+2)*natoms]
        # Correct force constants for the moved atom
        # Cf. Eq. (13) in Frederiksen et al. PRB 75, 205413 (2007)
        for i in range(FClast-FCfirst+1):
            for j in range(3):
                self.m[i,j,FCfirst-1+i,:] = 0.0
                self.m[i,j,FCfirst-1+i,:] = -N.sum(self.m[i,j],axis=0)
                self.p[i,j,FCfirst-1+i,:] = 0.0
                self.p[i,j,FCfirst-1+i,:] = -N.sum(self.p[i,j],axis=0)
        # Determine TSHS files
        files = glob.glob(self.directory+'/%s*.TSHS'%self.systemlabel)
        files.sort()
        if (FClast-FCfirst+1)*6+1 != len(files):
            sys.exit('Phonons.GetFileLists: WARNING - Inconsistent number of *.TSHS files in %s'%self.directory)
        self.DynamicAtoms = range(FCfirst,FClast+1)
        # Build dictionary over TSHS files and corresponding displacement amplitudes
        self.TSHS = {}
        self.TSHS[0] = files[0] # Equilibrium TSHS
        for i,v in enumerate(self.DynamicAtoms):
            for j in range(3):
                # Shifted TSHS files (atom,axis,direction)
                self.TSHS[v,j,-1] = files[1+6*i+2*j]
                self.TSHS[v,j,1] = files[1+6*i+2*j+1]

    def GetOrbitalIndices(self):
        import Scientific.IO.NetCDF as nc
        # Determine snr (siesta number) for each label
        csl = SIO.GetFDFblock(self.fdf, KeyWord = 'ChemicalSpeciesLabel')
        csl2snr = {}
        for set in csl:
            csl2snr[set[2]] = set[0]
        # Determine nao (number of orbitals) for each snr
        ionNCfiles = glob.glob(self.directory+'/*.ion.nc*')
        snr2nao = {}
        for ionfile in ionNCfiles:
            if ionfile.endswith('.gz'):
                print 'Phonons.GetOrbitalIndices: Unzipping',ionfile
                os.system('gunzip '+ionfile)
                ionfile = ionfile[:-3]
            file = nc.NetCDFFile(ionfile,'r')
            thissnr = csl2snr[file.Label]
            snr2nao[int(thissnr)] = int(file.Number_of_orbitals[0])
            file.close()
        print 'Phonons.GetOrbitalIndices: Dictionary snr2nao =',snr2nao
        # Determine which orbital indices that belongs to a certain atom
        orbitalIndices = []
        tmpOrb = 0
        for num in self.geom.snr:
            nao = snr2nao[num]
            orbitalIndices.append([tmpOrb,tmpOrb+int(nao)-1])
            tmpOrb+=nao
        self.orbitalIndices = N.array(orbitalIndices)
        self.nao = tmpOrb # total number of orbitals
        self.snr2nao = snr2nao # orbitals per species
        return self.orbitalIndices,self.nao

class OSrun:
    def __init__(self,onlySdir,kpoint):
        print 'Phonons.GetOnlyS: Reading from', onlySdir
        onlySfiles = glob.glob(onlySdir+'/*.onlyS*')
        onlySfiles.sort()
        if len(onlySfiles)<1:
            sys.exit('Phonons.GetOnlyS: No .onlyS file found!')
        if len(onlySfiles)!=6:
            sys.exit('Phonons.GetOnlyS: Wrong number of onlyS files found!')
        else:
            onlyS = {}
            Displ = {}
            for file in onlySfiles:
                thisHS = SIO.HS(file)
                thisHS.setkpoint(kpoint)
                S = thisHS.S
                del thisHS
                nao = len(S)/2
                S0 = S[0:nao,0:nao].copy()
                dmat=S[0:nao,nao:nao*2].copy()
                if file.endswith('_1.onlyS'):
                    onlyS[0,-1] = dmat
                elif file.endswith('_2.onlyS'):
                    onlyS[0,1] = dmat
                elif file.endswith('_3.onlyS'):
                    onlyS[1,-1] = dmat
                elif file.endswith('_4.onlyS'):
                    onlyS[1,1] = dmat
                elif file.endswith('_5.onlyS'):
                    onlyS[2,-1] = dmat
                elif file.endswith('_6.onlyS'):
                    onlyS[2,1] = dmat
            # Loop over the 6 doubled geometries and determine the displacement
            for i in range(1,7):
                thisd = 1e10
                xyz = N.array(SIO.Getxyz(onlySdir+'/RUN_%i.fdf'%i))
                for j in range(1,len(xyz)):
                    thisd = min(thisd,(N.dot(xyz[0]-xyz[j],xyz[0]-xyz[j]))**.5)
                Displ[(i-1)/2,1-2*(i%2)] = thisd
                print 'Phonons.GetOnlyS: OnlyS-displacement (min) = %.5f Ang'%thisd
            # Construct dS array
            self.S0 = S0
            self.dS = N.empty((3,)+dmat.shape,dtype=dmat.dtype)
            for j in range(3):
                self.dS[j] = (onlyS[j,1]-onlyS[j,-1])/(Displ[j,-1]+Displ[j,1])
            self.Displ = Displ


class DynamicalMatrix():

    def __init__(self,fdfs,DynamicAtoms=None):
        self.fdfs = fdfs
        self.FCRs = [FCrun(f) for f in fdfs]
        self.geom = self.FCRs[0].geom # assume identical geometries
        if DynamicAtoms:
            self.SetDynamicAtoms(DynamicAtoms)

    def SetDynamicAtoms(self,DynamicAtoms):        
        self.DynamicAtoms = DynamicAtoms
        NN = len(DynamicAtoms)
        self.m = N.zeros((NN,3,self.geom.natoms,3),N.complex)
        self.p = N.zeros((NN,3,self.geom.natoms,3),N.complex)
        self.Displ = {}
        self.TSHS = {}
        self.TSHS[0] = self.FCRs[0].TSHS[0]
        for i,v in enumerate(DynamicAtoms):
            for fcr in self.FCRs:
                if v in fcr.DynamicAtoms:
                    j = fcr.DynamicAtoms.index(v)
                    print 'Reading FC data for dynamic atom %i from %s' %(v,fcr.fdf)
                    self.Displ[v] = fcr.Displ
                    self.m[i] = fcr.m[j]
                    self.p[i] = fcr.p[j]
                    for k in range(3):
                        self.TSHS[v,k,-1] = fcr.TSHS[v,k,-1]
                        self.TSHS[v,k,1] = fcr.TSHS[v,k,1]
                    break
            # Check that we really found the required atom
            if len(self.Displ)<=i:
                sys.exit('Error: Did not find FC data for a dynamic atom %i'%v)
        self.mean = (self.m+self.p)/2

    def SetMasses(self,Isotopes=[]):
        self.Masses = []
        # Set default masses
        for i,v in enumerate(self.DynamicAtoms):
            self.Masses.append( PC.AtomicMass[self.geom.anr[v-1]] )
        # Override with specified Isotopes?
        for ii,anr in Isotopes:
            if ii in self.DynamicAtoms:
                j = self.DynamicAtoms.index(ii)
                print 'Phonons.Analyse: Isotope substitution for atom %i (SIESTA numbering):'%ii
                print '  ... atom type %i --> %i'%(self.geom.anr[ii-1],anr)
                print '  ... atom mass %.4f --> %.4f'%(PC.AtomicMass[self.geom.anr[ii-1]],\
                                                       PC.AtomicMass[anr])
                self.Masses[j] = PC.AtomicMass[anr]

    def ApplySumRule(self,FC):
        FC0 = FC.copy()
        # Correct force constants for the moved atom
        # Cf. Eq. (13) in Frederiksen et al. PRB 75, 205413 (2007)
        for i,v in enumerate(self.DynamicAtoms):
            for j in range(3):
                FC[i,j,v-1,:] = 0.0
                FC[i,j,v-1,:] = -N.sum(FC[i,j],axis=0)
        print 'Total sumrule change in FC: %.3e eV/Ang' % N.sum(abs(FC0)-abs(FC))
        return FC

    def ComputePhononModes(self,FC):
        dyn = len(self.DynamicAtoms)
        FCtilde = N.zeros((dyn,3,dyn,3),N.complex)
        # Symmetrize and mass-scale
        for i,v in enumerate(self.DynamicAtoms):
            for j,w in enumerate(self.DynamicAtoms):
                FCtilde[i,:,j,:] = 0.5*(FC[i,:,w-1,:]+MM.dagger(FC[j,:,v-1,:]))\
                                   /(self.Masses[i]*self.Masses[j])**0.5
        # Solve eigenvalue problem with symmetric FCtilde
        FCtilde = FCtilde.reshape((3*dyn,3*dyn),order='C')
        evalue,evec = LA.eigh(FCtilde)
        #evalue,evec = LA.eig(FCtilde)
        evec = N.transpose(evec)
        evalue = N.array(evalue,N.complex)
        # Calculate frequencies
        const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
        hw = const*evalue**0.5 # Units in eV
        for i in range(3*dyn):
            # Real eigenvalues are defined as positive, imaginary eigenvalues as negative
            hw[i] = hw[i].real - abs(hw[i].imag)
        hw = hw.real
        # Normalize eigenvectors
        U = evec.copy()
        for i in range(3*dyn):
            U[i] = U[i]/(N.dot(N.conjugate(U[i]),U[i])**0.5)
        # Sort in order descending mode energies
        #hw = hw[::-1] # reverse array
        #U = U[::-1] # reverse array
        indx = N.argsort(hw)[::-1] # reverse
        hw = hw[indx]
        U = U[indx]
        # Print mode frequencies
        print 'Phonons.CalcPhonons: Frequencies in meV:'
        for i in range(3*dyn):
            print string.rjust('%.3f'%(1000*hw[i]),9),
            if (i-5)%6==0: print
        if (i-5)%6!=0: print
        #print 'Phonons.CalcPhonons: Frequencies in cm^-1:'
        #for i in range(3*dyn):
        #    print string.rjust('%.3f'%(hw[i]/PC.invcm2eV),9),
        #    if (i-5)%6==0: print
        #if (i-5)%6!=0: print

        # Compute real displacement vectors
        Udisp = U.copy()
        for i in range(3*dyn):
            # Eigenvectors after division by sqrt(mass)
            Udisp[:,i] = U[:,i]/self.Masses[i/3]**.5

        # Expand vectors to full geometry
        UU = N.zeros((len(hw),self.geom.natoms,3),N.complex)
        UUdisp = N.zeros((len(hw),self.geom.natoms,3),N.complex)
        for i in range(len(hw)):
            for j,v in enumerate(self.DynamicAtoms):
                UU[i,v-1,:] = [U[i,3*j],U[i,3*j+1],U[i,3*j+2]]
                UUdisp[i,v-1,:] = [Udisp[i,3*j],Udisp[i,3*j+1],Udisp[i,3*j+2]]
        self.hw = hw        
        self.U = U
        self.Udisp = Udisp
        self.UU = UU
        self.UUdisp = UUdisp

        
    def PrepareGradients(self,onlySdir,kpoint,DeviceFirst,DeviceLast,AbsEref):
        print '\nPhonons.PrepareGradients: Setting up various arrays'
        self.kpoint = kpoint
        self.OrbIndx,nao = self.FCRs[0].GetOrbitalIndices()
        OS = OSrun(onlySdir,kpoint)
        self.dS = OS.dS
        self.TSHS0 = SIO.HS(self.TSHS[0])
        self.TSHS0.setkpoint(kpoint)
        self.invS0H0 = N.empty((2,)+self.TSHS0.H.shape,dtype=self.TSHS0.H.dtype)
        invS0 = LA.inv(OS.S0)
        self.nspin = len(self.TSHS0.H)
        for iSpin in range(self.nspin):
            self.invS0H0[0,iSpin,:,:] = MM.mm(invS0,self.TSHS0.H[iSpin,:,:])
            self.invS0H0[1,iSpin,:,:] = MM.mm(self.TSHS0.H[iSpin,:,:],invS0)
        del invS0
        # don't re-create the array every time... too expensive    
        self.dSdij = N.zeros((nao,nao),N.complex)

        # Take Device region 
        self.DeviceAtoms = range(DeviceFirst,DeviceLast+1)
        first,last = self.OrbIndx[DeviceFirst-1][0],self.OrbIndx[DeviceLast-1][1]
        self.h0 = self.TSHS0.H[:,first:last+1,first:last+1]
        self.s0 = self.TSHS0.S[first:last+1,first:last+1]
        self.DeviceFirst = DeviceFirst
        self.DeviceLast = DeviceLast
        self.AbsEref = AbsEref
        
    def GetGradient(self,Atom,Axis,AbsEref=False):
        print 'Phonons.GetGradient: Computing dH[%i,%i]'%(Atom,Axis)
        # Read TSHS files
        TSHSm = SIO.HS(self.TSHS[Atom,Axis,-1])
        TSHSm.setkpoint(self.kpoint)
        TSHSp = SIO.HS(self.TSHS[Atom,Axis,1])
        TSHSp.setkpoint(self.kpoint)
        # Use Fermi energy of equilibrium calculation as energy reference?
        if self.AbsEref:
            print 'Computing gradient with absolute energy reference'
            for iSpin in range(self.nspin):
                TSHSm.H[iSpin,:,:] += (TSHSm.ef-self.TSHS0.ef)*TSHSm.S
                TSHSp.H[iSpin,:,:] += (TSHSp.ef-self.TSHS0.ef)*TSHSp.S
        # Compute direct gradient
        dH = (TSHSp.H-TSHSm.H)/(2*self.Displ[Atom])
        del TSHSm, TSHSp
        # Orbital range for the displaced atom:
        f,l = self.OrbIndx[Atom-1]
        self.dSdij[:,f:l+1] = self.dS[Axis,:,f:l+1]
        # Apply Pulay-type corrections
        for iSpin in range(self.nspin):
            dH[iSpin,:,:] -= MM.mm(MM.dagger(self.dSdij),self.invS0H0[0,iSpin,:,:]) \
                             + MM.mm(self.invS0H0[1,iSpin,:,:],self.dSdij)
        self.dSdij[:,f:l+1] = 0. # reset
        return dH

    def ComputeEPHcouplings(self,PBCFirst,PBCLast):
        first,last = self.OrbIndx[self.DeviceFirst-1][0],self.OrbIndx[self.DeviceLast-1][1]
        rednao = last+1-first
        const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
        Heph = N.zeros((len(self.hw),self.nspin,rednao,rednao),N.complex)
        # Loop over dynamic atoms
        for i,v in enumerate(self.DynamicAtoms):
            # Loop over axes
            for j in range(3):
                # Compute gradient
                dH = self.GetGradient(v,j)
                # Remove Periodic Boundary terms
                nuo = len(dH[0])
                pbcf = self.OrbIndx[PBCFirst-1][0]
                pbcl  = self.OrbIndx[PBCLast-1][1]
                if v < PBCFirst:
                    # we have something to remove...
                    print 'Warning: Setting certain elements in dH[%i,%i] to zero because %i<PBCFirst'%(v,j,v)
                    #bb = (PBCFirst - FCfirst) * 3
                    dH[:,pbcl+1:nuo,:] = 0.0
                    dH[:,:,pbcl+1:nuo] = 0.0
                if PBCLast < v:
                    print 'Warning: Setting certain elements in dH[%i,%i] to zero because PBCLast<%i'%(v,j,v)
                    #aa = (PBCLast - FCfirst) * 3
                    dH[:,:pbcf-1,:] = 0.0
                    dH[:,:,:pbcf-1] = 0.0
                # Device part
                dh = dH[:,first:last+1,first:last+1]
                # Loop over modes and throw away the gradient (to save memory)
                for m in range(len(self.hw)):
                    Heph[m] += const*dh*self.UU[m,v-1,j]/(2*self.Masses[i]*self.hw[m])**.5
        del dH, dh
        self.heph = Heph

    def WriteOutput(self,label):
        ### Write MKL- and xyz-files
        natoms = self.geom.natoms
        hw = self.hw
        # Write only real part of eigenvectors
        UU = self.UU.reshape(len(hw),3*natoms).real 
        UUdisp = self.UUdisp.reshape(len(hw),3*natoms).real
        
        SIO.WriteMKLFile('%s.mkl'%label,self.geom.anr,self.geom.xyz,hw,UU,1,natoms)
        SIO.WriteMKLFile('%s.real-displ.mkl'%label,self.geom.anr,self.geom.xyz,hw,UUdisp,1,natoms)
        SIO.WriteXYZFile('%s.xyz'%label,self.geom.anr,self.geom.xyz)
        WriteFreqFile('%s.freq'%label,hw)
        WriteVibDOSFile('%s.Gfdos'%label,hw,type='Gaussian')
        WriteVibDOSFile('%s.Lfdos'%label,hw,type='Lorentzian')
        WriteAXSFFiles('%s.mol.axsf'%label,self.geom.xyz,self.geom.anr,hw,UU,1,natoms)
        WriteAXSFFiles('%s.mol.real-displ.axsf'%label,self.geom.xyz,self.geom.anr,hw,UUdisp,1,natoms)
        WriteAXSFFilesPer('%s.per.axsf'%label,self.geom.pbc,self.geom.xyz,self.geom.anr,hw,UU,1,natoms)
        WriteAXSFFilesPer('%s.per.real-displ.axsf'%label,self.geom.pbc,self.geom.xyz,self.geom.anr,\
                             hw,UUdisp,1,natoms)
        # Netcdf format
        NCDF.write('%s.nc'%label,hw,'hw')
        NCDF.write('%s.nc'%label,UU,'U')
        NCDF.write('%s.nc'%label,UUdisp,'Udisp')
        NCDF.write('%s.nc'%label,self.geom.pbc,'CellVectors')
        NCDF.write('%s.nc'%label,self.geom.xyz,'GeometryXYZ')
        NCDF.write('%s.nc'%label,self.geom.anr,'AtomNumbers')
        NCDF.write('%s.nc'%label,self.geom.snr,'SpeciesNumbers')
        NCDF.write('%s.nc'%label,self.Masses,'Masses')
        NCDF.write('%s.nc'%label,self.DynamicAtoms,'DynamicAtoms')
        try:
            NCDF.write('%s.nc'%label,self.DeviceAtoms,'DeviceAtoms')
            NCDF.write('%s.nc'%label,self.kpoint,'kpoint')
            NCDF.write('%s.nc'%label,self.h0.real,'H0')
            NCDF.write('%s.nc'%label,self.s0.real,'S0')
            GammaPoint = N.dot(self.kpoint,self.kpoint)<1e-7
            if not GammaPoint:
                NCDF.write('%s.nc'%label,self.h0.imag,'ImH0')
                NCDF.write('%s.nc'%label,self.s0.imag,'ImS0')
        except:
            print 'Hamiltonian etc not found'
        try:
            NCDF.write('%s.nc'%label,self.heph.real,'He_ph')
            if not GammaPoint:
                NCDF.write('%s.nc'%label,self.heph.imag,'ImHe_ph')
        except:
            print 'EPH couplings etc not found'

def WriteFreqFile(filename,hw):
    print 'Phonons.WriteFreqFile: Writing',filename
    file = open(filename,'w')
    file.write('# ')
    for i in range(len(hw)):
        file.write(' %f '%(1000*hw[i]))
    file.write('\n# index   freq/meV \n')
    for i in range(len(hw)):
        file.write('%i  %f \n'%(i,1000*hw[i]))
    file.close()

def WriteVibDOSFile(filename,hw,type='Gaussian'):
    'Vibrational DOS with Gaussian or Lorentzian broadening'
    fmax = max(hw)
    erng = N.linspace(0,1.2*fmax,1001)
    eta = N.linspace(0.001,0.01,11)
    ERNG = N.outer(erng,0*eta+1.)
    ETA = N.outer(0*erng+1,eta)
    spectrum = N.zeros((len(erng),len(eta)),N.float)
    for i in range(len(hw)):
        if type=='Gaussian':
            spectrum += (2*N.pi)**(-.5)/ETA*N.exp(N.clip(-1.0*(hw[i]-ERNG)**2/(2*ETA**2),-300,300))
            spectrum -= (2*N.pi)**(-.5)/ETA*N.exp(N.clip(-1.0*(-hw[i]-ERNG)**2/(2*ETA**2),-300,300))
        elif type=='Lorentzian':
            spectrum += 1/N.pi*ETA/((hw[i]-ERNG)**2+ETA**2)
            spectrum -= 1/N.pi*ETA/((-hw[i]-ERNG)**2+ETA**2)
    # Write data to file
    print 'Phonons.WriteVibDOSFile: Writing', filename
    f = open(filename,'w')
    f.write('\n# energy/eV  DOS/atom (eta=1,2,3,...,10meV) \n')
    for i in range(len(erng)):
        f.write('%.5e  '%(erng[i]))
        for j in range(len(eta)):
            f.write('%.5e '%(spectrum[i,j]*3/len(hw)))
        f.write('\n')
    f.close()                                                                

def WriteAXSFFiles(filename,xyz,anr,hw,U,FCfirst,FClast):
    'Writes the vibrational normal coordinates in xcrysden axsf-format (isolated molecule)'
    print 'Phonons.WriteAXSFFile: Writing',filename
    f = open(filename,'w')
    f.write('ANIMSTEPS %i\n'%len(hw))    
    for i in range(len(hw)):
        f.write('ATOMS %i\n'%(i+1))
        for j in range(len(xyz)):
            ln = ' %i'%anr[j]
            for k in range(3):
                ln += ' %.6f'%xyz[j][k]
            if j < FCfirst-1 or j > FClast-1:
                ln += ' %.6f %.6f %.6f'%(0,0,0)
            else:
                for k in range(3):
                    try:
                        ln += ' %.6f'%U[i][3*(j+1-FCfirst)+k]
                    except:
                        ln += ' %.6f'%U[i,j,k]
            ln += '\n'
            f.write(ln)
    f.close()

def WriteAXSFFilesPer(filename,vectors,xyz,anr,hw,U,FCfirst,FClast):
    'Writes the vibrational normal coordinates in xcrysden axsf-format (periodic structure)'
    print 'Phonons.WriteAXSFFilePer: Writing',filename
    VEC = N.zeros((len(hw),3*len(xyz)),N.float)
    VEC[:,3*(FCfirst-1):3*FClast] = U
    f = open(filename,'w')
    f.write('ANIMSTEPS %i\nCRYSTAL\n'%len(hw))
    for i in range(len(hw)):
        f.write('PRIMVEC %i\n'%(i+1))
        f.write('%.6f %.6f %.6f\n'%(vectors[0][0],vectors[0][1],vectors[0][2]))
        f.write('%.6f %.6f %.6f\n'%(vectors[1][0],vectors[1][1],vectors[1][2]))
        f.write('%.6f %.6f %.6f\n'%(vectors[2][0],vectors[2][1],vectors[2][2]))
        f.write('PRIMCOORD %i\n'%(i+1))
        f.write('%i 1\n'%(len(xyz)))
        for j in range(len(xyz)):
            ln = ' %i'%anr[j]
            for k in range(3):
                ln += ' %.6f'%xyz[j][k]
            for k in range(3):
                ln += ' %.6f'%VEC[i][3*j+k]
            ln += '\n'
            f.write(ln)
    f.close()
    
def main(options):
    CF.CreatePipeOutput(options.DestDir+'/'+options.Logfile)
    CF.PrintMainHeader('Phonons',vinfo,options)

    # Determine SIESTA input fdf files in FCruns
    fdf = glob.glob(options.FCwildcard+'/RUN.fdf')

    print 'Phonons.Analyze: This run uses'
    FCfirst,FClast = min(options.DynamicAtoms),max(options.DynamicAtoms)
    print '  ... FCfirst     = %4i, FClast     = %4i, Dynamic atoms = %4i'\
          %(FCfirst,FClast,len(options.DynamicAtoms))
    print '  ... DeviceFirst = %4i, DeviceLast = %4i, Device atoms  = %4i'\
          %(options.DeviceFirst,options.DeviceLast,options.DeviceLast-options.DeviceFirst+1)
    print '  ... PBC First   = %4i, PBC Last   = %4i, Device atoms  = %4i'\
          %(options.PBCFirst,options.PBCLast,options.PBCLast-options.PBCFirst+1)

    # Build Dynamical Matrix
    DM = DynamicalMatrix(fdf,options.DynamicAtoms)
    DM.SetMasses(options.Isotopes)
    # Compute modes
    DM.ComputePhononModes(DM.mean)
    # Compute e-ph coupling
    DM.PrepareGradients(options.onlySdir,options.kpoint,
                        options.DeviceFirst,options.DeviceLast,options.AbsEref)
    DM.ComputeEPHcouplings(options.PBCFirst,options.PBCLast)
    # Write data to files
    DM.WriteOutput(options.DestDir+'/Output')

    CF.PrintMainFooter('Phonons')

    return DM.h0,DM.s0,DM.hw,DM.heph
