print "SVN $Id$"

"""
Routines to calculate vibrational frequencies and e-ph coupling.

"""

import os, os.path, glob, string, time, sys
import Scientific.IO.NetCDF as nc
import numpy as N
import numpy.linalg as LA
import SiestaIO as SIO
import gzip
import copy
import PhysicalConstants as PC
import WriteXMGR as XMGR
import Symmetry
import MiscMath as MM
import ValueCheck as VC

mm = MM.mm
dagger = MM.dagger

def Analyze(FCwildcard,
            onlySdir='../onlyS',
            DeviceFirst=1,DeviceLast=1e3,
            FCfirst=1,FClast=1e3,
            PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
            outlabel='Out',
            CalcCoupl=True,
            PrintSOrbitals=False,
            AbsEref=False,
            AuxNCfile=None,
            Isotopes=[],
            kpoint=[0,0,0],
            PhBandStruct="None",
            PhBandRadie=0.0,
            AuxNCUseSinglePrec=False):
    """
    Determine vibrations and electron-phonon couplings from SIESTA calculations
    Needs two types of SIESTA calculations :
    -  FC calculations, may be in several subdirectories ...
    -  onlyS calculation in subdirectory onlyS (to get derivative of
        overlap...)
    Input:
    -  FCwildcard : String specifying all FC directories to be used
    -  DeviceFirst: Restrict Heph etc to the basis orbitals of these atoms
    -  DeviceLast : - (can be a subset of the SIESTA basis)   
    -  FCfirst    : Restrict FC matrix to these atoms
    -  FClast     : - (can be a subset of the SIESTA FC calculations)
    -  PerBoundCorrFirst : Prevent interactions through periodic boundary
                    conditions (defaults to DeviceFirst)
    -  PerBoundCorrLast : Defaults to DeviceLast
    -  CalcCoupl  : Whether or not to calculate e-ph couplings
    -  PrintSOrbitals : Print e-ph couplings in the s-orbital space
    -  CorrectFermiShifts : Use instantaneous Fermi energy as reference in finite-difference scheme
    -  AuxNCfile  : An optional auxillary netcdf-file used for read/writing dH matrix arrays
                    (useful for large systems where loading dH matrices to the memory becomes
                    an issue)
    -  Isotopes   : [[ii1, anr1], ...] substitute atom number ii1 to be of type anr1 etc.,
                    e.g., to substitute hydrogen with deuterium (anr 1001).
    -  kpoint     : Electronic k-point where e-ph couplings are evaluated (Gamma default).
    -  PhBandStruct: = != None -> Calculate bandstructure. 
                      String "AUTO", "BCC", "FCC", "CUBIC", "GRAPHENE", 
                      "POLYMER" etc
    -  PhBandRadie : Optional max distance for forces. 0.0 -> Automatically choosen
    """
    
    print '=========================================================================='
    print '                         Calculating Phonons'
    print '=========================================================================='

    ### Get file lists
    print '\nPhonons.Analyze: Searching file structure.'
    tree,XVfiles,FCfirstMIN,FClastMAX = GetFileLists(FCwildcard)
    FCfirst = max(FCfirst,FCfirstMIN)
    FClast  = min(FClast,FClastMAX)

    ### Correct and check geometry in XV-files
    print '\nPhonons.Analyze: Checking XV-files:'
    corrXVfiles = []
    for file in XVfiles:
        file, displacement = CorrectXVfile(file)
        corrXVfiles.append(file)
    if len(corrXVfiles)>1:
        CheckForIdenticalXVfiles(corrXVfiles)

    ### Read geometry
    print '\nPhonons.Analyze: Reading geometry:'
    vectors,speciesnumber,atomnumber,xyz = SIO.ReadXVFile(corrXVfiles[0])

    # Determine correspondence between speciesnumber and orbital-indices
    orbitalIndices,nao = GetOrbitalIndices(tree[0][2],speciesnumber)

    # Make isotope substitutions 
    for ii,anr in Isotopes:
        print 'Phonons.Analyse: Isotope substitution for atom index %i:'%ii
        print '  ... atom type %i --> %i'%(atomnumber[ii-1],anr)
        print '  ... atom mass %.4f --> %.4f'%(PC.AtomicMass[atomnumber[ii-1]],\
                                               PC.AtomicMass[anr])
        atomnumber[ii-1] = anr
    
    DeviceFirst = max(DeviceFirst,1)
    DeviceLast = min(DeviceLast,len(xyz))
    print 'Phonons.Analyze: This run uses'
    print '  ... DeviceFirst = %4i, DeviceLast = %4i, Device atoms  = %4i'\
          %(DeviceFirst,DeviceLast,DeviceLast-DeviceFirst+1)
    print '  ... FCfirst     = %4i, FClast     = %4i, Dynamic atoms = %4i'\
          %(FCfirst,FClast,FClast-FCfirst+1)

    ### Make sure PerBoundCorr is sensible numbers
    if PerBoundCorrFirst<1:
        PerBoundCorrFirst=DeviceFirst
    if PerBoundCorrLast<1 or PerBoundCorrLast>DeviceLast:
        PerBoundCorrLast=DeviceLast
    print '  ... PBC First   = %4i, PBC Last   = %4i, Device atoms  = %4i'\
          %(PerBoundCorrFirst,PerBoundCorrLast,PerBoundCorrLast-PerBoundCorrFirst+1)

    ### Build FC-matrix
    print '\nPhonons.Analyze: Building FC matrix:'
    # First, build FC on the full [FCfirstMIN,FCfirstMAX] space
    FCm,FCp = GetFCMatrices(tree,FCfirstMIN,FClastMAX,len(xyz))
    # Correct for egg-box effect
    FCm = CorrectFCMatrix(FCm,FCfirstMIN,FClastMAX,len(xyz))
    FCp = CorrectFCMatrix(FCp,FCfirstMIN,FClastMAX,len(xyz))
    FCmean = (FCm+FCp)/2

    ### Write FC-matrix to file
    OutputFC(FCmean,filename='%s.FC'%outlabel)
    
    ### Calculate phonon bandstructure
    if PhBandStruct!="None":
        CalcBandStruct(vectors,speciesnumber,xyz,FCmean,\
                           FCfirst,FClast,PhBandStruct,atomnumber,\
                           PhBandRadie)
    
    ### Calculate phonon frequencies and modes
    print '\nPhonons.Analyze: Calculating phonons from FCmean, FCm, FCp:'
    # Mean
    FC2 = ReduceAndSymmetrizeFC(FCmean,FCfirstMIN,FClastMAX,FCfirst,FClast)
    OutputFC(FC2,filename='%s.reduced.FC'%outlabel)
    hw,U = CalcPhonons(FC2,atomnumber,FCfirst,FClast)
    # FCm
    FC2 = ReduceAndSymmetrizeFC(FCm,FCfirstMIN,FClastMAX,FCfirst,FClast)
    CalcPhonons(FC2,atomnumber,FCfirst,FClast)
    # FCp
    FC2 = ReduceAndSymmetrizeFC(FCp,FCfirstMIN,FClastMAX,FCfirst,FClast)
    CalcPhonons(FC2,atomnumber,FCfirst,FClast)

    ### Write MKL- and xyz-files
    print '\nPhonons.Analyze: Writing geometry and phonons to files.'
    SIO.WriteMKLFile('%s_FC%i-%i.mkl'%(outlabel,FCfirst,FClast),
                     atomnumber,xyz,hw,U,FCfirst,FClast)
    SIO.WriteXYZFile('%s.xyz'%outlabel,atomnumber,xyz)
    WriteFreqFile('%s.freq'%outlabel,hw)
    WriteVibDOSFile('%s.fdos'%outlabel,hw)
    WriteAXSFFiles('%s.axsf'%outlabel,xyz,atomnumber,hw,U,FCfirst, FClast)
    
    ### Write data to NC-file
    print '\nPhonons.Analyze: Writing results to netCDF-file'
    tmp1, tmp2 = [], []
    for ii in range(DeviceFirst-1,DeviceLast):
        tmp1.append(orbitalIndices[ii,0]-orbitalIndices[DeviceFirst-1,0])
        tmp2.append(orbitalIndices[ii,1]-orbitalIndices[DeviceFirst-1,0])
    naoDev = orbitalIndices[DeviceLast-1][1]-orbitalIndices[DeviceFirst-1][0]+1
    NCfile,newNCfile = OpenNetCDFFile('%s.nc'%outlabel,
                            naoDev,xyz,DeviceFirst,DeviceLast,FCfirst,FClast,AuxNCfile)
    if newNCfile:
        Write2NetCDFFile(NCfile,tmp1,'FirstOrbital',('NumDevAtoms',),
                         description='Orbital index for the first orbital on the atoms (counting from 0)')
        Write2NetCDFFile(NCfile,tmp2,'LastOrbital',('NumDevAtoms',),
                         description='Orbital index for the last orbital on the atoms (counting from 0)')
        Write2NetCDFFile(NCfile,hw,'hw',('PhononModes',),units='eV')
        Write2NetCDFFile(NCfile,U,'U',('PhononModes','PhononModes',),
                         description='U[i,j] where i is mode index and j atom displacement')
        Write2NetCDFFile(NCfile,vectors,'UnitCell',('dim3','dim3',),units='Ang')
        Write2NetCDFFile(NCfile,xyz,'GeometryXYZ',('NumTotAtoms','dim3',),units='Ang')
        Write2NetCDFFile(NCfile,atomnumber,'AtomNumbers',('NumTotAtoms',),units='Atomic Number')
        Write2NetCDFFile(NCfile,speciesnumber,'SpeciesNumbers',('NumTotAtoms',),units='Species Number')
        DeviceAtoms = range(DeviceFirst,DeviceLast+1)
        Write2NetCDFFile(NCfile,DeviceAtoms,'DeviceAtoms',('NumDevAtoms',),
                         description='Range of atomic indices (counting from 1)')
        del DeviceAtoms
        Write2NetCDFFile(NCfile,N.array(range(FCfirst,FClast+1),N.float),'DynamicAtoms',('NumFCAtoms',),
                         description='Range of atomic indices (counting from 1)')

    if CalcCoupl:
        print '\nPhonons.Analyze: AbsEref =',AbsEref
    if CalcCoupl and AuxNCfile:
        # Heph couplings utilizing netcdf-file
        if not os.path.isfile(AuxNCfile):
            GenerateAuxNETCDF(tree,FCfirst,FClast,orbitalIndices,nao,onlySdir,PerBoundCorrFirst,PerBoundCorrLast,
                              AuxNCfile,displacement,kpoint,AuxNCUseSinglePrec,AbsEref)
        else:
            print 'Phonons.Analyze: Reading from AuxNCfile =', AuxNCfile
        H0,S0,Heph = CalcHephNETCDF(orbitalIndices,FCfirst,FClast,atomnumber,DeviceFirst,DeviceLast,
                                    hw,U,NCfile,AuxNCfile,kpoint,AuxNCUseSinglePrec)
    elif CalcCoupl:
        # Old way to calculate Heph (reading everything into the memory)
        print '\nPhonons.Analyze: Reading (H0,S0,dH) from .TSHS and .onlyS files:'
        # Read electronic structure from files
        eF,H0,S0,dH = GetH0S0dH(tree,FCfirst,FClast,displacement,kpoint,AbsEref)
        # Correct dH for change in basis states with displacement
        dH = CorrectdH(onlySdir,orbitalIndices,nao,eF,H0,S0,dH,FCfirst,displacement,kpoint)
        # Downfold matrices to the subspace of the device atoms
        H0,S0,dH = Downfold2Device(orbitalIndices,H0,S0,dH,DeviceFirst,DeviceLast,
                                   FCfirst,FClast,PerBoundCorrFirst,PerBoundCorrLast)
        # Calculate e-ph couplings
        print '\nPhonons.Analyze: Calculating electron-phonon couplings:'
        Heph = CalcHeph(dH,hw,U,atomnumber,FCfirst)
        # If CalcCoupl, count the actual number of atomic orbitals
        NCfile.createDimension('NSpin',len(H0))
        print '\nPhonons.Analyze: Setting kpoint =',kpoint
        Write2NetCDFFile(NCfile,kpoint,'kpoint',('dim3',))
        # Write real part
        Write2NetCDFFile(NCfile,H0.real,'H0',('NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
        Write2NetCDFFile(NCfile,S0.real,'S0',('AtomicOrbitals','AtomicOrbitals',),units='eV')
        Write2NetCDFFile(NCfile,Heph.real,'He_ph',('PhononModes','NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
        # Write imaginary part
        GammaPoint = N.dot(kpoint,kpoint)<1e-7
        if not GammaPoint:
            Write2NetCDFFile(NCfile,H0.imag,'ImH0',('NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
            Write2NetCDFFile(NCfile,S0.imag,'ImS0',('AtomicOrbitals','AtomicOrbitals',),units='eV')
            Write2NetCDFFile(NCfile,Heph.imag,'ImHe_ph',('PhononModes','NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
        # Ensure it is complete
        NCfile.CurrentHWidx = len(hw)
        NCfile.sync()

    if CalcCoupl and PrintSOrbitals:
        # Print e-ph coupling matrices in s-orbital subspace
        for iSpin in range(len(H0)):
            print '\nPhonons.Analyze: Hamiltonian H0.real (in s-orbital subspace) Spin=',iSpin
            ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,H0[iSpin,:,:].real)
            print '\nPhonons.Analyze: Hamiltonian H0.imag (in s-orbital subspace) Spin=',iSpin
            ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,H0[iSpin,:,:].imag)
        print '\nPhonons.Analyze: Overlap matrix S0.real (in s-orbital subspace)'
        ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,S0.real)
        print '\nPhonons.Analyze: Overlap matrix S0.imag (in s-orbital subspace)'
        ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,S0.imag)
        for iSpin in range(len(H0)):
            for i in range(len(hw)):
                print '\nPhonons.Analyze: Coupling matrix Heph[%i].real (in s-orbital subspace) Spin=%i'%(i,iSpin)
                ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,
                                       DeviceFirst,DeviceLast,Heph[i,iSpin,:,:].real)
                if N.dot(kpoint,kpoint)>1e-7:
                    print '\nPhonons.Analyze: Coupling matrix Heph[%i].imag (in s-orbital subspace) Spin=%i'%(i,iSpin)
                    ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,
                                           DeviceFirst,DeviceLast,Heph[i,iSpin,:,:].imag)
    NCfile.close()

    print '=========================================================================='
    print '  Program finished:  %s '%time.ctime()
    print '=========================================================================='
    if CalcCoupl:
        return hw,Heph
    else:
        return hw,0.0

def OutputFC(FC,filename='FC.matrix'):
    print 'Phonons.OutputFC: Writing',filename
    f = open(filename,'w')
    s = N.shape(FC)
    for i in range(s[0]):
        for j in range(s[1]):
            if len(s)==2:
                f.write(string.rjust('%.4e'%FC[i,j],12))
            elif len(s)==3:
                for k in range(s[2]):
                    f.write(string.rjust('%.4e'%FC[i,j,k],12))
        f.write('\n')
    f.close()

def ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,A):
    # Print matrix A defined on the electronic subspace [FCfirst,FClast]
    # FirstOrbital refers to the element A[0,0]
    FirstOrbital = orbitalIndices[DeviceFirst-1][0]
    sIndices = []
    for i in range(FCfirst,FClast+1):
        sIndices.append(orbitalIndices[i-1][0]-FirstOrbital)
    #print 'sIndices =',sIndices
    for i in sIndices:
        print '    ',
        for j in sIndices[:10]:
            print string.rjust('%.6f '%A[i,j].real,10),
        if len(sIndices)>10: print ' ...'
        else: print


def GetOrbitalIndices(dirname,speciesnumber):
    # Determine snr (siesta number) for each label
    csl = SIO.GetFDFblock(dirname+'/RUN.fdf', KeyWord = 'ChemicalSpeciesLabel')
    csl2snr = {}
    for set in csl:
        csl2snr[set[2]] = set[0]
    # Determine nao (number of orbitals) for each snr
    ionNCfiles = glob.glob(dirname+'/*.ion.nc*')
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
    for num in speciesnumber:
        nao = snr2nao[num]
        orbitalIndices.append([tmpOrb,tmpOrb+int(nao)-1])
        tmpOrb+=nao
    return N.array(orbitalIndices),tmpOrb

         
def Downfold2Device(orbitalIndices,H0,S0,dH,DeviceFirst,DeviceLast,FCfirst,FClast,PBCFirst,PBCLast):
    # Remove Periodic Boundary terms
    PBCorbFirst = orbitalIndices[PBCFirst-1][0]
    PBCorbLast  = orbitalIndices[PBCLast-1][1]
    if FCfirst < PBCFirst:
        # we have something to remove...
        bb = (PBCFirst - FCfirst) * 3
        dH[:bb,:,PBCorbLast+1:TSHS0.nuo,:] = 0.0
        dH[:bb,:,:,PBCorbLast+1:TSHS0.nuo] = 0.0
    if PBCLast < FClast:
        aa = (PBCLast - FCfirst) * 3
        dH[aa:,:,:PBCorbFirst-1,:] = 0.0
        dH[aa:,:,:,:PBCorbFirst-1] = 0.0

    # Downfold to device subspace        
    first,last = orbitalIndices[DeviceFirst-1][0],orbitalIndices[DeviceLast-1][1]
    h0 = H0[:,first:last+1,first:last+1].copy()
    s0 = S0[first:last+1,first:last+1].copy()
    dh = []
    for i in range(len(dH)):
        dh.append(dH[i][:,first:last+1,first:last+1])
    return h0,s0,N.array(dh)


def CalcHeph(dH,hw,U,atomnumber,FCfirst):
    print 'Phonons.CalcHeph: Calculating...\n',
    const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
    Heph = N.zeros(dH.shape,dH.dtype)
    # we extend the atomic mass and the other "constants" to ease the calculation
    UcOam = N.empty((len(hw),len(hw)),N.float)
    # create transformed phonon-vector displacements
    for i in range(len(hw)):
        UcOam[:,i] = U[:,i]*const/(2*PC.AtomicMass[atomnumber[FCfirst-1+i/3]]*hw)**.5
    for i in range(len(hw)):
        # Loop over modes
        SIO.printDone(i, len(hw),'Calculating Heph')
        if hw[i]>0:
            for j in range(len(hw)):
                # Loop over atomic coordinates
                Heph[i] += UcOam[i,j]*dH[j]
        else:
            print 'Phonons.CalcHeph: Nonpositive frequency --> Zero-valued coupling matrix' 
            # already zero
        #Check that Heph is Hermitian
        for iSpin in range(len(dH[0])):
            if not N.allclose(Heph[i,iSpin,:,:],dagger(Heph[i,iSpin,:,:]),atol=1e-6):
                print 'Phonons.CalcHeph: WARNING: Coupling matrix Heph[%i,%i] not Hermitian!'%(i,iSpin)
    print '  ... Done!'
    return Heph


def CorrectdH(onlySdir,orbitalIndices,nao,eF,H0,S0,dH,FCfirst,displacement,kpoint):
    print 'Phonons.CorrectdH: Applying correction to dH...'
    onlyS0,dS = GetOnlyS(onlySdir,nao,displacement,kpoint)
    invS0H0 = N.empty((2,)+H0.shape,dtype=H0.dtype)
    invS0 = LA.inv(S0)
    for iSpin in range(len(H0)):
        invS0H0[0,iSpin,:,:] = mm(invS0,H0[iSpin,:,:])
        invS0H0[1,iSpin,:,:] = mm(H0[iSpin,:,:],invS0)
    del invS0
    # don't re-create the array every time... too expensive
    dSdij = N.zeros((nao,nao),N.complex)
    for i in range(len(dH)):
        SIO.printDone(i, len(dH),'Correcting dH')
        first,last = orbitalIndices[FCfirst-1+i/3]
        dSdij[:,first:last+1] = dS[i%3,:,first:last+1]
        for iSpin in range(len(H0)):
            dH[i,iSpin,:,:] -= mm(dagger(dSdij),invS0H0[0,iSpin,:,:]) + mm(invS0H0[1,iSpin,:,:],dSdij)
        dSdij[:,first:last+1] = 0. # reset
    del invS0H0
    print '  ... Done!'
    return dH


def GetOnlyS(onlySdir,nao,displacement,kpoint):
    print 'Phonons.GetOnlyS: Reading from', onlySdir
    onlySfiles = glob.glob(onlySdir+'/*.onlyS*')
    onlySfiles.sort()
    if len(onlySfiles)<1:
        sys.exit('Phonons.GetOnlyS: No .onlyS file found!')
    # nao is number of orbitals in the original (not-doubled) structure
    # S0 = N.zeros((nao,nao),N.float)
    # dxm,dxp = N.zeros((nao,nao),N.float),N.zeros((nao,nao),N.float)
    # dym,dyp = N.zeros((nao,nao),N.float),N.zeros((nao,nao),N.float)
    # dzm,dzp = N.zeros((nao,nao),N.float),N.zeros((nao,nao),N.float) 
    if len(onlySfiles)!=6:
        sys.exit('Phonons.GetOnlyS: Wrong number of onlyS files found!')
    else:
        for file in onlySfiles:
            #S = SIO.ReadOnlyS(file)
            thisHS = SIO.HS(file)
            thisHS.setkpoint(kpoint)
            S = thisHS.S
            del thisHS
            orbitals = len(S)/2
            if orbitals!=nao:
                print 'nao=%i,  orbitals=%i'%(nao,orbitals)
                sys.exit('Phonons.GetOnlyS: Error assigning orbitals to atoms!')    
            S0 = S[0:nao,0:nao].copy()
            dmat=S[0:nao,nao:nao*2].copy()
            if file.endswith('_1.onlyS'):
                dxm = dmat
            elif file.endswith('_2.onlyS'):
                dxp = dmat
            elif file.endswith('_3.onlyS'):
                dym = dmat
            elif file.endswith('_4.onlyS'):
                dyp = dmat
            elif file.endswith('_5.onlyS'):
                dzm = dmat
            elif file.endswith('_6.onlyS'):
                dzp = dmat
        thisd = 1e10
        for i in range(1,7):
            xyz = N.array(SIO.Getxyz(onlySdir+'/RUN_%i.fdf'%i))
            for j in range(1,len(xyz)):
                thisd = min(thisd,(N.dot(xyz[0]-xyz[j],xyz[0]-xyz[j]))**.5)
    # Check that onlyS-directory also corresponds to the same displacement
    print 'Phonons.GetOnlyS: OnlyS-displacement (min) = %.5f Ang'%thisd
    print 'Phonons.GetOnlyS: FC-displacement          = %.5f Ang'%displacement
    VC.Check("displacement-tolerance",abs(thisd-displacement)/displacement,
             "Phonons.GetOnlyS: OnlyS-displacement different from FC-displacement")
    dS = N.empty((3,)+dxp.shape,dtype=dxp.dtype)
    dS[0] = (dxp-dxm)/(2.*displacement)
    dS[1] = (dyp-dym)/(2.*displacement)
    dS[2] = (dzp-dzm)/(2.*displacement)
    return S0 , dS


def GetH0S0dH(tree,FCfirst,FClast,displacement,kpoint,AbsEref):
    dH = []
    for i in range(len(tree)):
        # Go through each subdirectory
        first,last,dir = tree[i][0],tree[i][1],tree[i][2]
        HSfiles = tree[i][4]
        if len(HSfiles)!=(last-first+1)*6+1:
            sys.exit('Phonons.GetH0S0dH: Wrong number of *.TSHS files in %s\n'+\
                     'PROGRAM STOPPED!!!'%dir)
        # The first file is the one corresponding to no displacement
        TSHS0 = SIO.HS(HSfiles[0])
        TSHS0.setkpoint(kpoint) # Here eF is moved to zero
        # Check that H0 is Hermitian
        for iSpin in range(len(TSHS0.H)):
            if not N.allclose(TSHS0.H[iSpin,:,:],dagger(TSHS0.H[iSpin,:,:]),atol=1e-6):
                print 'Phonons.GetH0S0dH: WARNING: Hamiltonian H0 not Hermitian!'
        if not N.allclose(TSHS0.S,dagger(TSHS0.S),atol=1e-6):
            print 'Phonons.GetH0S0dH: WARNING: Overlap matrix S0 not Hermitian!'
        if TSHS0.istep!=0: # the first TSHS file should have istep=0
            print "Phonons.GetH0S0dH: Assumption on file order not right ",HSfiles[0]
            raise IOError("Files not in order")
        for j in range(len(HSfiles)/2):
            if TSHS0.ia1+j/3 >= FCfirst and TSHS0.ia1+j/3 <= FClast:
                # Read TSHS file since it is within (FCfirst,FClast)
                TSHSm = SIO.HS(HSfiles[2*j+1])
                TSHSm.setkpoint(kpoint) # Here eF is moved to zero of HSm
                TSHSp = SIO.HS(HSfiles[2*j+2])
                TSHSp.setkpoint(kpoint) # Here eF is moved to zero of HSp
                if AbsEref:
                    # Measure energies wrt. TSHS0.ef (absolute energy ref).
                    for iSpin in range(len(TSHS0.H)):
                        TSHSm.H[iSpin,:,:] += (TSHSm.ef-TSHS0.ef)*TSHSm.S 
                        TSHSp.H[iSpin,:,:] += (TSHSp.ef-TSHS0.ef)*TSHSp.S
                dH.append((TSHSp.H-TSHSm.H)/(2*displacement))
                # don't wait for the garbage collector...
                del TSHSm, TSHSp
    return TSHS0.ef,TSHS0.H,TSHS0.S,N.array(dH)

    
def OpenNetCDFFile(filename,nao,xyz,DeviceFirst,DeviceLast,FCfirst,FClast,AuxNCfile):
    print 'Phonons.WriteNetCDFFile: Writing', filename
    # If no auxillary file, we *MUST* create the regular one
    i_f = False
    if AuxNCfile: i_f = os.path.isfile(filename)
    if i_f:
        # We will try and append to the file instead of reading it
        file = nc.NetCDFFile(filename,'a')
        i_f = i_f and file.dimensions['AtomicOrbitals'] == int(nao)
        i_f = i_f and file.dimensions['PhononModes'] == (FClast-FCfirst+1)*3
        i_f = i_f and file.dimensions['NumTotAtoms'] == len(xyz)
        i_f = i_f and file.dimensions['NumDevAtoms'] == DeviceLast-DeviceFirst+1
        i_f = i_f and file.dimensions['NumFCAtoms']  == FClast-FCfirst+1
        if i_f: return file,False
    try:
        # Use 64 bit format
        file = nc.NetCDFFile(filename,'wl',
                             'Created '+time.ctime(time.time()))
    except:
        file = nc.NetCDFFile(filename,'w','Created '+time.ctime(time.time()))
    file.title = 'Output from Phonons.py'
    file.createDimension('AtomicOrbitals',int(nao))
    file.createDimension('dim3',3)
    file.createDimension('dim1',1)
    file.createDimension('PhononModes',(FClast-FCfirst+1)*3)
    file.createDimension('NumTotAtoms',len(xyz))
    file.createDimension('NumDevAtoms',DeviceLast-DeviceFirst+1)
    file.createDimension('NumFCAtoms',FClast-FCfirst+1)
    file.CurrentHWidx = 0
    return file,True


def Write2NetCDFFile(file,var,varLabel,dimensions,units=None,description=None):
    print 'Phonons.Write2NetCDFFile:', varLabel, dimensions
    tmp = file.createVariable(varLabel,'d',dimensions)
    tmp[:] = var
    if units: tmp.units = units
    if description: tmp.description = description


def CalcPhonons(FC,atomnumber,FCfirst,FClast):
    NumberOfAtoms = len(atomnumber)
    DynamicAtoms = FClast-FCfirst+1
    FCtilde = N.zeros((DynamicAtoms*3,DynamicAtoms*3),N.complex)
    for i in range(DynamicAtoms*3):
        for j in range(DynamicAtoms):
            for k in range(3):
                FCtilde[i,3*j+k] = FC[i,j,k]/ \
                  ( PC.AtomicMass[atomnumber[FCfirst-1+i/3]] * PC.AtomicMass[atomnumber[FCfirst-1+j]] )**0.5
    # Solve eigenvalue problem with symmetric FCtilde
    eval,evec = LA.eigh(FCtilde)
    evec=N.transpose(evec)
    eval=N.array(eval,N.complex)
    # Calculate frequencies
    const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
    hw = const*eval**0.5 # Units in eV
    for i in range(DynamicAtoms*3):
        # Real eigenvalues are defined as positive, imaginary eigenvalues as negative
        hw[i] = hw[i].real - abs(hw[i].imag)
    hw = hw.real
    # Normalize eigenvectors
    U = evec.real.copy()
    for i in range(DynamicAtoms*3):
        U[i] = U[i]/(N.dot(U[i],U[i])**0.5)
    # Sort in order descending mode energies
    hwrev = hw[::-1] # reversed array
    Urev = U[::-1] # reversed array
    print 'Phonons.CalcPhonons: Frequencies in meV:'
    for i in range(DynamicAtoms*3):
        print string.rjust('%.3f'%(1000*hwrev[i]),9),
        if (i-5)%6==0: print
    if (i-5)%6!=0: print
    print 'Phonons.CalcPhonons: Frequencies in cm^-1:'
    for i in range(DynamicAtoms*3):
        print string.rjust('%.3f'%(hwrev[i]/PC.invcm2eV),9),
        if (i-5)%6==0: print
    if (i-5)%6!=0: print
    return hwrev,Urev


def GetFCMatrices(tree,FCfirst,FClast,NumberOfAtoms):
    'Returns FC matrices (minus,plus) in units eV/Ang^2'
    FCm = N.zeros(((FClast-FCfirst+1)*3,NumberOfAtoms,3), N.float)
    FCp = N.zeros(((FClast-FCfirst+1)*3,NumberOfAtoms,3), N.float)
    # Positive/Negative displacement FC matrices
    for i in range(len(tree)): # Go through separate FC runs
        localFCfirst,localFClast = tree[i][0],tree[i][1]
        FC = N.array(SIO.ReadFCFile(tree[i][3]))
        LocalMoves = 3*(localFClast-localFCfirst+1)
        for j in range(LocalMoves):
            thisAtom = localFCfirst+j/3
            if thisAtom >= FCfirst and thisAtom <= FClast:
                FCm[3*(localFCfirst-FCfirst)+j] = FC[2*j*NumberOfAtoms:(2*j+1)*NumberOfAtoms]
                FCp[3*(localFCfirst-FCfirst)+j] = FC[(2*j+1)*NumberOfAtoms:(2*j+2)*NumberOfAtoms]
    return FCm,FCp
    
    
def CorrectFCMatrix(FC,FCfirst,FClast,NumberOfAtoms):
    TotalMoves = 3*(FClast-FCfirst+1)
    FCcorr = FC.copy()
    for i in range(TotalMoves): # Go through each displacement calc.
        subFC = FC[i].copy()
        movedAtomIndex = FCfirst-1+i/3 # Which atom was moved?
        subFC[movedAtomIndex] = N.array([0,0,0],N.float)
        subFC[movedAtomIndex] = -1*N.sum(N.array(subFC),axis=0)
        FCcorr[i] = subFC
    print 'Phonons.CorrectFCMatrix: Done'
    return FCcorr


def ReduceAndSymmetrizeFC(FC,FCfirstMIN,FClastMAX,FCfirst,FClast):
    # Reduce FC Matrix from Device-space to FC-space
    a = 3*(FCfirst-FCfirstMIN)
    b = 3*(FClast-FCfirstMIN+1)
    FC2 = FC[a:b,FCfirst-1:FClast].copy() # Python counts from 0
    FC3 = FC[a:b,FCfirst-1:FClast].copy() # Python counts from 0
    DynamicAtoms = len(FC2[0])
    # Symmetrize square FC2 Matrix
    for i in range(DynamicAtoms):
        for j in range(3):
            for k in range(DynamicAtoms):
                for l in range(3):
                    FC3[3*i+j,k,l] = (FC2[3*i+j,k,l]+FC2[3*k+l,i,j])/2
    print 'Phonons.ReduceAndSymmetrizeFC: Done'
    return FC3


def GetFileLists(FCwildcard):
    "Returns absolute paths to (FCfiles,TSHSfiles,XVfiles) matching the FCwildcard"
    dirs,tree,FCfiles,HSfiles,XVfiles = [],[],[],[],[]
    FCfirst,FClast = 1e10,-1e10
    # Find FCwildcard directories
    for elm in glob.glob(FCwildcard):
        if os.path.isdir(elm):
            dirs.append(elm)
    dirs.sort()
    for dir in dirs:
        localFCs,localHSs,localXVs = [],[],[]
        localFCfirst,localFClast = 1e10,-1e10
        # Find FCfiles in FCwildcard directories
        FCglob = glob.glob(dir+'/*.FC*')
        FCglob.sort()
        for elm in FCglob:
            if os.path.isfile(elm):
                FCfiles.append(elm)
                localFCs.append(elm)
        if len(FCglob)!=1:
            print 'Phonons.GetFileLists: Not exactly one *.FC file in directory',dir
        # Determine FCfirst and FClast without TSHS files
        localFCfirst = int(SIO.GetFDFline(dir+'/RUN.fdf',KeyWord='MD.FCfirst')[0])
        localFClast = int(SIO.GetFDFline(dir+'/RUN.fdf',KeyWord='MD.FClast')[0])
        FCfirst = min(FCfirst,localFCfirst)
        FClast = max(FClast,localFClast)
        # Find TSHSfiles in FCwildcard directories
        HSglob = glob.glob(dir+'/*.TSHS*')
        HSglob.sort()
        for elm in HSglob:
            if elm.endswith('.TSHS.gz'):
                # We are reading gzipped *.TSHS.gz files
                firstatom = string.atoi(elm[-16:-12])
                step = string.atoi(elm[-12:-8])
            elif elm.endswith('.TSHS'):
                # We are reading unzipped *.TSHS files
                firstatom = string.atoi(elm[-13:-9])
                step = string.atoi(elm[-9:-5])
            FCfirst = min(FCfirst,firstatom)
            FClast = max(FClast,firstatom+step/6-1)
            localFCfirst = min(localFCfirst,firstatom)
            localFClast = max(localFClast,firstatom+step/6-1)
            HSfiles.append(elm)
            localHSs.append(elm)
        if (localFClast-localFCfirst+1)*6+1 != len(HSglob):
            print 'Phonons.GetFileLists: WARNING - Inconsistent number of *.TSHS files in directory %s\n'%dir
        # Find XVfiles in FCwildcard directories
        for elm in glob.glob(dir+'/*.XV*'):
            if elm.endswith('.XV') or elm.endswith('.XV.gz'):
                XVfiles.append(elm)
        # Collect local directory information
        if len(localFCs)==1:
            tree.append([localFCfirst,localFClast,dir,localFCs[0],localHSs])
        else: print 'Phonons.GetFileLists: Two FC files in a directory encountered.'
    print 'Phonons.GetFileLists: %i folder(s) match FCwildcard %s' %(len(dirs),FCwildcard)
    print 'Phonons.GetFileLists: FCfirstMIN=%i, FClastMAX=%i, DynamicAtoms=%i' \
          %(FCfirst,FClast,FClast-FCfirst+1)
    tree.sort(),FCfiles.sort(),HSfiles.sort(),XVfiles.sort()
    return tree,XVfiles,FCfirst,FClast


def CorrectXVfile(XVfile):
    dir,file = os.path.split(XVfile)
    FCfirst, FClast = 1e10,-1e10
    # Determine FCfirst and FClast without TSHS files
    localFCfirst = int(SIO.GetFDFline(dir+'/RUN.fdf',KeyWord='MD.FCfirst')[0])
    localFClast = int(SIO.GetFDFline(dir+'/RUN.fdf',KeyWord='MD.FClast')[0])
    FCfirst = min(FCfirst,localFCfirst)
    FClast = max(FClast,localFClast)
    #Old way based on the HS filenames
    #for elm in glob.glob(dir+'/*.HS'):
    #    if os.path.isfile(elm):
    #        firstatom = string.atoi(elm[-9:-6])
    #        step = string.atoi(elm[-6:-3])
    #        FCfirst = min(FCfirst,firstatom)
    #        FClast = max(FClast,firstatom+step/6-1)
    vectors,speciesnumber,atomnumber,xyz = SIO.ReadXVFile(XVfile,InUnits='Bohr',OutUnits='Ang')
    # Determine the displacement
    try:
        list = SIO.GetFDFline(dir+'/RUN.fdf','MD.FCDispl')
        d, unit = list[0], list[1]
        if unit=='Bohr' or unit=='bohr':
            displacement = float(d)*PC.Bohr2Ang
        elif unit=='Ang' or unit=='ang':
            displacement = float(d)
        print 'Phonons.CorrectXVfile: displacement %s %s found from RUN.fdf'%(d,unit)
    except:
        print 'Phonons.CorrectXVfile: displacement set to 0.04*PC.Bohr2Ang (default)'
        displacement = 0.04*PC.Bohr2Ang
    print 'Phonons.CorrectXVfile: Correcting z-coord for atom %i (by %f Ang) \n  ... in %s' \
          %(FClast,displacement,XVfile)
    xyz[FClast-1][2] -= displacement # Python counts from 0
    NewXVfile = XVfile.replace('.XV','.XV2')
    SIO.WriteXVFile(NewXVfile,vectors,speciesnumber,atomnumber,xyz,InUnits='Ang',OutUnits='Bohr')
    return NewXVfile,displacement


def CheckForIdenticalXVfiles(XVfileList):
    err = 'Phonons.CheckForIdenticalXVfiles: Error encounted in\n'
    count = 0
    for i in range(len(XVfileList)-1):
        vectors1,speciesnumber1,atomnumber1,xyz1 = SIO.ReadXVFile(XVfileList[i])
        vectors2,speciesnumber2,atomnumber2,xyz2 = SIO.ReadXVFile(XVfileList[i+1])
        if not N.allclose(vectors1,vectors2,1e-7):
            count += 1
            print err, XVfileList[i],XVfileList[i+1], '(vectors) WARNING'
        if not N.allclose(speciesnumber1,speciesnumber2,1e-7):
            count += 1
            print err, XVfileList[i],XVfileList[i+1], '(speciesnumber) WARNING'
        if not N.allclose(atomnumber1,atomnumber2,1e-7):
            count += 1
            print err, XVfileList[i],XVfileList[i+1], '(atomnumber) WARNING'
        if not N.allclose(xyz1,xyz2,1e-7):
            count += 1
            print err, XVfileList[i],XVfileList[i+1], '(xyz) WARNING'
            print '... max(abs(N.array(xyz1)-N.array(xyz2))) =',N.max(N.absolute(N.array(xyz1)-N.array(xyz2)))
    if count == 0:
        print 'Phonons.CheckForIdenticalXVfiles: XV-files in list are identical'

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

def WriteVibDOSFile(filename,hw,gam=0.001,type='Gaussian'):
    'Writes the vibrational frequencies as a Gaussian or Lorentzian broadened spectrum'
    fmin = min(hw)
    fmax = max(hw)
    erange = N.arange(fmin-40*gam,fmax+40*gam,gam/10)
    spectrum = 0.0*erange
    for i in range(len(hw)):
        if type=='Gaussian':
            spectrum += (2*N.pi)**(-.5)/gam*N.exp(N.clip(-1.0*(hw[i]-erange)**2/(2*gam**2),-300,300))
        elif type=='Lorentzian':
            spectrum += 1/N.pi*gam/((hw[i]-erange)**2+gam**2)
    # Write data to file
    print 'Phonons.WriteVibDOSFile: Writing', filename
    f = open(filename,'w')
    f.write('\n# energy/meV  VibDOS \n')
    for i in range(len(erange)):
        f.write('%.5e   %.5e\n'%(1000*erange[i],spectrum[i]))
    f.close()                                                                

def WriteAXSFFiles(filename,xyz,anr,hw,U,FCfirst,FClast):
    'Writes the vibrational normal coordinates in xcrysden axsf-format'
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


def GenerateAuxNETCDF(tree,FCfirst,FClast,orbitalIndices,nao,onlySdir,PBCFirst,PBCLast,
                      AuxNCfile,displacement,kpoint,SinglePrec,AbsEref):
    if SinglePrec: precision = 'f'
    else: precision = 'd'
    SIO.HS.setkpoint2=SIO.HS.setkpoint
    def mysetk(self,kpt):
        if N.dot(kpt,kpt)<1e-7: 
            if SinglePrec: type = N.float32
            else: type = N.float
        else:
            if SinglePrec: type = N.complex64
            else: type = N.complex
        self.setkpoint2(kpt)
        try: self.H = N.array(self.H,type)
        except: pass
        self.S = N.array(self.S,type)
    SIO.HS.setkpoint=mysetk
    GammaPoint = N.dot(kpoint,kpoint)<1e-7

    index = 0
    # Read TSHSfiles
    for i in range(len(tree)):
        # Go through each subdirectory
        first,last,dir = tree[i][0],tree[i][1],tree[i][2]
        HSfiles = tree[i][4]
        if len(HSfiles)!=(last-first+1)*6+1:
            sys.exit('Phonons.GenerateAuxNETCDF: Wrong number of *.TSHS files in %s\n'%dir+\
                     ' PROGRAM STOPPED!!!')
        # The first TSHSfile in a folder is the one corresponding to no displacement
        TSHS0 = SIO.HS(HSfiles[0])
        TSHS0.setkpoint(kpoint) # Here eF is moved to zero
        if TSHS0.istep!=0: # the first TSHS file should have istep=0
            print "Phonons.GenerateAuxNETCDF: Assumption on file order not right ",HSfiles[0]
            raise IOError("File order not good")
        for j in range(len(HSfiles)/2):
            if TSHS0.ia1+j/3 >= FCfirst and TSHS0.ia1+j/3 <= FClast:
                # Read TSHS file since it is within (FCfirst,FClast)
                TSHSm = SIO.HS(HSfiles[2*j+1])
                TSHSm.setkpoint(kpoint) # Here eF is moved to zero of HSm
                TSHSp = SIO.HS(HSfiles[2*j+2])
                TSHSp.setkpoint(kpoint) # Here eF is moved to zero of HSp
                if AbsEref:
                    # Measure energies wrt. TSHS0.ef (absolute energy ref).
                    for iSpin in range(len(TSHS0.H)):
                        TSHSm.H[iSpin,:,:] += (TSHSm.ef-TSHS0.ef)*TSHSm.S
                        TSHSp.H[iSpin,:,:] += (TSHSp.ef-TSHS0.ef)*TSHSp.S
                # Calculate differences
                tmpdH = (TSHSp.H-TSHSm.H)/(2*displacement)
                del TSHSm, TSHSp
                try:
                    RedH[index,:] = tmpdH.real
                    if not GammaPoint: ImdH[index,:] = tmpdH.imag
                    index += 1
                except:
                    # First time a netcdf-file is created
                    print 'Phonons.GenerateAuxNETCDF: Generating', AuxNCfile
                    try:
                        # Use 64 bit format
                        NCfile2 = nc.NetCDFFile(AuxNCfile,'wl','Created '+time.ctime(time.time()))
                    except:
                        NCfile2 = nc.NetCDFFile(AuxNCfile,'w','Created '+time.ctime(time.time()))
                    NCfile2.createDimension('Index',None)
                    NCfile2.createDimension('NSpin',TSHS0.nspin)
                    NCfile2.createDimension('AtomicOrbitals',TSHS0.nuo)
                    NCfile2.createDimension('dim1',1)
                    NCfile2.createDimension('dim3',3)
                    FCfirsttmp = NCfile2.createVariable('FCfirst','d',('dim1',))
                    FCfirsttmp[:] = FCfirst
                    FClasttmp = NCfile2.createVariable('FClast','d',('dim1',))
                    FClasttmp[:] = FClast
                    kpointtmp = NCfile2.createVariable('kpoint','d',('dim3',))
                    kpointtmp[:] = kpoint
                    # Write real part of matrices
                    ReH = NCfile2.createVariable('H0',precision,('NSpin','AtomicOrbitals','AtomicOrbitals',))
                    print 'Phonons.GenerateAuxNETCDF: Precision=', precision
                    print 'Phonons.GenerateAuxNETCDF: TSHS0.H.dtype=', TSHS0.H.dtype
                    ReH[:] = TSHS0.H.real
                    ReS = NCfile2.createVariable('S0',precision,('AtomicOrbitals','AtomicOrbitals',))
                    ReS[:] = TSHS0.S.real
                    RedH = NCfile2.createVariable('dH',precision,('Index','NSpin','AtomicOrbitals','AtomicOrbitals',))
                    RedH[index,:] = tmpdH.real
                    # Write imag part of matrices
                    if not GammaPoint:
                        ImH = NCfile2.createVariable('ImH0',precision,('NSpin','AtomicOrbitals','AtomicOrbitals',))
                        ImH[:] = TSHS0.H.imag
                        ImS = NCfile2.createVariable('ImS0',precision,('AtomicOrbitals','AtomicOrbitals',))
                        ImS[:] = TSHS0.S.imag
                        ImdH = NCfile2.createVariable('ImdH',precision,('Index','NSpin','AtomicOrbitals','AtomicOrbitals',))
                        ImdH[index,:] = tmpdH.imag
                    index += 1
    NCfile2.sync()
    print 'Phonons.GenerateAuxNETCDF: len(dH) =',index
    
    # Correct dH
    onlyS0,dS = GetOnlyS(onlySdir,nao,displacement,kpoint)
    invS0H0 = N.empty((2,)+TSHS0.H.shape,dtype=TSHS0.H.dtype)
    invS0 = LA.inv(TSHS0.S)
    for iSpin in range(TSHS0.nspin):
        invS0H0[0,iSpin,:,:] = mm(invS0,TSHS0.H[iSpin,:,:])
        invS0H0[1,iSpin,:,:] = mm(TSHS0.H[iSpin,:,:],invS0)
    del invS0
    # don't re-create the array every time... too expensive
    if SinglePrec:
        dSdij = N.zeros((nao,nao),N.complex64)
    else:
        dSdij = N.zeros((nao,nao),N.complex)
    for i in range(index):
        SIO.printDone(i,index,'Correcting dH')
        first,last = orbitalIndices[FCfirst-1+i/3]
        dSdij[:,first:last+1] = dS[i%3,:,first:last+1]
        for iSpin in range(TSHS0.nspin):
            tmp = mm(dagger(dSdij),invS0H0[0,iSpin,:,:]) + mm(invS0H0[1,iSpin,:,:],dSdij)
            RedH[i,iSpin,:] -= tmp.real
            if not GammaPoint: 
                ImdH[i,iSpin,:] -= tmp.imag
        dSdij[:,first:last+1] = 0. # reset
    del tmp
    print '  ... Done!'
    NCfile2.sync()
    
    # Remove Periodic Boundary terms
    PBCorbFirst = orbitalIndices[PBCFirst-1][0]
    PBCorbLast  = orbitalIndices[PBCLast-1][1]
    if FCfirst < PBCFirst:
        # we have something to remove...
        bb = (PBCFirst - FCfirst) * 3
        RedH[:bb,:,PBCorbLast+1:TSHS0.nuo,:] = 0.
        RedH[:bb,:,:,PBCorbLast+1:TSHS0.nuo] = 0.
        if not GammaPoint:
            ImdH[aa:bb,:,PBCorbLast+1:TSHS0.nuo,:] = 0.
            ImdH[aa:bb,:,:,PBCorbLast+1:TSHS0.nuo] = 0.
    if PBCLast < FClast:
        aa = (PBCLast - FCfirst) * 3
        RedH[aa:,:,:PBCorbFirst-1,:] = 0.
        RedH[aa:,:,:,:PBCorbFirst-1] = 0.
        if not GammaPoint:
            ImdH[aa:,:,:PBCorbFirst-1,:] = 0.
            ImdH[aa:,:,:,:PBCorbFirst-1] = 0.

    NCfile2.close()
    print 'Phonons.GenerateAuxNETCDF: File %s written.'%AuxNCfile


def CalcHephNETCDF(orbitalIndices,FCfirst,FClast,atomnumber,DeviceFirst,DeviceLast,
                   hw,U,NCfile,AuxNCfile,kpoint,SinglePrec):
    if SinglePrec:
        precision = 'f'
    else:
        precision = 'd'
    GammaPoint = N.dot(kpoint,kpoint)<1e-7
    # Read AuxNCfile and downfold to device
    first,last = orbitalIndices[DeviceFirst-1][0],orbitalIndices[DeviceLast-1][1]
    NCfile2 = nc.NetCDFFile(AuxNCfile,'r')
    ReH0 = N.array(NCfile2.variables['H0'])[:,first:last+1,first:last+1]
    ReS0 = N.array(NCfile2.variables['S0'])[first:last+1,first:last+1]
    RedH = NCfile2.variables['dH']
    if not GammaPoint:
        ImH0 = N.array(NCfile2.variables['ImH0'])[:,first:last+1,first:last+1]
        ImS0 = N.array(NCfile2.variables['ImS0'])[first:last+1,first:last+1]
        ImdH = NCfile2.variables['ImdH']
    auxkpoint = NCfile2.variables['kpoint'][:]
    print '   ... kpoint from %s ='%AuxNCfile,auxkpoint
    VC.Check("same-kpoint",N.abs(auxkpoint-kpoint),
             "Aux. file does not match specified k-point.",
             "Delete {0} and try again!".format(AuxNCfile))
    newNCfile = NCfile.CurrentHWidx == 0
    if newNCfile:
        NCfile.createDimension('NSpin',len(ReH0))
        Write2NetCDFFile(NCfile,ReH0,'H0',
                         ('NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
        Write2NetCDFFile(NCfile,ReS0,'S0',
                         ('AtomicOrbitals','AtomicOrbitals',),units='eV')
        if not GammaPoint:
            Write2NetCDFFile(NCfile,ImH0,'ImH0',
                             ('NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
            Write2NetCDFFile(NCfile,ImS0,'ImS0',
                             ('AtomicOrbitals','AtomicOrbitals',),units='eV')
            Write2NetCDFFile(NCfile,kpoint,'kpoint',('dim3',))
    else:
        if NCfile.dimensions['NSpin'] != len(ReH0):
            sys.exit('Phonons.CalcHephNETCDF: Current %s does not correspond to correct spin-component, please delete file'%(NCfile))


    # Check AuxNCfile
    if FCfirst != int(NCfile2.variables['FCfirst'][0]):
        sys.exit('Phonons.CalcHephNETCDF: AuxNCfile %s does not correspond to FCfirst = %i'%(AuxNCfile,FCfirst))
    if FClast != int(NCfile2.variables['FClast'][0]):
        sys.exit('Phonons.CalcHephNETCDF: AuxNCfile %s does not correspond to FClast = %i'%(AuxNCfile,FClast))
        
    # CalcHeph
    print 'Phonons.CalcHephNETCDF: Calculating...\n',
    const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
    if newNCfile:
        ReHeph  = NCfile.createVariable('He_ph',precision,('PhononModes','NSpin','AtomicOrbitals','AtomicOrbitals',))
        if not GammaPoint: 
            ImHeph  = NCfile.createVariable('ImHe_ph',precision,('PhononModes','NSpin','AtomicOrbitals','AtomicOrbitals',))
        NCfile.sync()
    else:
        # Retrieve the data points as we are appending
        ReHeph = NCfile.variables['He_ph']
        if not GammaPoint: 
            ImHeph = NCfile.variables['ImHe_ph']

    # we extend the atomic mass and the other "constants" to ease the calculation
    UcOam = N.empty((len(hw),len(hw)),N.float)
    # create transformed phonon-vector displacements
    for i in range(len(hw)):
        UcOam[:,i] = U[:,i]*const/(2*PC.AtomicMass[atomnumber[FCfirst-1+i/3]]*hw)**.5

    for i in range(len(hw)):
        # skip already calculated contributions
        if i < NCfile.CurrentHWidx: continue
        # Loop over modes
        SIO.printDone(i, len(hw),'Calculating Heph')
        if hw[i]>0:
            # initialize
            ReHeph[i,:] = UcOam[i,0] * RedH[0][:,first:last+1,first:last+1]
            if not GammaPoint:
                ImHeph[i,:] = UcOam[i,0] * ImdH[0][:,first:last+1,first:last+1]
            for j in range(1,len(hw)):
                # Loop over atomic coordinates
                ReHeph[i,:] += UcOam[i,j] * RedH[j][:,first:last+1,first:last+1]
                if not GammaPoint:
                    ImHeph[i,:] += UcOam[i,j] * ImdH[j][:,first:last+1,first:last+1]

        else:
            ReHeph[i,:] = 0.
            if not GammaPoint: ImHeph[i,:] = 0.
            print 'Phonons.CalcHephNETCDF: Nonpositive frequency --> Zero-valued coupling matrix' 

        #Check that Heph is Hermitian
        for iSpin in range(len(ReH0)):
            if not GammaPoint:
                if not N.allclose(ReHeph[i,iSpin]+1j*ImHeph[i,iSpin],
                                  N.transpose(ReHeph[i,iSpin]-1j*ImHeph[i,iSpin]),
                                  atol=1e-5):
                    print 'Phonons.CalcHephNETCDF: WARNING: Coupling matrix Heph[%i,%i] not Hermitian!'%(i,iSpin)
            else:
                if not N.allclose(ReHeph[i,iSpin],N.transpose(ReHeph[i,iSpin]),atol=1e-5):
                    print 'Phonons.CalcHephNETCDF: WARNING: Coupling matrix Heph[%i,%i] not Hermitian!'%(i,iSpin)
        # Update ref-counter
        NCfile.CurrentHWidx = i + 1
        NCfile.sync()
    print '  ... Done!'
    #NCfile2.close() # Leave open for later access...
    return ReH0,ReS0,ReHeph

########### Phonon bandstructure ##########################################
def CalcBandStruct(vectors,speciesnumber,xyz,FCmean,FCfirst,FClast,\
                       LatticeType,atomnumber,PhBandRadie):
    print """
    PhBandStruct:
    """    
    # Symmetrize forces
    Sym = Symmetry.Symmetry()
    Sym.setupGeom(vectors, speciesnumber, atomnumber, xyz)
    FCsym = Sym.symmetrizeFC(FCmean,FCfirst,FClast,radi=PhBandRadie)
    FCsym = CorrectFCMatrix(FCsym,FCfirst,FClast,Sym.NN)

    if Sym.basis.NN!=FClast-FCfirst+1:
        raise ValueError("Phonons: ERROR: FCfirst/last do not fit with the number of atoms in the basis (%i)."%Sym.basis.NN)
    
    # a_i : real space lattice vectors
    a1,a2,a3 = Sym.a1,Sym.a2,Sym.a3
    numBasis =  Sym.basis.NN
    
    # Use automatically determined lattice type or override
    if LatticeType != "AUTO":
        Sym.latticeType = LatticeType
        print "Symmetry: Using lattice \"%s\" for phonon bands."%LatticeType

    # Calculate lattice vectors for phase factors
    # The closest cell might be different depending on which atom is moved
    sxyz = Sym.xyz.copy()
    xyz = N.zeros((FClast-FCfirst+1,Sym.NN,3),N.float)
    for i0 in range(FClast-FCfirst+1):
        ii = FCfirst-1+i0
        micxyz = Symmetry.moveIntoClosest(\
                sxyz-sxyz[ii,:],Sym.pbc[0],Sym.pbc[1],Sym.pbc[2])
        for jj in range(Sym.NN):
            xyz[i0,jj,:] = micxyz[jj,:]+sxyz[ii,:]-\
                sxyz[FCfirst-1+Sym.basisatom[jj],:]

    # Mass scaled dynamical matrix
    MSFC = FCsym.copy()
    for ii in range(len(MSFC[:,0])):
        for jj in range(len(MSFC[0,:])):
            MSFC[ii,jj,:]=MSFC[ii,jj,:]/N.sqrt(\
                PC.AtomicMass[atomnumber[FCfirst-1+ii/3]]*\
                    PC.AtomicMass[atomnumber[jj]])

    # Get bandlines
    what = Sym.what()

    bands, datasets, kdist = [], [], []
    for elem in what:
        kdist += [N.dot(elem[1],elem[1])**.5]
        # Calculate bands
        bands+=[ calcPhBands(MSFC, a1, a2, a3, xyz,\
                                 elem[1], elem[2], elem[3], numBasis,\
                                 Sym.basisatom)]
        f=open(elem[0]+'.dat','w')
        for ii,data in enumerate(bands[-1]):
            f.write("%i "%ii)
            for jj in data:
                f.write("%e "%(jj.real))
            f.write("\n")
        f.close()
        xx = N.array(range(elem[3]),N.float)/(elem[3]-1.0)
        #datasets += [[XMGR.XYDYset(xx,bands[-1][:,ii].real,bands[-1][:,ii].imag/2,Lcolor=ii+1) for ii in range(len(bands[-1][0,:]))]]
        datasets += [[XMGR.XYset(xx,bands[-1][:,ii].real,Lcolor=1,Lwidth=2) for ii in range(len(bands[-1][0,:]))]]

    maxEnergy = N.max(bands).real
    munit = 50.0
    maxEnergy = int(1.05*maxEnergy/munit)*munit+munit

    # Compute x-axis widths in the final plot
    widths = N.array(kdist)/N.sum(kdist)

    gg=[]
    xmin = 0.15
    for jj, ii in enumerate(datasets):
        g =XMGR.Graph()
        for data in ii:
            g.AddDatasets(data)
        g.SetSubtitle(what[jj][0])
        g.SetXaxis(label='',majorUnit=1.0,minorUnit=1.0,max=1,min=0,useticklabels=False)     
        if jj==0:
            g.SetYaxis(label='Phonon energy (meV)',majorUnit=munit,minorUnit=munit/5,max=maxEnergy,min=0)
        else:
            g.SetYaxis(label='',majorUnit=munit,minorUnit=munit/5,max=maxEnergy,min=0,useticklabels=False)
        # Place graph on canvas
        xmax = xmin + widths[jj]
        g.SetView(xmin=xmin,xmax=xmax,ymin=0.15,ymax=0.85)    
        xmin = 1*xmax
        gg+=[g]

    p = XMGR.Plot('PhononBands_%.2f.agr'%PhBandRadie,gg[0])
    
    for ii in range(1,len(gg)):
        p.AddGraphs(gg[ii])
    #p.ArrangeGraphs(nx=len(gg),ny=1,hspace=0.0,vspace=0.0)

    # Finally, write the plot file
    #p.ShowTimestamp()
    p.WriteFile()
    p.Print2File('PhononBands_%.2f.eps'%PhBandRadie)


def calcPhBands(FCmean, a1, a2, a3, xyz, kdir, korig, Nk, Nbasis, basisatom):
    NN = 3*Nbasis
    Band = N.zeros((Nk, NN),N.complex)
    units = PC.eV2Joule*1e20/PC.amu2kg*(PC.hbar2SI/PC.eV2Joule)**2*1e6
    for ik in range(Nk):
        kvec = kdir*ik/float(Nk-1)+korig
        #k1, k2, k3 = N.dot(a1,kvec), N.dot(a2,kvec), N.dot(a3,kvec)
        FCk = N.zeros((NN,NN),N.complex)
        for ii in range(Nbasis):
            for jj in range(len(xyz[0,:])):
                tmp=FCmean[ii*3:(ii+1)*3,jj,:]*N.exp(2.0j*N.pi*\
                       (N.dot(kvec,xyz[ii,jj,:])))
                FCk[ii*3:(ii+1)*3,(basisatom[jj])*3:(basisatom[jj]+1)*3]+=N.transpose(tmp)
        
        #print N.max(N.abs((FCk-N.conjugate(N.transpose(FCk)))))/N.max(N.abs(FCk))
        FCk = (FCk+dagger(FCk))/2
        a,b = LA.eig(FCk*units)
        Band[ik,:] = N.sort(N.sqrt(a))
    return Band

