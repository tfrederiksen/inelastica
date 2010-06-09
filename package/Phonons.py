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


def Analyze(dirname,wildcard,
            onlySdir='../onlyS',
            newPHrun='PhononCalc',
            DeviceFirst=1,DeviceLast=1e3,
            FCfirst=1,FClast=1e3,
            BulkAtomsLeft=-1,BulkAtomsRight=-1,
            output2file=False,outlabel='Out',
            CorrPotentialShift=True,
            TryMirrorSymmetry=False,
            CalcCoupl=True,
            PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
            PrintSOrbitals=True,
            AuxNCfile=None,Isotopes=[]):
    """
    Calculate electron phonon coupling from siesta calculations
    Needs two types of siesta calculations :
    -  FC calculation, may be in several subdirectories ...
    -  onlyS calculation in subdirectory onlyS (to get derivative of
        overlap...)
    Input:
    -  dirname    : main directory
    -  wildcard   : e.g. 'Au_*', gives all subdirs for FC calculations
    -  Ef         : Fermi energy used in correcting the Heph in the "old way"
    -  DeviceFirst: Restrict Heph etc to the basis orbitals of these
    -  DeviceLast :   atoms
    -  FCfirst    : restrict FC matrix to these atoms
    -  FClast     : (may be different from the siesta FC calculation)
    -  Which_eph  : 0 (new corrected), 1 (old corrected, needs Ef),
                    2 (uncorrected)
    -  PerBoundCorrFirst : Prevent interactions through periodic
                           boundary conditions defaults to DeviceFirst
    -  PerBoundCorrLast  : Defaults to DeviceLast
    -  AuxNCfile  : An optional auxillary netcdf-file used for read/writing dH matrix arrays
    -  Isotopes   : [[ii1, anr1], ...] substitute to atom type anr for siesta numbering atom ii
                    I.e., use to substitute deuterium (anr 1001) for hydrogens
    """
    
    ### Make directory for output files etc.
    phononDirectory = dirname+'/'+newPHrun
    if not os.path.isdir(phononDirectory):
        os.mkdir(phononDirectory)

    if output2file:
        # Redirecting outputs
        outfile = phononDirectory+'/'+outlabel+'.log'
        sys.stdout = open(outfile,'w')
        sys.stderr = open(outfile,'w')
    
    print '=========================================================================='
    print '                         Calculating Phonons'
    print '=========================================================================='

    ### Get file lists
    print '\nPhonons.Analyze: Searching file structure.'
    tree,XVfiles,FCfirstMIN,FClastMAX = GetFileLists(dirname,wildcard)
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
        atomnumber[ii-1]=anr
        print "Phonons.Analyse: Isotope substitution for atom index %i to atom type %i"%(ii,anr)
        
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

    ### Check for mirror symmetries
    if TryMirrorSymmetry:
        print '\nPhonons.Analyze: Checking for mirror symmetry in [FCfirstMIN,FCfirstMAX]:'
        pairs = FindMirrorSymmetry(FCfirstMIN,FClastMAX,xyz)
    else:
        pairs = {}

    ### Build FC-matrix
    print '\nPhonons.Analyze: Building FC matrix:'
    # First, build FC on the full [FCfirstMIN,FCfirstMAX] space
    FCm,FCp = GetFCMatrices(tree,FCfirstMIN,FClastMAX,len(xyz))
    # Correct for egg-box effect
    FCm = CorrectFCMatrix(FCm,FCfirstMIN,FClastMAX,len(xyz))
    FCp = CorrectFCMatrix(FCp,FCfirstMIN,FClastMAX,len(xyz))
    if len(pairs)>0:
        # Apply mirror symmetry correction
        FCm = MirrorCorrection(FCm,FCfirstMIN,FClastMAX,pairs)
        FCp = MirrorCorrection(FCp,FCfirstMIN,FClastMAX,pairs)
    FCmean = (FCm+FCp)/2

    ### Write mass-scaled FC-matrix to file
    OutputFC(FCmean,filename=phononDirectory+'/%s.MSFC'%outlabel)
    
    ### Calculate phonon frequencies and modes
    print '\nPhonons.Analyze: Calculating phonons from FCmean, FCm, FCp:'
    hw,U = [],[]    
    for FC in [FCmean,FCm,FCp]:
        FC2 = ReduceAndSymmetrizeFC(FC,FCfirstMIN,FClastMAX,FCfirst,FClast)
        a,b = CalcPhonons(FC2,atomnumber,FCfirst,FClast)
        hw.append(a),U.append(b)
    hw,U = hw[0],U[0] # Keep only those from FCmean

    ### Calculate phonon spectrum
    a,b,c = N.shape(FCmean)
    if a==b*c and BulkAtomsLeft>0 and BulkAtomsRight>0:
        # Force constants available for whole unit cell
        print '\nPhonons.Analyze: Calculating phonon spectrum.'
        CalcPhononDOS(FCmean,U,atomnumber,FCfirst,FClast,BulkAtomsLeft,BulkAtomsRight)

    ### Write MKL- and xyz-files
    print '\nPhonons.Analyze: Writing geometry and phonons to files.'
    SIO.WriteMKLFile(phononDirectory+'/%s_FC%i-%i.mkl'%(outlabel,FCfirst,FClast),
                     atomnumber,xyz,hw,U,FCfirst,FClast)
    SIO.WriteXYZFile(phononDirectory+'/%s.xyz'%outlabel,atomnumber,xyz)
    WriteFreqFile(phononDirectory+'/%s.freq'%outlabel,hw)
    WriteVibDOSFile(phononDirectory+'/%s.fdos'%outlabel,hw)
    WriteAXSFFiles(phononDirectory+'/%s.axsf'%outlabel,xyz,atomnumber,hw,U,FCfirst, FClast)
    
    ### Write data to NC-file
    print '\nPhonons.Analyze: Writing results to netCDF-file'
    tmp1, tmp2 = [], []
    for ii in range(DeviceFirst-1,DeviceLast):
        tmp1.append(orbitalIndices[ii,0]-orbitalIndices[DeviceFirst-1,0])
        tmp2.append(orbitalIndices[ii,1]-orbitalIndices[DeviceFirst-1,0])
    naoDev = orbitalIndices[DeviceLast-1][1]-orbitalIndices[DeviceFirst-1][0]+1
    NCfile = OpenNetCDFFile(phononDirectory+'/%s.nc'%outlabel,
                            naoDev,xyz,DeviceFirst,DeviceLast,FCfirst,FClast)
    Write2NetCDFFile(NCfile,tmp1,'FirstOrbital',('NumDevAtoms',),
                     description='Orbital index for the first orbital on the atoms (counting from 0)')
    Write2NetCDFFile(NCfile,tmp2,'LastOrbital',('NumDevAtoms',),
                     description='Orbital index for the last orbital on the atoms (counting from 0)')
    Write2NetCDFFile(NCfile,hw,'hw',('PhononModes',),units='eV')
    Write2NetCDFFile(NCfile,U,'U',('PhononModes','PhononModes',),
                     description='U[i,j] where i is mode index and j atom displacement')
    Write2NetCDFFile(NCfile,vectors,'UnitCell',('XYZ','XYZ',),units='Ang')
    Write2NetCDFFile(NCfile,xyz,'GeometryXYZ',('NumTotAtoms','XYZ',),units='Ang')
    Write2NetCDFFile(NCfile,atomnumber,'AtomNumbers',('NumTotAtoms',),units='Atomic Number')
    Write2NetCDFFile(NCfile,speciesnumber,'SpeciesNumbers',('NumTotAtoms',),units='Species Number')
    DeviceAtoms = range(DeviceFirst,DeviceLast+1)
    Write2NetCDFFile(NCfile,DeviceAtoms,'DeviceAtoms',('NumDevAtoms',),
                     description='Range of atomic indices (counting from 1)')
    Write2NetCDFFile(NCfile,range(FCfirst,FClast+1),'DynamicAtoms',('NumFCAtoms',),
                     description='Range of atomic indices (counting from 1)')

    if CalcCoupl and AuxNCfile:
        # Heph couplings utilizing netcdf-file
        if not os.path.isfile(AuxNCfile):
            GenerateAuxNETCDF(tree,FCfirst,FClast,orbitalIndices,nao,onlySdir,PerBoundCorrFirst,PerBoundCorrLast,
                              AuxNCfile,displacement,CorrPotentialShift=CorrPotentialShift)
        else:
            print 'Phonons.Analyze: Reading from AuxNCfile =', AuxNCfile
        H0,S0,Heph = CalcHephNETCDF(orbitalIndices,FCfirst,FClast,atomnumber,DeviceFirst,DeviceLast,
                                  hw,U,NCfile,AuxNCfile)
    elif CalcCoupl:
        # Old way to calculate Heph (reading everything into the memory)
        print '\nPhonons.Analyze: Reading (H0,S0,dH) from .TSHS and .onlyS files:'
        # Read electronic structure from files
        eF,H0,S0,dH = GetH0S0dH(tree,FCfirst,FClast,displacement,CorrPotentialShift=CorrPotentialShift)
        # Correct dH for change in basis states with displacement
        dH = CorrectdH(onlySdir,orbitalIndices,nao,eF,H0,S0,dH,FCfirst,displacement)
        # Downfold matrices to the subspace of the device atoms
        H0,S0,dH = Downfold2Device(orbitalIndices,H0,S0,dH,DeviceFirst,DeviceLast,
                                   FCfirst,FClast,PerBoundCorrFirst,PerBoundCorrLast)
        # Calculate e-ph couplings
        print '\nPhonons.Analyze: Calculating electron-phonon couplings:'
        Heph = CalcHeph(dH,hw,U,atomnumber,FCfirst)
        # If CalcCoupl, count the actual number of atomic orbitals
        NCfile.createDimension('NSpin',len(H0))
        Write2NetCDFFile(NCfile,H0,'H0',('NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
        Write2NetCDFFile(NCfile,S0,'S0',('AtomicOrbitals','AtomicOrbitals',),units='eV')
        Write2NetCDFFile(NCfile,Heph,'He_ph',('PhononModes','NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')

    if CalcCoupl and PrintSOrbitals:
        # Print e-ph coupling matrices in s-orbital subspace
        for iSpin in range(len(H0)):
            print '\nPhonons.Analyze: Hamiltonian H0 (in s-orbital subspace) Spin=',iSpin
            ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,H0[iSpin,:,:])
        print '\nPhonons.Analyze: Overlap matrix S0 (in s-orbital subspace)'
        ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,DeviceFirst,DeviceLast,S0)
        for iSpin in range(len(H0)):
            for i in range(len(hw)):
                print '\nPhonons.Analyze: Coupling matrix Heph[%i] (in s-orbital subspace) Spin=%i'%(i,iSpin)
                ShowInSOrbitalSubspace(orbitalIndices,FCfirst,FClast,
                                       DeviceFirst,DeviceLast,Heph[i,iSpin,:,:])
    NCfile.close()

    print '=========================================================================='
    print '  Program finished:  %s '%time.ctime()
    print '=========================================================================='
    

def OutputFC(FC,filename='FC.matrix'):
    print 'Phonons.OutputFC: Writing',filename
    f = open(filename,'w')
    s = N.shape(FC)
    for i in range(s[0]):
        for j in range(s[1]):
            if len(s)==2:
                f.write(string.rjust('%.2e'%FC[i,j],10))
            elif len(s)==3:
                for k in range(s[2]):
                    f.write(string.rjust('%.2e'%FC[i,j,k],10))
        f.write('\n')
    f.close()

def MirrorCorrection(FC,FCfirstMIN,FClastMAX,pairs):
    FC2 = 0.0*FC.copy()
    T = N.identity(3,N.float)
    T[2,2] = -1.0
    mm = N.dot
    for i in pairs:
        for j in pairs:
            #print i,j,' (%i,%i)'%(pairs[i],pairs[j])
            a,Ta = 3*(i-FCfirstMIN),3*(pairs[i]-FCfirstMIN)
            fc1 = FC[a:a+3,j-1]
            fc2 = mm(mm(T,FC[Ta:Ta+3,pairs[j]-1]),T)
            fcm = (fc1+fc2)/2.0
            FC2[a:a+3,j-1] = fcm
    print 'Phonons.MirrorCorrection: Applied to FC'
    return FC2


def FindMirrorSymmetry(First,Last,xyz):
    MirrorTol=0.1 # Ang
    Zmirror = 0.0
    for i in range(First,Last+1):
        Zmirror += xyz[i-1][2] # Add up z-coord among dynamic atoms
    Zmirror = Zmirror/(Last-First+1) # Calculate mean z-coord
    # Go through device atoms and find pairs
    pairs = {}
    for i in range(First,Last+1):
        for j in range(First,Last+1):
            dev = [abs(xyz[i-1][0]-xyz[j-1][0]),
                   abs(xyz[i-1][1]-xyz[j-1][1]),
                   abs(xyz[i-1][2]+xyz[j-1][2]-2*Zmirror)]
            if max(dev) < MirrorTol:
                # (i,j) forms a mirror pair
                print 'Phonons.FindMirrorSymmetry: Atoms (%i,%i) mirrors around z=%f Ang'%(i,j,Zmirror)
                print '... abs([x1-x2,y1-y2,z1+z2-2*Zmirror]) = [%.3f,%.3f,%.3f] Ang'%(dev[0],dev[1],dev[2])
                pairs[i] = j # Add pair to dictionary
    if len(pairs)==0:
        print 'Phonons.FindMirrorSymmetry: No mirror symmetries was found (tol=%f Ang).'%tol
    return pairs


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
            print string.rjust('%.6f '%A[i,j],10),
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
    print 'snr2nao =',snr2nao
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
    PBCorbFirst=orbitalIndices[PBCFirst-1][0]
    PBCorbLast=orbitalIndices[PBCLast-1][1]
    
    for ii in range(FCfirst-1,FClast):
        if ii+1<PBCFirst:
            # Left contact, remove contribution to right contact
            print "Phonons.Downfold2Device: Removing He-ph, Left atom %i" % (ii+1)
            for aa in range(PBCorbLast+1,len(S0)):
                for bb in range(len(S0)):
                    for dir in range(3): # x,y,z
                        for iSpin in range(len(H0)):
                            dH[(ii-FCfirst+1)*3+dir,iSpin,aa,bb]=0.0
                            dH[(ii-FCfirst+1)*3+dir,iSpin,bb,aa]=0.0
        if ii+1>PBCLast:
            # Right contact, remove contribution to left contact
            print "Phonons.Downfold2Device: Removing He-ph, Right atom %i" % (ii+1)
            for aa in range(PBCorbFirst-1):
                for bb in range(len(S0)):
                    for dir in range(3): # x,y,z
                        for iSpin in range(len(H0)):
                            dH[(ii-FCfirst+1)*3+dir,iSpin,aa,bb]=0.0
                            dH[(ii-FCfirst+1)*3+dir,iSpin,bb,aa]=0.0
                    
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
    mm = N.dot
    Heph = 0.0*dH
    for i in range(len(hw)):
        # Loop over modes
        SIO.printDone(i, len(hw),'Calculating Heph')
        if hw[i]>0:
            for j in range(len(hw)):
                # Loop over atomic coordinates
                Heph[i] += const*dH[j]*U[i,j]/(2*PC.AtomicMass[atomnumber[FCfirst-1+j/3]]*hw[i])**.5
        else:
            print 'Phonons.CalcHeph: Nonpositive frequency --> Zero-valued coupling matrix' 
            Heph[i] = 0.0*dH[0]
        #Check that Heph is symmetric???
        for iSpin in range(len(dH[0])):
            if not N.allclose(Heph[i,iSpin,:,:],N.transpose(Heph[i,iSpin,:,:])):
                print 'Phonons.CalcHeph: WARNING: Coupling matrix not symmetric!'
    print '  ... Done!'
    return N.array(Heph)


def CorrectdH(onlySdir,orbitalIndices,nao,eF,H0,S0,dH,FCfirst,displacement):
    print 'Phonons.CorrectdH: Applying correction to dH...'
    dHnew = dH.copy()
    mm = N.dot
    onlyS0,dSx,dSy,dSz = GetOnlyS(onlySdir,nao,displacement)
    invS0 = LA.inv(S0)
    for i in range(len(dH)):
        SIO.printDone(i, len(dH),'Correcting dH')
        # Explicit correction
        dSdij = N.zeros((nao,nao),N.float)
        first,last = orbitalIndices[FCfirst-1+i/3]
        if i%3==0:   dSdij[:,first:last+1] = dSx[:,first:last+1]  # x-move
        elif i%3==1: dSdij[:,first:last+1] = dSy[:,first:last+1]  # y-move
        elif i%3==2: dSdij[:,first:last+1] = dSz[:,first:last+1]  # z-move
        for iSpin in range(len(H0)):
            dHnew[i,iSpin,:,:] = dH[i,iSpin,:,:] - mm(N.transpose(dSdij),mm(invS0,H0[iSpin,:,:])) - mm(H0[iSpin,:,:],mm(invS0,dSdij))
    print '  ... Done!'
    return dHnew


def GetOnlyS(onlySdir,nao,displacement):
    print 'Phonons.GetOnlyS: Reading from', onlySdir
    onlySfiles = glob.glob(onlySdir+'/*.onlyS*')
    onlySfiles.sort()
    if len(onlySfiles)<1:
        sys.exit('Phonons.GetOnlyS: No .onlyS file found!')
    # nao is number of orbitals in the original (not-doubled) structure
    S0 = N.zeros((nao,nao),N.float)
    dxm,dxp = N.zeros((nao,nao),N.float),N.zeros((nao,nao),N.float)
    dym,dyp = N.zeros((nao,nao),N.float),N.zeros((nao,nao),N.float)
    dzm,dzp = N.zeros((nao,nao),N.float),N.zeros((nao,nao),N.float) 
    if len(onlySfiles)!=6:
        sys.exit('Phonons.GetOnlyS: Wrong number of onlyS files found!')
    else:
        for file in onlySfiles:
            S = SIO.ReadOnlyS(file)
            orbitals = len(S)/2
            if orbitals!=nao:
                print 'nao=%i,  orbitals=%i'%(nao,orbitals)
                sys.exit('Phonons.GetOnlyS: Error assigning orbitals to atoms!')    
            S0=S[0:nao,0:nao].copy()
            dmat=S[0:nao,nao:nao*2].copy()
            if file.endswith('_1.onlyS'):
                dxm=dmat
            elif file.endswith('_2.onlyS'):
                dxp=dmat
            elif file.endswith('_3.onlyS'):
                dym=dmat
            elif file.endswith('_4.onlyS'):
                dyp=dmat
            elif file.endswith('_5.onlyS'):
                dzm=dmat
            elif file.endswith('_6.onlyS'):
                dzp=dmat
        thisd = 1e10
        for i in range(1,7):
            xyz = N.array(SIO.Getxyz(onlySdir+'/RUN_%i.fdf'%i))
            for j in range(1,len(xyz)):
                thisd = min(thisd,(N.dot(xyz[0]-xyz[j],xyz[0]-xyz[j]))**.5)
    # Check that onlyS-directory also corresponds to the same displacement
    print 'Phonons.GetOnlyS: OnlyS-displacement (min) = %.5f Ang'%thisd
    print 'Phonons.GetOnlyS: FC-displacement          = %.5f Ang'%displacement    
    if abs(thisd-displacement)/displacement > 0.05:  # Tolerate 5 percent off...
        sys.exit('Phonons.GetOnlyS: OnlyS-displacement different from FC-displacement!')    
    dSx = (dxp-dxm)/(2.*displacement)
    dSy = (dyp-dym)/(2.*displacement)
    dSz = (dzp-dzm)/(2.*displacement)
    return S0,dSx,dSy,dSz


def GetH0S0dH(tree,FCfirst,FClast,displacement,CorrPotentialShift=True):
    dH = []
    for i in range(len(tree)):
        # Go through each subdirectory
        first,last,dir = tree[i][0],tree[i][1],tree[i][2]
        HSfiles = tree[i][4]
        if len(HSfiles)!=(last-first+1)*6+1:
            sys.exit('Phonons.GetH0S0dH: Wrong number of *.TSHS files in %s\n'+\
                     'PROGRAM STOPPED!!!'%dir)
        # The first file is the one corresponding to no displacement
        eF,H0,S0,ia1,istep = SIO.ReadTSHSFile(HSfiles[0])
        if istep!=0:
            print "Phonons::GetH0S0dH Assumption on file order not right ",HSfiles[0]
            kuk
        for j in range(len(HSfiles)/2):
            if ia1+j/3 >= FCfirst and ia1+j/3 <= FClast:
                # Read TSHS file since it is within (FCfirst,FClast)
                eFm,tmpHm,tmpSm,ia1,istep = SIO.ReadTSHSFile(HSfiles[2*j+1])
                eFp,tmpHp,tmpSp,ia1,istep = SIO.ReadTSHSFile(HSfiles[2*j+2])
                if CorrPotentialShift:
                    for iSpin in range(len(H0)):
                        tmpHm[iSpin,:,:] -= (eFm-eF)*S0 # NB eF-shift multiplies to SO...
                        tmpHp[iSpin,:,:] -= (eFp-eF)*S0 
                dH.append((tmpHp-tmpHm)/(2*displacement))
    return eF,H0,S0,N.array(dH)

    
def OpenNetCDFFile(filename,nao,xyz,DeviceFirst,DeviceLast,FCfirst,FClast):
    print 'Phonons.WriteNetCDFFile: Writing', filename
    file = nc.NetCDFFile(filename,'w','Created '+time.ctime(time.time()))
    file.title = 'Output from Phonons.py'
    file.createDimension('AtomicOrbitals',int(nao))
    file.createDimension('XYZ',3)
    file.createDimension('One',1)
    file.createDimension('PhononModes',(FClast-FCfirst+1)*3)
    file.createDimension('NumTotAtoms',len(xyz))
    file.createDimension('NumDevAtoms',DeviceLast-DeviceFirst+1)
    file.createDimension('NumFCAtoms',FClast-FCfirst+1)
    return file


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
    tmp = []
    for i in range(len(hw)):
        tmp.append((hw[i],U[i].copy()))
    tmp.sort()
    tmp.reverse()
    print 'Phonons.CalcPhonons: Frequencies in meV:'
    for i in range(DynamicAtoms*3):
        hw[i],U[i] = tmp[i]
        print string.rjust('%.3f'%(1000*hw[i]),9),
        if (i-5)%6==0: print
    if (i-5)%6!=0: print
    print 'Phonons.CalcPhonons: Frequencies in cm^-1:'
    for i in range(DynamicAtoms*3):
        hw[i],U[i] = tmp[i]
        print string.rjust('%.3f'%(hw[i]/PC.invcm2eV),9),
        if (i-5)%6==0: print
    if (i-5)%6!=0: print


    return hw,U

def CalcPhononDOS(FC,U,atomnumber,FCfirst,FClast,BulkAtomsLeft,BulkAtomsRight,outputfile='PhDOS.dat'):
    import RecursiveMethods as RM
    mm = N.dot
    NumberOfAtoms = len(atomnumber)
    const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
    
    # Change FC to normal coordinates and include constants to make
    # the frequencies come out in eV.
    FCtilde = N.zeros((len(FC),len(FC)),N.float)
    for i in range(len(FC)):
        for j in range(NumberOfAtoms):
            for k in range(3):
                FCtilde[i,3*j+k] = (const**2)*FC[i,j,k]/ \
                  ( PC.AtomicMass[atomnumber[i/3]] * PC.AtomicMass[atomnumber[j]] )**0.5

    # Output the dense 3D-periodic FCtilde
    OutputFC(FCtilde)

    # Pick out 9 separate blocks from FCtilde
    iL = 3*BulkAtomsLeft
    iR = 3*(NumberOfAtoms-BulkAtomsRight)
    print iL,iR
    
    LL = FCtilde[:iL,:iL]
    LC = FCtilde[:iL,iL:iR]
    LR = FCtilde[:iL,iR:]
    CL = FCtilde[iL:iR,:iL]
    CC = FCtilde[iL:iR,iL:iR]
    CR = FCtilde[iL:iR,iR:]
    RL = FCtilde[iR:,:iL]
    RC = FCtilde[iR:,iL:iR]
    RR = FCtilde[iR:,iR:]

    OutputFC(LL,filename='FC_LL.matrix')
    OutputFC(LR,filename='FC_LR.matrix')
    OutputFC(RL,filename='FC_RL.matrix')
    OutputFC(RR,filename='FC_RR.matrix')
    
    ID = N.identity(len(FC))
    IDL = N.identity(len(LL))
    IDR = N.identity(len(RR))
    
    print N.shape(FC),N.shape(LL),N.shape(CC),N.shape(RR)


    # Build new FCmatrix from 7 of the 9 blocks, where
    # direct electrode couplings not included (LR=RL=0)
    FC2 = N.zeros((len(FC),len(FC)),N.float)
    
    FC2[:iL,:iL]     = LL
    FC2[:iL,iL:iR]   = LC
    FC2[iL:iR,:iL]   = CL
    FC2[iL:iR,iL:iR] = CC
    FC2[iL:iR,iR:]   = CR
    FC2[iR:,iL:iR]   = RC
    FC2[iR:,iR:]     = RR

    # Output this corrected FCmatrix of the finite system (in the z-direction)
    OutputFC(FC2,filename='FC2.matrix')
    
    # Calculate phonon self-energies (energy dependent)
    eta = 1e-6
    f = open(outputfile,'w')
    f2 = open('pdos.dat','w')
    FC2 = N.array(FC2,N.complex)
    for z in N.arange(0.0,30e-3,1e-4):
        z += eta*1j
        # Calculate self-energies
        L00,L01,L10,L11 = RM.IterateAlternating(z**2*IDL-LL,-LR,-RL,z**2*IDR-RR,
                                                iter=1e2,method='Fast',eta=0.0,mix=0.5)
        R00,R01,R10,R11 = RM.IterateAlternating(z**2*IDR-RR,-RL,-LR,z**2*IDL-LL,
                                                iter=1e2,method='Fast',eta=0.0,mix=0.5)
        # Insert self-energies in FC2
        FC2[:iL,:iL] = z**2*IDL-L00
        FC2[iR:,iR:] = z**2*IDR-R00

        # Calculate phonon density-of-states
        D = LA.inv(z**2*ID-FC2)  # The Phonon Greens function
        A = -4.0*z.real*D.imag       # Spectral density
        
        DOS = N.trace(A) # DOS projected on whole unit cell
        
        #dosZ = 0.0
        #for i in range(4):
        #    dosZ += -A[3*(27+i)+2,3*(27+i)+2]

        # Projection onto eigenmodes
        f2.write('%.7e  '%(1000*z.real))
        for vec in U:
            VEC = N.zeros(len(FC2),N.float)
            VEC[3*(FCfirst-1):3*FClast] = vec
            pdos = mm(VEC,mm(A,VEC)) 
            f2.write('%.3e  '%pdos)
        f2.write('\n')

        # Bulk and surface DOS    
        D = LA.inv(L00)
        dosSURFL = -4.0*z.real*N.trace(D.imag)
        D = LA.inv(L11)
        dosBULKL = -4.0*z.real*N.trace(D.imag)
        D = LA.inv(R00)
        dosSURFR = -4.0*z.real*N.trace(D.imag)
        D = LA.inv(R11)
        dosBULKR = -4.0*z.real*N.trace(D.imag)
        f.write('%.7f    %.3e    %.3e    %.3e    %.3e    %.3e \n'%(z.real,DOS,dosSURFL,dosBULKL,dosSURFR,dosBULKR))
    f.write('\n')
    f.close()
    f2.close()


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


def GetFileLists(dirname,wildcard):
    "Returns absolute paths to (FCfiles,TSHSfiles,XVfiles) matching the wildcard"
    dirs,tree,FCfiles,HSfiles,XVfiles = [],[],[],[],[]
    FCfirst,FClast = 1e10,-1e10
    # Find wildcard directories
    for elm in glob.glob(dirname+'/'+wildcard):
        if os.path.isdir(elm):
            dirs.append(elm)
    dirs.sort()
    for dir in dirs:
        localFCs,localHSs,localXVs = [],[],[]
        localFCfirst,localFClast = 1e10,-1e10
        # Find FCfiles in wildcard directories
        FCglob = glob.glob(dir+'/*.FC*')
        FCglob.sort()
        for elm in FCglob:
            if os.path.isfile(elm):
                FCfiles.append(elm)
                localFCs.append(elm)
        if len(FCglob)!=1:
            print 'Phonons.GetFileLists: Not exactly one *.FC file in directory',dir
        # Determine FCfirst and FClast without TSHS files
        runfdf = SIO.SIO_open(dir+'/RUN.fdf','r')
        lines = runfdf.readlines()
        for line in lines:
            if line.find('MD.FCfirst')>-1:
                localFCfirst = int(string.split(line)[1])
            if line.find('MD.FClast')>-1:
                localFClast = int(string.split(line)[1])
        runfdf.close()
        FCfirst = min(FCfirst,localFCfirst)
        FClast = max(FClast,localFClast)
        # Find TSHSfiles in wildcard directories
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
            sys.exit('Phonons.GetFileLists: Inconsistent number of *.TSHS files in directory %s'%dir)
        # Find XVfiles in wildcard directories
        for elm in glob.glob(dir+'/*.XV*'):
            if elm.endswith('.XV') or elm.endswith('.XV.gz'):
                XVfiles.append(elm)
        # Collect local directory information
        if len(localFCs)==1:
            tree.append([localFCfirst,localFClast,dir,localFCs[0],localHSs])
        else: print 'Phonons.GetFileLists: Two FC files in a directory encountered.'
    print 'Phonons.GetFileLists: %i folder(s) match wildcard %s' %(len(dirs),wildcard)
    print 'Phonons.GetFileLists: FCfirstMIN=%i, FClastMAX=%i, DynamicAtoms=%i' \
          %(FCfirst,FClast,FClast-FCfirst+1)
    tree.sort(),FCfiles.sort(),HSfiles.sort(),XVfiles.sort()
    return tree,XVfiles,FCfirst,FClast


def CorrectXVfile(XVfile):
    dir,file = os.path.split(XVfile)
    FCfirst, FClast = 1e10,-1e10
    # Determine FCfirst and FClast without TSHS files
    runfdf = SIO.SIO_open(dir+'/RUN.fdf','r')
    lines = runfdf.readlines()
    for line in lines:
        if line.find('MD.FCfirst')>-1:
            localFCfirst = int(string.split(line)[1])
        if line.find('MD.FClast')>-1:
            localFClast = int(string.split(line)[1])
    runfdf.close()
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
        if vectors1!=vectors2:
            count += 1
            print err, XVfileList[i],XVfileList[i+1], '(vectors) WARNING'
        if speciesnumber1!=speciesnumber2:
            count += 1
            print err, XVfileList[i],XVfileList[i+1], '(speciesnumber) WARNING'
        if atomnumber1!=atomnumber2:
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
                    ln += ' %.6f'%U[i][3*(j+1-FCfirst)+k]
            ln += '\n'
            f.write(ln)
    f.close()


def GenerateAuxNETCDF(tree,FCfirst,FClast,orbitalIndices,nao,onlySdir,PBCFirst,PBCLast,AuxNCfile,
                      displacement,CorrPotentialShift=True):
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
        eF,H0,S0,ia1,istep = SIO.ReadTSHSFile(HSfiles[0])
        # Read the rest of TSHSfiles
        for j in range(len(HSfiles)/2):
            if ia1+j/3 >= FCfirst and ia1+j/3 <= FClast:
                # Read TSHS file since it is within (FCfirst,FClast)
                eFm,tmpHm,tmpSm,ia1,istep = SIO.ReadTSHSFile(HSfiles[2*j+1])
                eFp,tmpHp,tmpSp,ia1,istep = SIO.ReadTSHSFile(HSfiles[2*j+2])
                if CorrPotentialShift:
                    tmpHm -= (eFm-eF)*S0 # NB eF-shift multiplies to SO...
                    tmpHp -= (eFp-eF)*S0
                try:
                    dH[index,:] = (tmpHp-tmpHm)/(2*displacement)
                    index += 1
                except:
                    # First time a netcdf-file is created
                    print 'Phonons.GenerateAuxNETCDF: Generating', AuxNCfile
                    NCfile2 = nc.NetCDFFile(AuxNCfile,'w','Created '+time.ctime(time.time()))
                    NCfile2.createDimension('Index',None)
                    NCfile2.createDimension('NSpin',len(H0))
                    NCfile2.createDimension('AtomicOrbitals',len(H0[0,:,:]))
                    NCfile2.createDimension('One',1)
                    FCfirsttmp = NCfile2.createVariable('FCfirst','d',('One',))
                    FCfirsttmp[:] = FCfirst
                    FClasttmp = NCfile2.createVariable('FClast','d',('One',))
                    FClasttmp[:] = FClast
                    H0tmp = NCfile2.createVariable('H0','d',('NSpin','AtomicOrbitals','AtomicOrbitals',))
                    print H0.shape
                    print H0tmp.shape
                    H0tmp[:] = H0
                    S0tmp = NCfile2.createVariable('S0','d',('AtomicOrbitals','AtomicOrbitals',))
                    S0tmp[:] = S0
                    dH = NCfile2.createVariable('dH','d',('Index','NSpin','AtomicOrbitals','AtomicOrbitals',))
                    dH[index,:] = (tmpHp-tmpHm)/(2*displacement)
                    index += 1
    NCfile2.sync()
    print 'Phonons.GenerateAuxNETCDF: len(dH) =',index
    
    # Correct dH
    mm = N.dot
    onlyS0,dSx,dSy,dSz = GetOnlyS(onlySdir,nao,displacement)
    invS0 = LA.inv(S0)
    for i in range(index):
        SIO.printDone(i,index,'Correcting dH')
        dSdij = N.zeros((nao,nao),N.float)
        first,last = orbitalIndices[FCfirst-1+i/3]
        if i%3==0:   dSdij[:,first:last+1] = dSx[:,first:last+1]  # x-move
        elif i%3==1: dSdij[:,first:last+1] = dSy[:,first:last+1]  # y-move
        elif i%3==2: dSdij[:,first:last+1] = dSz[:,first:last+1]  # z-move
        for iSpin in range(len(H0)):
            dH[i,iSpin,:] = dH[i,iSpin,:] - mm(N.transpose(dSdij),mm(invS0,H0[iSpin,:,:])) - mm(H0[iSpin,:,:],mm(invS0,dSdij))
    print '  ... Done!'
    NCfile2.sync()
    
    # Remove Periodic Boundary terms
    PBCorbFirst = orbitalIndices[PBCFirst-1][0]
    PBCorbLast  = orbitalIndices[PBCLast-1][1]    
    for ii in range(FCfirst-1,FClast):
        if ii+1<PBCFirst:
            # Left contact, remove contribution to right contact
            print "Phonons.GenerateAuxNETCDF: Removing He-ph, Left atom %i" % (ii+1)
            for aa in range(PBCorbLast+1,len(H0)):
                for bb in range(len(H0)):
                    for dir in range(3): # x,y,z
                        for iSpin in range(len(H0)):
                            dH[(ii-FCfirst+1)*3+dir,iSpin,aa,bb]=0.0
                            dH[(ii-FCfirst+1)*3+dir,iSpin,bb,aa]=0.0
        if ii+1>PBCLast:
            # Right contact, remove contribution to left contact
            print "Phonons.GenerateAuxNETCDF: Removing He-ph, Right atom %i" % (ii+1)
            for aa in range(PBCorbFirst-1):
                for bb in range(len(H0)):
                    for dir in range(3): # x,y,z
                        for iSpin in range(len(H0)):
                            dH[(ii-FCfirst+1)*3+dir,iSpin,aa,bb]=0.0
                            dH[(ii-FCfirst+1)*3+dir,iSpin,bb,aa]=0.0
    NCfile2.close()
    print 'Phonons.GenerateAuxNETCDF: File %s written.'%AuxNCfile


def CalcHephNETCDF(orbitalIndices,FCfirst,FClast,atomnumber,DeviceFirst,DeviceLast,
                   hw,U,NCfile,AuxNCfile):
    # Read AuxNCfile and downfold to device
    first,last = orbitalIndices[DeviceFirst-1][0],orbitalIndices[DeviceLast-1][1]
    NCfile2 = nc.NetCDFFile(AuxNCfile,'r')
    H0 = N.array(NCfile2.variables['H0'])[:,first:last+1,first:last+1]
    S0 = N.array(NCfile2.variables['S0'])[first:last+1,first:last+1]
    dH = NCfile2.variables['dH']

    NCfile.createDimension('NSpin',len(H0))
    Write2NetCDFFile(NCfile,H0,'H0',('NSpin','AtomicOrbitals','AtomicOrbitals',),units='eV')
    Write2NetCDFFile(NCfile,S0,'S0',('AtomicOrbitals','AtomicOrbitals',),units='eV')

    # Check AuxNCfile
    if FCfirst != int(NCfile2.variables['FCfirst'][0]):
        sys.exit('Phonons.CalcHephNETCDF: AuxNCfile %s does not correspond to FCfirst = %i'%(AuxNCfile,FCfirst))
    if FClast != int(NCfile2.variables['FClast'][0]):
        sys.exit('Phonons.CalcHephNETCDF: AuxNCfile %s does not correspond to FClast = %i'%(AuxNCfile,FClast))
        
    # CalcHeph
    print 'Phonons.CalcHephNETCDF: Calculating...\n',
    const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
    Heph  = NCfile.createVariable('He_ph','d',('PhononModes','NSpin','AtomicOrbitals','AtomicOrbitals',))
    for i in range(len(hw)):
        Heph[i,:] = 0.0*H0
    NCfile.sync()

    for i in range(len(hw)):
        # Loop over modes
        SIO.printDone(i, len(hw),'Calculating Heph')
        for j in range(len(hw)):
            # Loop over atomic coordinates
            if hw[i]>0:
                for iSpin in range(len(H0)):
                    Heph[i,iSpin,:] += const*N.array(dH[j,iSpin])[first:last+1,first:last+1]*U[i,j] \
                                       /(2*PC.AtomicMass[atomnumber[FCfirst-1+j/3]]*hw[i])**.5
            else:
                print 'Phonons.CalcHephNETCDF: Nonpositive frequency --> Zero-valued coupling matrix' 
                #Heph[i,:] = 0.0*N.array(dH[0])[first:last+1,first:last+1]
        #Check that Heph is symmetric???
        for iSpin in range(len(H0)):
            if not N.allclose(Heph[i,iSpin],N.transpose(Heph[i,iSpin])):
                print 'Phonons.CalcHephNETCDF: WARNING: Coupling matrix not symmetric!'
        NCfile.sync()
    print '  ... Done!'
    NCfile2.close()
    return H0,S0,Heph

