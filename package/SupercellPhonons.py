version = "SVN $Id$"
print version

"""
A simple interface to evaluate electron and phonon bands on 
a set of points in reciprocal space.

The input file format with N points is simply:

  kx(1) ky(1) kz(1) [label]
  kx(2) ky(2) kz(2) [label]
...
  kx(N) ky(N) kz(N) [label]

Units: 2pi/Ang.

Phase factors defined as: exp(2pi i k.r)

Thomas Frederiksen, March 2015
"""

import SiestaIO as SIO
import Symmetry
import CommonFunctions as CF
import Phonons as PH
import PhysicalConstants as PC
import MiscMath as MM
import WriteNetCDF as NCDF
import ValueCheck as VC

import numpy as N
import numpy.linalg as LA
import glob, os,sys,string
import scipy.linalg as SLA
import Scientific.IO.NetCDF as NC
    
vinfo = [version,SIO.version,Symmetry.version,CF.version,
         PH.version,PC.version,MM.version,NCDF.version,
         VC.version]

def GetOptions(argv,**kwargs):
    # if text string is specified, convert to list
    if type(argv)==type(''): argv = argv.split()

    import optparse as o

    usage = "usage: %prog [options] DestinationDirectory"
    description = "Methods to calculate electron and phonon band structures from finite-difference calculations"

    p = o.OptionParser(description=description,usage=usage)

    p.add_option("--FCwildcard",dest="FCwildcard",
                 help="Wildcard for FC directories [default=%default]",
                 type="str",default="./FC*")
    
    p.add_option("--OSdir",dest="onlySdir",
                 help="Location of OnlyS directory [default=%default]",
                 type="str",default="./OSrun")
    
    p.add_option("-r","--radius", dest="radius",
                 help="Force cutoff radius in Angstroms [default=%default]" ,
                 type="float",default=0.)

    p.add_option("--AtomicMass", dest='AtomicMass', default='[]',
                 help="Option to add to (or override!) existing dictionary of atomic masses. Format is a list [[anr1,mass1(,label)],...] [default=%default]")
    
    p.add_option("-k","--kpointfile",dest='kfile',default=None,
                 help="Input file with electronic k-points to be evaluated [default=%default]")

    p.add_option("-q","--qpointfile",dest='qfile',default=None,
                 help="Input file with phonon q-points to be evaluated [default=%default]")

    p.add_option('-s','--steps',dest='steps',default=100,type="int",
                 help="Number of points on path between high-symmetry k-points [default=%default]")

    p.add_option('--mesh',dest='mesh',default='[0,0,0]',type="str",
                 help="Mesh sampling of 1BZ (powers of 2) [default=%default]")

    (options, args) = p.parse_args(argv)

    # Get the last positional argument
    options.DestDir = VC.GetPositional(args,"You need to specify a destination directory")

    # With this one can overwrite the logging information
    if "log" in kwargs:
        options.Logfile = kwargs["log"]
    else:
        options.Logfile = 'Bandstructures.log'

    # Check if AtomicMasses are specified
    if options.AtomicMass!='[]':
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
  

class Supercell_DynamicalMatrix(PH.DynamicalMatrix):

    def CheckSymmetries(self):
        print '\nPerforming symmetry analysis'
        # Check if a smaller basis is present:
        Sym = Symmetry.Symmetry()
        # Find lattice symmetries
        Sym.setupGeom(self.geom.pbc,self.geom.snr,self.geom.anr,self.geom.xyz,onlyLatticeSym=True)
        # A primitive cell was found
        Sym.pointGroup()
        #Sym.findIrreducible()
        Sym.what()
        self.SetDynamicAtoms(range(1,Sym.basis.NN+1))
        # Calculate lattice vectors for phase factors
        # The closest cell might be different depending on which atom is moved
        sxyz = Sym.xyz.copy()
        latticevectors = N.zeros((Sym.NN,Sym.NN,3),N.float)
        for ii in range(Sym.NN):
            micxyz = Symmetry.moveIntoClosest(sxyz-sxyz[ii],Sym.pbc[0],Sym.pbc[1],Sym.pbc[2])
            for jj in range(Sym.NN):
                latticevectors[ii,jj] = micxyz[jj]+sxyz[Sym.basisatom[ii]]-sxyz[Sym.basisatom[jj]]           
        self.latticevectors = latticevectors
        self.supercell = True
        self.Sym = Sym

    def SymmetrizeFC(self,radius): 
        print '\nComputing symmetrized force constants'
        self.mean_sym = self.Sym.symmetrizeFC(self.mean,1,self.Sym.basis.NN,radi=radius)
        self.mean_sym = self.ApplySumRule(self.mean_sym)

    def ComputePhononModes_q(self,qpoint,verbose=True):
        # Compute phonon vectors
        if verbose:
            print '\nSupercellPhonons.ComputePhononModes_q: Computing force constants at q = ',qpoint,'(2pi/Ang)'
        NN = self.Sym.basis.NN
        self.q = N.zeros((NN,3,NN,3),N.complex)
        # Loop over basis atoms
        for n in range(NN):
            # Loop over all atoms in the supercell
            for m in range(len(self.latticevectors[n])):
                #print 'CV=',self.latticevectors[n,m]
                # exp( - i q R0m )
                R0m = self.latticevectors[n,m]
                phase = N.exp(-2.0j*N.pi*N.dot(qpoint,R0m))
                self.q[n,:,self.Sym.basisatom[m],:] += phase*self.mean_sym[n,:,m,:]
        # Now compute modes using standard module
        self.ComputePhononModes(self.q,verbose)
        # Expand U and Udispl to the whole supercell
        for n in range(len(self.latticevectors)):
            j = self.Sym.basisatom[n]
            R0n = self.latticevectors[0,n]
            phase = N.exp(-2.0j*N.pi*N.dot(qpoint,R0n))
            self.UU[:,n,:] = phase*self.U[:,3*j:3*j+3]
            self.UUdisp[:,n,:] = phase*self.Udisp[:,3*j:3*j+3]
        return self.hw,self.U

    def Fold2PrimitiveCell(self,H,kpoint):
        # Folding down H (or S) to specified kpoint
        # in the primitive cell
        sh = list(H.shape)
        sh[-1],sh[-2] = self.rednao,self.rednao
        H_k = N.zeros(tuple(sh),N.complex)
        # Loop over basis atoms
        for n in range(self.Sym.basis.NN):
            # Loop over all atoms in the supercell
            fn,ln = self.OrbIndx[n]
            fnb,lnb = self.OrbIndx[self.Sym.basisatom[n]] 
            #print 'n=%i, fn=%i, ln=%i, fnb=%i, lnb=%i'%(n,fn,ln,fnb,lnb)
            for m in range(len(self.latticevectors[n])):
                fm,lm = self.OrbIndx[m]
                fmb,lmb = self.OrbIndx[self.Sym.basisatom[m]]
                #print 'm=%i, fm=%i, lm=%i, fmb=%i, lmb=%i'%(m,fm,lm,fmb,lmb)
                # exp( i k R0m )
                R0m = self.latticevectors[n,m]
                phase = N.exp(2.0j*N.pi*N.dot(kpoint,R0m))
                H_k[...,fnb:lnb+1,fmb:lmb+1] += phase*H[...,fn:ln+1,fm:lm+1]
        return H_k

    def ComputeElectronStates(self,kpoint,verbose=True):
        if verbose:
            print 'SupercellPhonons.ComputeElectronStates: k = ',kpoint,'(2pi/Ang)'
        # Fold onto primitive cell
        self.h0_k = self.Fold2PrimitiveCell(self.h0,kpoint)
        self.s0_k = self.Fold2PrimitiveCell(self.s0,kpoint)
        ev, evec = SLA.eigh(self.h0_k[0],self.s0_k)
        #evecd = MM.dagger(evec)
        return ev,evec

def ReadKpoints(filename):
    print 'SupercellPhonons.ReadKpoints: Reading from',filename
    if filename.endswith('.nc'):
        return ReadKpoints_netcdf(filename)
    else:
        return ReadKpoints_ascii(filename)

def ReadKpoints_netcdf(filename):
    ncf = NC.NetCDFFile(filename,'r')
    kpts = ncf.variables['grid'][:]
    ncf.close()
    dk = N.empty(len(kpts))
    dk[0] = 0.0
    tmp = kpts[1:]-kpts[:-1]
    dk[1:] = N.diagonal(N.dot(tmp,tmp.T))**.5
    dk = N.cumsum(dk)
    labels,ticks = None,None
    return kpts,dk,labels,ticks
    
def ReadKpoints_ascii(filename):
    klist = []
    dklist = []
    labels = []
    ticks = []
    f = open(filename,'r')
    for ln in f.readlines():
        i = len(klist)
        s = ln.split()
        if len(s)==3 or len(s)==4:
            klist += [N.array([N.float(s[0]),N.float(s[1]),N.float(s[2])])]
            if i==0:
                dk = N.zeros(3,N.float)
                dklist += [N.dot(dk,dk)**.5]
            else:
                dk = klist[i]-klist[i-1]
                dklist += [dklist[i-1]+N.dot(dk,dk)**.5]
            if len(s)==4:
                labels += [s[3]]
                ticks += [[dklist[i],s[3]]]
            else:
                labels += ['']
    f.close()
    return N.array(klist),N.array(dklist),labels,ticks

def WriteKpoints(filename,klist,labels=None):
    f = open(filename,'w')
    for i in range(len(klist)):
        k = klist[i]
        for j in range(3):
            f.write('%.8e '%k[j])
        if labels:
            f.write(labels[i])
        f.write('\n')
    f.close()
    print 'SupercellPhonons.WriteKpoints: Wrote %i points to %s'%(len(klist),filename)

def WritePath(filename,path,steps):
    # Write path between high-symmetry points
    kpts = []
    labels = []
    for i,k in enumerate(path):
        if i<len(path)-1:
            k1 = path[i][0]
            k2 = path[i+1][0]
            for j in range(steps):
                kj = k1 + (k2-k1)*j/steps
                kpts.append(kj)
                if j==0:
                    labels.append(path[i][1])
                else:
                    labels.append('')
        else:
        # last k-point in path
            k1 = path[i][0]
            kpts.append(k1)
            labels.append(path[i][1])
    print 'High-symmetry path:'
    for k in path:
        print k[1],k[0]
    WriteKpoints(filename,kpts,labels)

def PlotElectronBands(filename,dk,elist,ticks):
    print elist.shape,dk.shape
    # Make xmgrace plots
    import WriteXMGR as XMGR
    if len(dk)>1:
        x = N.array([dk])
    else:
        x = N.array([[0.0]])
    print x.shape,elist.shape
    e = N.concatenate((x,elist.T)).T    
    es = XMGR.Array2XYsets(e,Lwidth=2,Lcolor=1)
    ge = XMGR.Graph(es)
    ge.SetXaxisSpecialTicks(ticks)
    ge.SetXaxis(max=dk[-1],majorGridlines=True)
    ge.SetYaxis(min=-20,max=20,label='E-E\sF\N (eV)',majorUnit=5.0)
    pe = XMGR.Plot(filename,ge)
    pe.WriteFile()

def PlotPhononBands(filename,dq,phlist,ticks):
    # Make xmgrace plots
    import WriteXMGR as XMGR
    if len(dq)>1:
        x = N.array([dq])
    else:
        x = N.array([[0.0]])
    p = N.concatenate((x,1000*phlist.T)).T
    ps = XMGR.Array2XYsets(p,Lwidth=2,Lcolor=1)
    gp = XMGR.Graph(ps)
    gp.SetXaxisSpecialTicks(ticks)
    gp.SetXaxis(max=dq[-1],majorGridlines=True)
    maxy = 1000*N.amax(phlist)
    if maxy<20: mu,mx = 5,20
    elif maxy<30: mu,mx = 5,30
    elif maxy<40: mu,mx = 5,40
    elif maxy<50: mu,mx = 10,50
    elif maxy<75: mu,mx = 10,75
    elif maxy<100: mu,mx = 20,100
    elif maxy<125: mu,mx = 25,125
    elif maxy<150: mu,mx = 25,150
    elif maxy<175: mu,mx = 25,175
    elif maxy<200: mu,mx = 25,200
    elif maxy<220: mu,mx = 25,220
    elif maxy<250: mu,mx = 25,250
    elif maxy<500: mu,mx = 100,500
    gp.SetYaxis(label='\\f{Symbol}w\\f{} (meV)',majorUnit=mu,min=0.0,max=mx)
    pp = XMGR.Plot(filename,gp)
    pp.WriteFile()
    
def ComputeDOS(ncfile,outfile,emin=0.0,emax=1.0,pts=1001,smear=1e-3):
    ncf = NC.NetCDFFile(ncfile,'r')
    bands = ncf.variables['eigenvalues'][:]
    egrid = N.linspace(emin,emax,pts)
    id1 = N.ones(bands.shape,N.float)
    id2 = N.ones(egrid.shape,N.float)
    dE = N.outer(egrid,id1)-N.outer(id2,bands) # [e,kn]
    w = N.exp(-dE**2/(2*smear**2))/(smear*(2*N.pi)**.5) # [e,b]
    dos = N.sum(w,axis=1)/len(bands) # sum over bands
    ncf.close()
    # Write plot
    import WriteXMGR as XMGR
    ps = XMGR.XYset(egrid,dos,Lwidth=2,Lcolor=1)
    gp = XMGR.Graph(ps)
    gp.SetXaxis(label='E (eV)',min=emin,max=emax,majorUnit=emax/5)
    ymax = N.max(dos)
    gp.SetYaxis(label='DOS (states/eV)',min=0,max=ymax,majorUnit=ymax/5)
    #gp.SetSubtitle(ncfile)
    pp = XMGR.Plot(outfile,gp)
    pp.PutText('smear = %.3f meV'%(1e3*smear),0.20,0.75)
    pp.PutText('kpts = %i'%len(bands),0.20,0.70)
    pp.WriteFile()

def main(options):
    CF.CreatePipeOutput(options.DestDir+'/'+options.Logfile)
    #VC.OptionsCheck(options,'Phonons')

    CF.PrintMainHeader('Bandstructures',vinfo,options)

    fdf = glob.glob(options.FCwildcard+'/RUN.fdf')  
    SCDM = Supercell_DynamicalMatrix(fdf)
    SCDM.CheckSymmetries()

    # Write high-symmetry path
    WritePath(options.DestDir+'/symmetry-path',SCDM.Sym.path,options.steps)

    # Write mesh
    k1,k2,k3 = eval(options.mesh)
    rvec = N.array([SCDM.Sym.b1,SCDM.Sym.b2,SCDM.Sym.b3])
    import Kmesh
    # Full mesh
    kmesh = Kmesh.kmesh(2**k1,2**k2,2**k3,meshtype=['LIN','LIN','LIN'],invsymmetry=False)
    WriteKpoints(options.DestDir+'/mesh_%ix%ix%i'%tuple(kmesh.Nk),N.dot(kmesh.k,rvec))
    # Mesh reduced by inversion symmetry 
    kmesh = Kmesh.kmesh(2**k1,2**k2,2**k3,meshtype=['LIN','LIN','LIN'],invsymmetry=True)
    WriteKpoints(options.DestDir+'/mesh_%ix%ix%i_invsym'%tuple(kmesh.Nk),N.dot(kmesh.k,rvec))
    
    # Evaluate electron k-points
    if options.kfile:
        # Prepare Hamiltonian etc in Gamma for whole supercell
        natoms = SIO.GetFDFlineWithDefault(fdf[0],'NumberOfAtoms',int,-1,'Error')
        SCDM.PrepareGradients(options.onlySdir,N.array([0.,0.,0.]),1,natoms,AbsEref=False,SinglePrec=True)
        SCDM.nao = SCDM.h0.shape[-1]
        SCDM.FirstOrb = SCDM.OrbIndx[0][0] # First atom = 1
        SCDM.LastOrb = SCDM.OrbIndx[SCDM.Sym.basis.NN-1][1] # Last atom = Sym.NN
        SCDM.rednao = SCDM.LastOrb+1-SCDM.FirstOrb
        # Compute electron eigenvalues
        kpts,dk,klabels,kticks = ReadKpoints(options.kfile)
        if klabels:
            # Only write ascii if labels exist
            WriteKpoints(options.DestDir+'/kpoints',kpts,klabels)
        # Loop over kpoints
        ncfn = options.DestDir+'/Electrons.nc'
        ncf = NC.NetCDFFile(ncfn, 'w')
        ncf.createDimension('gridpts',len(kpts))
        ncf.createDimension('kpt',3)
        grid = ncf.createVariable('grid','d',('gridpts','kpt'))
        grid[:] = kpts
        for i,k in enumerate(kpts):
            if i<100: # Print only for the first 100 points
                ev,evec = SCDM.ComputeElectronStates(k,verbose=True)
            else:
                ev,evec = SCDM.ComputeElectronStates(k,verbose=False)
                # otherwise something simple 
                if i%100==0: print '%i out of %i k-points computed'%(i,len(kpts))
            if i==0:
                # Initiate array
                ncf.createDimension('bands',len(ev))
                evals = ncf.createVariable('eigenvalues','d',('gridpts','bands'))
                evecsRe = ncf.createVariable('eigenvectors.re','d',('gridpts','bands','bands'))
                evecsIm = ncf.createVariable('eigenvectors.im','d',('gridpts','bands','bands'))
            evals[i] = ev
            evecsRe[i] = evec.real
            evecsIm[i] = evec.imag
        ncf.sync()
        # Produce nice plots if labels exist
        if klabels:
            PlotElectronBands(options.DestDir+'/Electrons.agr',dk,N.array(evals[:]),kticks)
        ncf.close()


    # Compute phonon eigenvalues
    if options.qfile:
        SCDM.SymmetrizeFC(options.radius)
        SCDM.SetMasses()
        qpts,dq,qlabels,qticks = ReadKpoints(options.qfile)
        if qlabels:
            # Only write ascii if labels exist
            WriteKpoints(options.DestDir+'/qpoints',qpts,qlabels)
        # Loop over qpoints
        ncfn = options.DestDir+'/Phonons.nc'
        ncf = NC.NetCDFFile(ncfn, 'w')
        ncf.createDimension('gridpts',len(qpts))
        ncf.createDimension('kpt',3)
        grid = ncf.createVariable('grid','d',('gridpts','kpt'))
        grid[:] = qpts        
        for i,q in enumerate(qpts):
            if i<100: # Print only for the first 100 points
                hw,U = SCDM.ComputePhononModes_q(q,verbose=True)
            else:
                hw,U = SCDM.ComputePhononModes_q(q,verbose=False)
                # otherwise something simple
                if i%100==0: print '%i out of %i q-points computed'%(i,len(qpts))
            if i==0:
                # Initiate array
                ncf.createDimension('bands',len(hw))
                evals = ncf.createVariable('eigenvalues','d',('gridpts','bands'))
                evecsRe = ncf.createVariable('eigenvectors.re','d',('gridpts','bands','bands'))
                evecsIm = ncf.createVariable('eigenvectors.im','d',('gridpts','bands','bands'))
            evals[i] = hw
            evecsRe[i] = U.real
            evecsIm[i] = U.imag
        ncf.sync()
        # Produce nice plots if labels exist
        if qlabels:
            PlotPhononBands(options.DestDir+'/Phonons.agr',dq,N.array(evals[:]),qticks)
        ncf.close()        
