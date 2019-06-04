"""

:mod:`Inelastica.SupercellPhonons`
==================================

A simple interface to evaluate electron and phonon bands on
a set of points in reciprocal space.

The input file format with `N` points is simply:

.. code-block:: bash

    kx(1) ky(1) kz(1) [label]
    kx(2) ky(2) kz(2) [label]
    ...
    kx(N) ky(N) kz(N) [label]

**Units:** 1/Ang.

Phase factors defined as: `exp(i k.r)`

Thomas Frederiksen, March 2015.

Classes
-------

.. autosummary::
   :toctree:

   Supercell_DynamicalMatrix

.. currentmodule:: Inelastica.SupercellPhonons

"""
from __future__ import print_function

import numpy as N
import glob
import scipy.linalg as SLA
import netCDF4 as NC4
import ast
import Inelastica.io.siesta as SIO
import Inelastica.Symmetry as Symmetry
import Inelastica.io.log as Log
import Inelastica.Phonons as PH
import Inelastica.physics.constants as PC
import Inelastica.math as MM
import Inelastica.misc.valuecheck as VC
import Inelastica.io.xmgrace as XMGR

__all__ = ['Supercell_DynamicalMatrix']


def GetOptions(argv):
    # if text string is specified, convert to list
    if isinstance(argv, VC.string_types):
        argv = argv.split()

    import argparse

    p = argparse.ArgumentParser(description='Methods to calculate electron and phonon band structures from finite-difference calculations')
    p.add_argument('DestDir', help='Destination directory')
    p.add_argument('--FCwildcard', dest='FCwildcard', type=str, default='./FC*',
                 help='Wildcard for FC directories [default=%(default)s]')
    p.add_argument('--OSdir', dest='onlySdir', type=str, default='./OSrun',
                 help='Location of OnlyS directory [default=%(default)s]')
    p.add_argument('-r', '--radius', dest='radius', type=float, default=0.,
                 help='Force cutoff radius in Angstroms [default=%(default)s]')
    p.add_argument('--AtomicMass', dest='AtomicMass', default='[]',
                 help='Option to add to (or override!) existing dictionary of atomic masses. Format is a list [[anr1,mass1(,label)],...] [default=%(default)s]')
    p.add_argument('-k', '--kpointfile', dest='kfile', default=None,
                 help='Input file with electronic k-points to be evaluated [default=%(default)s]')
    p.add_argument('-q', '--qpointfile', dest='qfile', default=None,
                 help='Input file with phonon q-points to be evaluated [default=%(default)s]')
    p.add_argument('-s', '--steps', dest='steps', default=100, type=int,
                 help='Number of points on path between high-symmetry k-points [default=%(default)s]')
    p.add_argument('--mesh', dest='mesh', default='[0,0,0]', type=str,
                 help='Mesh sampling over one BZ (powers of 2) [default=%(default)s]')
    p.add_argument('--sort', dest='sorting', action='store_true', default=False,
                 help='Sort eigenvalues along k-mesh for nice plots? [default=%(default)s]')
    p.add_argument('--TSdir', dest='onlyTSdir', type=str, default=None,
                 help='Location of TranSIESTA calculation directory (will ignore FC and OnlyS directories) [default=%(default)s]')
    p.add_argument('--nbands', dest='nbands', type=int, default=None,
                 help='Number of electronic bands to be included in netCDF output (lower energy bands) [default=%(default)s]')
    p.add_argument('--FermiSurface', dest='FermiSurface', action='store_true', default=False,
                 help='Write FermiSurface.BXSF file for visualization of the Fermi surface? [default=%(default)s]')

    options = p.parse_args(argv)

    # Set module name
    options.module = 'Bandstructures'

    # Check if AtomicMasses are specified
    if options.AtomicMass != '[]':
        masslist = ast.literal_eval(options.AtomicMass.replace('\n', '').replace(' ', ''))
        for elm in masslist:
            anr = int(elm[0])
            mass = float(elm[1])
            PC.AtomicMass[anr] = mass
            if len(elm) == 3:
                label = elm[2]
                PC.PeriodicTable[anr] = label
                PC.PeriodicTable[label] = anr
        print('AtomicMass =', PC.AtomicMass)
        print('PeriodicTable =', PC.PeriodicTable)

    return options


class Supercell_DynamicalMatrix(PH.DynamicalMatrix):

    def __init__(self, fdf, TSrun=False):
        PH.DynamicalMatrix.__init__(self, fdf, DynamicAtoms=None, TSrun=TSrun)
        # Find basis and symmetry operations
        self.CheckSymmetries(TSrun=TSrun)

    # PHONONIC PART

    def CheckSymmetries(self, TSrun=False):
        print('\nPerforming symmetry analysis')
        # Check if a smaller basis is present:
        Sym = Symmetry.Symmetry()
        # Find lattice symmetries
        Sym.setupGeom(self.geom.pbc, self.geom.snr, self.geom.anr, self.geom.xyz, onlyLatticeSym=True)
        if not TSrun:
            Sym.pointGroup() # Actually this call is only needed if phonons are computed
            #Sym.findIrreducible()
            self.SetDynamicAtoms(list(range(1, Sym.basis.NN+1)))
        Sym.what()
        # Calculate lattice vectors for phase factors
        # The closest cell might be different depending on which atom is moved
        sxyz = Sym.xyz.copy()
        latticevectors = N.zeros((Sym.NN, Sym.NN, 3), N.float)
        for ii in range(Sym.NN):
            micxyz = Symmetry.moveIntoClosest(sxyz-sxyz[ii], Sym.pbc[0], Sym.pbc[1], Sym.pbc[2])
            for jj in range(Sym.NN):
                latticevectors[ii, jj] = micxyz[jj]+sxyz[Sym.basisatom[ii]]-sxyz[Sym.basisatom[jj]]
        self.latticevectors = latticevectors
        self.supercell = True
        self.Sym = Sym

    def SymmetrizeFC(self, radius):
        print('\nComputing symmetrized force constants')
        self.mean_sym = self.Sym.symmetrizeFC(self.mean, 1, self.Sym.basis.NN, radi=radius)
        self.mean_sym = self.ApplySumRule(self.mean_sym)

    def ComputePhononModes_q(self, qpoint, verbose=True):
        # Compute phonon vectors
        if verbose:
            print('\nSupercellPhonons.ComputePhononModes_q: Computing force constants at q = ', qpoint, '(1/Ang)')
        NN = self.Sym.basis.NN
        self.q = N.zeros((NN, 3, NN, 3), N.complex)
        # Loop over basis atoms
        for n in range(NN):
            # Loop over all atoms in the supercell
            for m in range(len(self.latticevectors[n])):
                #print 'CV=',self.latticevectors[n,m]
                # exp( - i q R0m )
                R0m = self.latticevectors[n, m]
                phase = N.exp(-1.0j*N.dot(qpoint, R0m))
                self.q[n, :, self.Sym.basisatom[m], :] += phase*self.mean_sym[n, :, m, :]
        # Now compute modes using standard module
        self.ComputePhononModes(self.q, verbose)
        # Expand U and Udispl to the whole supercell
        for n in range(len(self.latticevectors)):
            j = self.Sym.basisatom[n]
            R0n = self.latticevectors[0, n]
            phase = N.exp(-1.0j*N.dot(qpoint, R0n))
            self.UU[:, n, :] = phase*self.U[:, 3*j:3*j+3]
            self.UUdisp[:, n, :] = phase*self.Udisp[:, 3*j:3*j+3]
        return self.hw, self.U

    # ELECTRONIC PART

    def Fold2PrimitiveCell(self, H, kpoint):
        # Folding down H (or S) to specified kpoint
        # in the primitive cell
        sh = list(H.shape)
        sh[-1], sh[-2] = self.rednao, self.rednao
        H_k = N.zeros(tuple(sh), N.complex)
        # Loop over basis atoms
        for n in range(self.Sym.basis.NN):
            # Loop over all atoms in the supercell
            fn, ln = self.OrbIndx[n]
            fnb, lnb = self.OrbIndx[self.Sym.basisatom[n]]
            #print 'n=%i, fn=%i, ln=%i, fnb=%i, lnb=%i'%(n,fn,ln,fnb,lnb)
            for m in range(len(self.latticevectors[n])):
                fm, lm = self.OrbIndx[m]
                fmb, lmb = self.OrbIndx[self.Sym.basisatom[m]]
                #print 'm=%i, fm=%i, lm=%i, fmb=%i, lmb=%i'%(m,fm,lm,fmb,lmb)
                # exp( i k R0m )
                R0m = self.latticevectors[n, m]
                phase = N.exp(1.0j*N.dot(kpoint, R0m))
                #H_k[...,fnb:lnb+1,fmb:lmb+1] += phase*H[...,fn:ln+1,fm:lm+1]
                H_k[..., fmb:lmb+1, fnb:lnb+1] += phase*H[..., fm:lm+1, fn:ln+1]
        return H_k

    def ComputeElectronStates(self, kpoint, verbose=True, TSrun=False):
        if TSrun:
            # kpoint has unit of '2*pi/a'
            kpt = kpoint/(2*N.pi)
            kpt2 = MM.mm(N.array([self.Sym.a1, self.Sym.a2, self.Sym.a3]), kpt)
            self.TSHS0.setkpoint(kpt2, atype=N.complex, verbose=verbose)
            self.h0_k = self.TSHS0.H[:, :, :]
            self.s0_k = self.TSHS0.S[:, :]
        else:
            if verbose:
                print('SupercellPhonons.ComputeElectronStates: k = ', kpoint, '(1/Ang)')
            # Fold onto primitive cell
            self.h0_k = self.Fold2PrimitiveCell(self.h0, kpoint)
            self.s0_k = self.Fold2PrimitiveCell(self.s0, kpoint)
        ev = N.empty((self.nspin, self.rednao), N.float)
        evec = N.empty((self.nspin, self.rednao, self.rednao), N.complex)
        for ispin in range(self.nspin):
            ev[ispin], evec[ispin] = SLA.eigh(self.h0_k[ispin], self.s0_k)
        return ev, evec

    def ReadGradients(self):
        # Read in gradients once into memory (in the full supercell basis)
        # No folding yet onto k and q
        self.dH = {}
        # Loop over dynamic atoms
        for v in self.DynamicAtoms:
            # Loop over axes
            for j in range(3):
                # Compute gradient
                self.dH[v, j] = self.GetGradient(v, j)

    # ELECTRON-PHONON COUPLING PART

    def FoldMatrix_DoubleSum(self, H, kpoint, qpoint):
        sh = list(H.shape)
        sh[-1], sh[-2] = self.rednao, self.rednao
        H_kq = N.zeros(tuple(sh), N.complex)
        # Loop over atoms
        for n in range(len(self.latticevectors[0])):
            fn, ln = self.OrbIndx[n]
            fnb, lnb = self.OrbIndx[self.Sym.basisatom[n]]
            # Loop over atoms
            for m in range(len(self.latticevectors[0])):
                fm, lm = self.OrbIndx[m]
                fmb, lmb = self.OrbIndx[self.Sym.basisatom[m]]
                # exp( i q R0m )
                R0m = self.latticevectors[self.Sym.basisatom[n], m]
                phase1 = N.exp(1.0j*N.dot(qpoint, R0m))
                # exp( i k Rnm )
                Rnm = self.latticevectors[n, m]
                #print 'allclose',N.allclose(Rnm,-1.* self.latticevectors[m,n]) # this is always true
                phase2 = N.exp(1.0j*N.dot(kpoint, Rnm))
                #H_kq[...,fnb:lnb+1,fmb:lmb+1] += phase1*phase2*H[...,fn:ln+1,fm:lm+1]
                H_kq[..., fmb:lmb+1, fnb:lnb+1] += phase1*phase2*H[..., fm:lm+1, fn:ln+1]
        #return H_kq/self.geom.natoms*self.Sym.basis.NN
        return H_kq

    def ComputeEPHcouplings_kq(self, kpoint, qpoint, verbose=True):
        if verbose:
            print('SupercellPhonons.ComputeEPHcouplings_kq: k = ', kpoint, ' q = ', qpoint, '(1/Ang)')
        # This assumes that phonon modes UU has been computed at that q
        const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
        M = N.zeros((len(self.hw), self.nspin, self.nao, self.nao), N.complex)
        G = N.zeros((len(self.hw), self.nspin, self.nao, self.nao), N.complex)
        # Loop over dynamic atoms
        for i, v in enumerate(self.DynamicAtoms):
            # Loop over axes
            for j in range(3):
                # Loop over modes
                for m in range(len(self.hw)):
                    M[m] += self.dH[v, j]*self.UU[m, v-1, j]
                    if self.hw[m] > 0.0: # Only proceed with a postive frequency
                        G[m] += const*self.dH[v, j]*self.UU[m, v-1, j]/(2*self.Masses[i]*self.hw[m])**.5
        m = self.FoldMatrix_DoubleSum(M, kpoint, qpoint)
        g = self.FoldMatrix_DoubleSum(G, kpoint, qpoint)
        return m, g

# AUXILIARY FUNCTIONS


def ReadKpoints(filename):
    print('SupercellPhonons.ReadKpoints: Reading from', filename)
    if filename.endswith('.nc'):
        return ReadKpoints_netcdf(filename)
    else:
        return ReadKpoints_ascii(filename)


def ReadKpoints_netcdf(filename):
    ncf = NC4.Dataset(filename, 'r')
    kpts = ncf.variables['grid'][:]
    ncf.close()
    dk = N.empty(len(kpts))
    dk[0] = 0.0
    tmp = kpts[1:]-kpts[:-1]
    dk[1:] = N.diagonal(N.dot(tmp, tmp.T))**.5
    dk = N.cumsum(dk)
    labels, ticks = None, None
    return kpts, dk, labels, ticks


def ReadKpoints_ascii(filename):
    klist = []
    dklist = []
    labels = []
    ticks = []
    f = open(filename, 'r')
    # Initialize variables with the first read k-point
    ln = f.readlines()[0]
    s = ln.split()
    if len(s) == 3 or len(s) == 4:
        klist += [N.array([N.float(s[0]), N.float(s[1]), N.float(s[2])])]
        dk = N.zeros(3, N.float)
        dklist += [N.dot(dk, dk)**.5]
        if len(s) == 4:
            labels += [s[3]]
            ticks += [[dklist[0], s[3]]]
        else:
            labels += ['']
    # Continue reading the rest of the file
    f.seek(0)
    for ln in f.readlines()[1:]:
        i = len(klist)
        s = ln.split()
        if len(s) == 3 or len(s) == 4:
            klist += [N.array([N.float(s[0]), N.float(s[1]), N.float(s[2])])]
            dk = klist[i]-klist[i-1]
            dklist += [dklist[i-1]+N.dot(dk, dk)**.5]
            if len(s) == 4:
                labels += [s[3]]
                ticks += [[dklist[i], s[3]]]
            else:
                labels += ['']
    f.close()
    return N.array(klist), N.array(dklist), labels, ticks


def WriteKpoints(filename, klist, labels=None):
    f = open(filename, 'w')
    for i, kval in enumerate(klist):
        k = kval
        for j in range(3):
            f.write('%.8e '%k[j])
        if labels:
            f.write(labels[i])
        f.write('\n')
    f.close()
    print('SupercellPhonons.WriteKpoints: Wrote %i points to %s'%(len(klist), filename))


def WritePath(filename, path, steps):
    # Write path between high-symmetry points
    kpts = []
    labels = []
    for i, k in enumerate(path):
        if i < len(path)-1:
            k1 = path[i][0]
            k2 = path[i+1][0]
            for j in range(steps):
                kj = k1 + (k2-k1)*j/steps
                kpts.append(kj)
                if j == 0:
                    labels.append(path[i][1])
                else:
                    labels.append('')
        else:
        # last k-point in path
            k1 = path[i][0]
            kpts.append(k1)
            labels.append(path[i][1])
    print('High-symmetry path:')
    for k in path:
        print(k[1], k[0])
    WriteKpoints(filename, kpts, labels)


def SortBands(ev):
    # Sort bands by curvature minimization
    print('SupercellPhonons.SortBands: Minimizing curvature in band structure')
    kpts, bands = ev.shape
    # loop over k-points, starting from the third
    for i in range(2, kpts):
        # loop over band index
        for j in range(bands):
            # loop over higher-lying bands
            for k in range(j+1, bands):
                d2ev = abs(ev[i, j]+ev[i-2, j]-2*ev[i-1, j])+abs(ev[i, k]+ev[i-2, k]-2*ev[i-1, k])
                d2evfl = abs(ev[i, k]+ev[i-2, j]-2*ev[i-1, j])+abs(ev[i, j]+ev[i-2, k]-2*ev[i-1, k])
                if d2ev > d2evfl:
                    # flip bands
                    tmp = ev[i:, j].copy()
                    ev[i:, j] = ev[i:, k]
                    ev[i:, k] = tmp
                    print("   ...@ k(%i): flip indices %i-%i"%(i, j, k))
    return ev


def PlotElectronBands(filename, dk, elist, ticks):
    # Make xmgrace plots
    if len(dk) > 1:
        x = N.array([dk])
    else:
        x = N.array([[0.0]])
    e = N.concatenate((x, elist.T)).T
    es = XMGR.Array2XYsets(e, Lwidth=2, Lcolor=1)
    ge = XMGR.Graph(es)
    ge.SetXaxisSpecialTicks(ticks)
    ge.SetXaxis(vmax=dk[-1], majorGridlines=True)
    ge.SetYaxis(vmin=-20, vmax=20, label=r'E-E\sF\N (eV)', majorUnit=5.0)
    pe = XMGR.Plot(filename, ge)
    pe.WriteFile()


def PlotPhononBands(filename, dq, phlist, ticks):
    # Make xmgrace plots
    if len(dq) > 1:
        x = N.array([dq])
    else:
        x = N.array([[0.0]])
    p = N.concatenate((x, 1000*phlist.T)).T
    ps = XMGR.Array2XYsets(p, Lwidth=2, Lcolor=1)
    gp = XMGR.Graph(ps)
    gp.SetXaxisSpecialTicks(ticks)
    gp.SetXaxis(vmax=dq[-1], majorGridlines=True)
    maxy = 1000*N.amax(phlist)
    if maxy < 20: mu, mx = 5, 20
    elif maxy < 30: mu, mx = 5, 30
    elif maxy < 40: mu, mx = 5, 40
    elif maxy < 50: mu, mx = 10, 50
    elif maxy < 75: mu, mx = 10, 75
    elif maxy < 100: mu, mx = 20, 100
    elif maxy < 125: mu, mx = 25, 125
    elif maxy < 150: mu, mx = 25, 150
    elif maxy < 175: mu, mx = 25, 175
    elif maxy < 200: mu, mx = 25, 200
    elif maxy < 220: mu, mx = 25, 220
    elif maxy < 250: mu, mx = 25, 250
    elif maxy < 500: mu, mx = 100, 500
    gp.SetYaxis(label='\\f{Symbol}w\\f{} (meV)', majorUnit=mu, vmin=0.0, vmax=mx)
    pp = XMGR.Plot(filename, gp)
    pp.WriteFile()


def ComputeDOS(ncfile, outfile, emin=0.0, emax=1.0, pts=1001, smear=1e-3):
    ncf = NC4.Dataset(ncfile, 'r')
    ev = ncf.variables['eigenvalues'][:]
    if len(ev.shape) == 2: # Phonons.nc (gridpts, bands)
        WriteDOS(outfile, ev, emin, emax, pts, smear)
    elif len(ev.shape) == 3: # Electrons.nc (gridpts, nspin, bands)
        if len(ev[0]) == 1: # nspin==1
            WriteDOS(outfile, ev[:, 0], emin, emax, pts, smear)
        elif len(ev[0]) == 2: # nspin==2
            WriteDOS(outfile+'.UP', ev[:, 0], emin, emax, pts, smear)
            WriteDOS(outfile+'.DOWN', ev[:, 1], emin, emax, pts, smear)
    ncf.close()


def WriteDOS(outfile, bands, emin, emax, pts, smear):
    egrid = N.linspace(emin, emax, pts)
    id1 = N.ones(bands.shape, N.float)
    id2 = N.ones(egrid.shape, N.float)
    dE = N.outer(egrid, id1)-N.outer(id2, bands) # [e,kn]
    w = N.exp(-dE**2/(2*smear**2))/(smear*(2*N.pi)**.5) # [e,b]
    dos = N.sum(w, axis=1)/len(bands) # sum over bands
    # Write plot
    ps = XMGR.XYset(egrid, dos, Lwidth=2, Lcolor=1)
    gp = XMGR.Graph(ps)
    gp.SetXaxis(label='E (eV)', vmin=emin, vmax=emax, majorUnit=emax/5)
    ymax = N.max(dos)
    gp.SetYaxis(label='DOS (states/eV)', vmin=0, vmax=ymax, majorUnit=ymax/5)
    #gp.SetSubtitle(ncfile)
    pp = XMGR.Plot(outfile, gp)
    pp.PutText('smear = %.3f meV'%(1e3*smear), 0.20, 0.75)
    pp.PutText('kpts = %i'%len(bands), 0.20, 0.70)
    pp.WriteFile()


def main(options):
    Log.CreatePipeOutput(options)
    #VC.OptionsCheck(options)
    Log.PrintMainHeader(options)

    try:
        fdf = glob.glob(options.onlyTSdir+'/RUN.fdf')
        TSrun = True
    except:
        fdf = glob.glob(options.FCwildcard+'/RUN.fdf') # This should be made an input flag
        TSrun = False
    SCDM = Supercell_DynamicalMatrix(fdf, TSrun)

    # Write high-symmetry path
    WritePath(options.DestDir+'/symmetry-path', SCDM.Sym.path, options.steps)

    # Write mesh
    k1, k2, k3 = ast.literal_eval(options.mesh)
    rvec = 2*N.pi*N.array([SCDM.Sym.b1, SCDM.Sym.b2, SCDM.Sym.b3])
    import Inelastica.physics.mesh as Kmesh
    # Full mesh
    kmesh = Kmesh.kmesh(2**k1, 2**k2, 2**k3, meshtype=['LIN', 'LIN', 'LIN'], invsymmetry=False)
    WriteKpoints(options.DestDir+'/mesh_%ix%ix%i'%tuple(kmesh.Nk), N.dot(kmesh.k, rvec))
    # Mesh reduced by inversion symmetry
    kmesh = Kmesh.kmesh(2**k1, 2**k2, 2**k3, meshtype=['LIN', 'LIN', 'LIN'], invsymmetry=True)
    WriteKpoints(options.DestDir+'/mesh_%ix%ix%i_invsym'%tuple(kmesh.Nk), N.dot(kmesh.k, rvec))

    # Evaluate electron k-points
    if options.kfile:
        # Prepare Hamiltonian etc in Gamma for whole supercell
        natoms = SIO.GetFDFlineWithDefault(fdf[0], 'NumberOfAtoms', int, -1, 'Error')
        SCDM.PrepareGradients(options.onlySdir, N.array([0., 0., 0.]), 1, natoms, AbsEref=False, atype=N.complex, TSrun=TSrun)
        SCDM.nao = SCDM.h0.shape[-1]
        SCDM.FirstOrb = SCDM.OrbIndx[0][0] # First atom = 1
        SCDM.LastOrb = SCDM.OrbIndx[SCDM.Sym.basis.NN-1][1] # Last atom = Sym.NN
        SCDM.rednao = SCDM.LastOrb+1-SCDM.FirstOrb
        # Read kpoints
        kpts, dk, klabels, kticks = ReadKpoints(options.kfile)
        if klabels:
            # Only write ascii if labels exist
            WriteKpoints(options.DestDir+'/kpoints', kpts, klabels)
        # Prepare netcdf
        ncfn = options.DestDir+'/Electrons.nc'
        ncf = NC4.Dataset(ncfn, 'w')
        # Grid
        ncf.createDimension('gridpts', len(kpts))
        ncf.createDimension('vector', 3)
        grid = ncf.createVariable('grid', 'd', ('gridpts', 'vector'))
        grid[:] = kpts
        grid.units = '1/Angstrom'
        # Geometry
        ncf.createDimension('atoms', SCDM.Sym.basis.NN)
        xyz = ncf.createVariable('xyz', 'd', ('atoms', 'vector'))
        xyz[:] = SCDM.Sym.basis.xyz
        xyz.units = 'Angstrom'
        pbc = ncf.createVariable('pbc', 'd', ('vector', 'vector'))
        pbc.units = 'Angstrom'
        pbc[:] = [SCDM.Sym.a1, SCDM.Sym.a2, SCDM.Sym.a3]
        rvec1 = ncf.createVariable('rvec', 'd', ('vector', 'vector'))
        rvec1.units = '1/Angstrom (incl. factor 2pi)'
        rvec1[:] = rvec
        ncf.sync()
        # Loop over kpoints
        for i, k in enumerate(kpts):
            if i < 100: # Print only for the first 100 points
                ev, evec = SCDM.ComputeElectronStates(k, verbose=True, TSrun=TSrun)
            else:
                ev, evec = SCDM.ComputeElectronStates(k, verbose=False, TSrun=TSrun)
                # otherwise something simple
                if i%100 == 0: print('%i out of %i k-points computed'%(i, len(kpts)))
            if i == 0:
                ncf.createDimension('nspin', SCDM.nspin)
                ncf.createDimension('orbs', SCDM.rednao)
                if options.nbands and options.nbands < SCDM.rednao:
                    nbands = options.nbands
                else:
                    nbands = SCDM.rednao
                ncf.createDimension('bands', nbands)
                evals = ncf.createVariable('eigenvalues', 'd', ('gridpts', 'nspin', 'bands'))
                evals.units = 'eV'
                evecsRe = ncf.createVariable('eigenvectors.re', 'd', ('gridpts', 'nspin', 'orbs', 'bands'))
                evecsIm = ncf.createVariable('eigenvectors.im', 'd', ('gridpts', 'nspin', 'orbs', 'bands'))
                # Check eigenvectors
                print('SupercellPhonons: Checking eigenvectors at', k)
                for j in range(SCDM.nspin):
                    ev2 = N.diagonal(MM.mm(MM.dagger(evec[j]), SCDM.h0_k[j], evec[j]))
                    print(' ... spin %i: Allclose='%j, N.allclose(ev[j], ev2, atol=1e-5, rtol=1e-3))
                ncf.sync()
            # Write to NetCDF
            evals[i, :] = ev[:, :nbands]
            evecsRe[i, :] = evec[:, :, :nbands].real
            evecsIm[i, :] = evec[:, :, :nbands].imag
        ncf.sync()
        # Include basis orbitals in netcdf file
        if SCDM.Sym.basis.NN == len(SCDM.OrbIndx):
            lasto = N.zeros(SCDM.Sym.basis.NN+1, N.float)
            lasto[:SCDM.Sym.basis.NN] = SCDM.OrbIndx[:SCDM.Sym.basis.NN, 0]
            lasto[SCDM.Sym.basis.NN] = SCDM.OrbIndx[SCDM.Sym.basis.NN-1, 1]+1
        else:
            lasto = SCDM.OrbIndx[:SCDM.Sym.basis.NN+1, 0]
        orbbasis = SIO.BuildBasis(fdf[0], 1, SCDM.Sym.basis.NN, lasto)
        # Note that the above basis is for the geometry with an atom FC-moved in z.
        #print dir(orbbasis)
        #print orbbasis.xyz # Hence, this is not the correct geometry of the basis atoms!
        center = ncf.createVariable('orbcenter', 'i', ('orbs',))
        center[:] = N.array(orbbasis.ii-1, dtype='int32')
        center.description = 'Atom index (counting from 0) of the orbital center'
        nn = ncf.createVariable('N', 'i', ('orbs',))
        nn[:] = N.array(orbbasis.N, dtype='int32')
        ll = ncf.createVariable('L', 'i', ('orbs',))
        ll[:] = N.array(orbbasis.L, dtype='int32')
        mm = ncf.createVariable('M', 'i', ('orbs',))
        mm[:] = N.array(orbbasis.M, dtype='int32')
        # Cutoff radius and delta
        Rc = ncf.createVariable('Rc', 'd', ('orbs',))
        Rc[:] = orbbasis.coff
        Rc.units = 'Angstrom'
        delta = ncf.createVariable('delta', 'd', ('orbs',))
        delta[:] = orbbasis.delta
        delta.units = 'Angstrom'
        # Radial components of the orbitals
        ntb = len(orbbasis.orb[0])
        ncf.createDimension('ntb', ntb)
        rii = ncf.createVariable('rii', 'd', ('orbs', 'ntb'))
        rii[:] = N.outer(orbbasis.delta, N.arange(ntb))
        rii.units = 'Angstrom'
        radialfct = ncf.createVariable('radialfct', 'd', ('orbs', 'ntb'))
        radialfct[:] = orbbasis.orb
        # Sort eigenvalues to connect crossing bands?
        if options.sorting:
            for i in range(SCDM.nspin):
                evals[:, i, :] = SortBands(evals[:, i, :])
        # Produce nice plots if labels exist
        if klabels:
            if SCDM.nspin == 1:
                PlotElectronBands(options.DestDir+'/Electrons.agr', dk, evals[:, 0, :], kticks)
            elif SCDM.nspin == 2:
                PlotElectronBands(options.DestDir+'/Electrons.UP.agr', dk, evals[:, 0, :], kticks)
                PlotElectronBands(options.DestDir+'/Electrons.DOWN.agr', dk, evals[:, 1, :], kticks)
        ncf.close()

    if TSrun: # only electronic calculation
        # Ugly hack to get my old code to work again. -Magnus
        if options.FermiSurface == True:
            from . import BandStruct as BS
            options.fdfFile = 'RUN.fdf'
            options.eMin, options.eMax = -10, 10
            options.NNk = 101
            BS.general = options
            BS.main()

        return SCDM.Sym.path

    # Compute phonon eigenvalues
    if options.qfile:
        SCDM.SymmetrizeFC(options.radius)
        SCDM.SetMasses()
        qpts, dq, qlabels, qticks = ReadKpoints(options.qfile)
        if qlabels:
            # Only write ascii if labels exist
            WriteKpoints(options.DestDir+'/qpoints', qpts, qlabels)
        # Prepare netcdf
        ncfn = options.DestDir+'/Phonons.nc'
        ncf = NC4.Dataset(ncfn, 'w')
        # Grid
        ncf.createDimension('gridpts', len(qpts))
        ncf.createDimension('vector', 3)
        grid = ncf.createVariable('grid', 'd', ('gridpts', 'vector'))
        grid[:] = qpts
        grid.units = '1/Angstrom'
        # Geometry
        ncf.createDimension('atoms', SCDM.Sym.basis.NN)
        xyz = ncf.createVariable('xyz', 'd', ('atoms', 'vector'))
        xyz[:] = SCDM.Sym.basis.xyz
        xyz.units = 'Angstrom'
        pbc = ncf.createVariable('pbc', 'd', ('vector', 'vector'))
        pbc.units = 'Angstrom'
        pbc[:] = [SCDM.Sym.a1, SCDM.Sym.a2, SCDM.Sym.a3]
        rvec1 = ncf.createVariable('rvec', 'd', ('vector', 'vector'))
        rvec1.units = '1/Angstrom (incl. factor 2pi)'
        rvec1[:] = rvec
        ncf.sync()
        # Loop over q
        for i, q in enumerate(qpts):
            if i < 100: # Print only for the first 100 points
                hw, U = SCDM.ComputePhononModes_q(q, verbose=True)
            else:
                hw, U = SCDM.ComputePhononModes_q(q, verbose=False)
                # otherwise something simple
                if i%100 == 0: print('%i out of %i q-points computed'%(i, len(qpts)))
            if i == 0:
                ncf.createDimension('bands', len(hw))
                ncf.createDimension('displ', len(hw))
                evals = ncf.createVariable('eigenvalues', 'd', ('gridpts', 'bands'))
                evals.units = 'eV'
                evecsRe = ncf.createVariable('eigenvectors.re', 'd', ('gridpts', 'bands', 'displ'))
                evecsIm = ncf.createVariable('eigenvectors.im', 'd', ('gridpts', 'bands', 'displ'))
                # Check eigenvectors
                print('SupercellPhonons.Checking eigenvectors at', q)
                tmp = MM.mm(N.conjugate(U), SCDM.FCtilde, N.transpose(U))
                const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5
                hw2 = const*N.diagonal(tmp)**0.5 # Units in eV
                print(' ... Allclose=', N.allclose(hw, N.absolute(hw2), atol=1e-5, rtol=1e-3))
                ncf.sync()
                # Write only AXSF files for the first q-point
                PH.WriteAXSFFiles(options.DestDir+'/q%i_re.axsf'%i, SCDM.Sym.basis.xyz, SCDM.Sym.basis.anr, hw, U.real, 1, SCDM.Sym.basis.NN)
                PH.WriteAXSFFiles(options.DestDir+'/q%i_im.axsf'%i, SCDM.Sym.basis.xyz, SCDM.Sym.basis.anr, hw, U.imag, 1, SCDM.Sym.basis.NN)
                PH.WriteFreqFile(options.DestDir+'/q%i.freq'%i, hw)
            evals[i] = hw
            evecsRe[i] = U.real
            evecsIm[i] = U.imag
        ncf.sync()
        # Sort eigenvalues to connect crossing bands?
        if options.sorting:
            evals = SortBands(evals)
        # Produce nice plots if labels exist
        if qlabels:
            PlotPhononBands(options.DestDir+'/Phonons.agr', dq, N.array(evals[:]), qticks)
        ncf.close()

    # Compute e-ph couplings
    if options.kfile and options.qfile:
        SCDM.ReadGradients()
        ncf = NC4.Dataset(options.DestDir+'/EPH.nc', 'w')
        ncf.createDimension('kpts', len(kpts))
        ncf.createDimension('qpts', len(qpts))
        ncf.createDimension('modes', len(hw))
        ncf.createDimension('nspin', SCDM.nspin)
        ncf.createDimension('bands', SCDM.rednao)
        ncf.createDimension('vector', 3)
        kgrid = ncf.createVariable('kpts', 'd', ('kpts', 'vector'))
        kgrid[:] = kpts
        qgrid = ncf.createVariable('qpts', 'd', ('qpts', 'vector'))
        qgrid[:] = qpts
        evalfkq = ncf.createVariable('evalfkq', 'd', ('kpts', 'qpts', 'nspin', 'bands'))
        # First (second) band index n (n') is the initial (final) state, i.e.,
        # Mkq(k,q,mode,spin,n,n') := < n',k+q | dV_q(mode) | n,k >
        MkqAbs = ncf.createVariable('Mkqabs', 'd', ('kpts', 'qpts', 'modes', 'nspin', 'bands', 'bands'))
        GkqAbs = ncf.createVariable('Gkqabs', 'd', ('kpts', 'qpts', 'modes', 'nspin', 'bands', 'bands'))
        ncf.sync()
        # Loop over k-points
        for i, k in enumerate(kpts):
            kpts[i] = k
            # Compute initial electronic states
            evi, eveci = SCDM.ComputeElectronStates(k, verbose=True)
            # Loop over q-points
            for j, q in enumerate(qpts):
                # Compute phonon modes
                hw, U = SCDM.ComputePhononModes_q(q, verbose=True)
                # Compute final electronic states
                evf, evecf = SCDM.ComputeElectronStates(k+q, verbose=True)
                evalfkq[i, j, :] = evf
                # Compute electron-phonon couplings
                m, g = SCDM.ComputeEPHcouplings_kq(k, q) # (modes,nspin,bands,bands)
                # Data to file
                # M (modes,spin,i,l) = m(modes,k,j) init(i,j) final(k,l)
                #                            0 1 2       0,1        0 1
                #                                ^-------^
                #                              ^----------------------^
                for ispin in range(SCDM.nspin):
                    evecfd = MM.dagger(evecf[ispin]) # (bands,bands)
                    M = N.tensordot(N.tensordot(m[:, ispin], eveci[ispin], axes=[2, 0]), evecfd, axes=[1, 1])
                    G = N.tensordot(N.tensordot(g[:, ispin], eveci[ispin], axes=[2, 0]), evecfd, axes=[1, 1])
                    MkqAbs[i, j, :, ispin] = N.absolute(M)
                    GkqAbs[i, j, :, ispin] = N.absolute(G)
                ncf.sync()
        ncf.close()
    return SCDM.Sym.path
