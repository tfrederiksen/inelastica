"""

:mod:`Inelastica.Phonons`
=========================

Overview
--------

* The dynamic atoms no longer need to be a complete range (FCfirst,FClast+1).
  Any arbitrary list of atoms (the variable ``options.DynamicAtoms``) can now be specified.

* The displacement amplitude in the various `FCrun` and `OSrun` directories
  may correspond to different values.

* The code accepts that some atoms may have been displaced in several `FCrun`
  directories. Only the first instance (first `FCrun` directory) encountered
  is read/used.

* The auxiliary NetCDF file has been eliminated by simply interchanging the
  loops over gradients and phonon modes. Thus, only one gradient need to be
  present in the memory at one time.

Thomas Frederiksen, August 2014.

Additional improvements to facilitate large-scale calcuations:

* The code allows to specify the range of dynamic atoms for which the
  electron-phonon coupling elements are evaluated
  (flags: ``--EPHfirst``, ``--EPHlast``)

* It is possible to restart the code (flag: ``--Restart``) by using the
  NetCDF output obtained from a previous calculation as checkpoint file
  (flag: ``--CheckPointNetCDF``)

* The matrix of electron-phonon coupling elements can be evaluated in single
  precision to reduce disk space usage (flag: ``--SinglePrec``).

Daniele Stradi, April 2015.

Scripting
---------

To run phonon calculations as a script one can execute:

.. code-block:: bash

    >>> import Inelatica.Phonons as P
    >>> my_argv = '--FCfirst=9 --FClast=9 --DeviceFirst=8 --DeviceLast=13 -c PHrun'
    >>> my_opts = P.GetOptions(my_argv)
    >>> P.main(my_opts)


Classes and methods
-------------------

.. autosummary::
   :toctree:

   GetOptions
   main
   FCrun
   OTSrun
   OSrun
   DynamicalMatrix

.. currentmodule:: Inelastica.Phonons

"""
from __future__ import print_function

import netCDF4 as NC4
import numpy as N
import numpy.linalg as LA
import warnings
import glob
import os
import sys
import ast
import Inelastica.io.siesta as SIO
import Inelastica.io.log as Log
import Inelastica.MakeGeom as MG
import Inelastica.physics.constants as PC
import Inelastica.math as MM
import Inelastica.misc.valuecheck as VC


def GetOptions(argv):
    """
    Returns an instance of ``options`` for the ``Phonons`` module

    Parameters
    ----------
    argv : string
        For example `-c test_dir`, which gives instructions to compute not only
        vibrational modes and frequencies, but also the corresponding electron-vibration
        couplings and to place the results in the output directory `test_dir`.
    """
    # if text string is specified, convert to list
    if isinstance(argv, VC.string_types):
        argv = argv.split()

    import argparse

    p = argparse.ArgumentParser(description='Methods to calculate vibrations and e-ph couplings from SIESTA output')
    p.add_argument('DestDir', help='Destination directory')
    p.add_argument('-f', '--fdf', dest='fdf', default='./RUN.fdf', type=str,
                 help='Input fdf-file for SIESTA/TranSIESTA calculation [%(default)s]')
    p.add_argument('-c', '--CalcCoupl', dest='CalcCoupl', action='store_true', default=False,
                   help='Calculate e-ph couplings [default=%(default)s]')
    p.add_argument('-r', '--Restart', dest='Restart', action='store_true', default=False,
                   help='Restart from a previous run [default=%(default)s]')
    p.add_argument('--CheckPointNetCDF', dest='CheckPointNetCDF', type=str, default='None',
                   help='Old NetCDF file used for restart [default=%(default)s]')
    p.add_argument('-s', '--SinglePrec', dest='SinglePrec', action='store_true', default=False,
                   help='Calculate e-ph couplings using single precision arrays [default=%(default)s]')
    p.add_argument('-F', '--DeviceFirst', dest='DeviceFirst', type=int, default=1,
                   help='First device atom index (in the electronic basis) [default=%(default)s]')
    p.add_argument('-L', '--DeviceLast', dest='DeviceLast', type=int, default=0,
                   help='Last device atom index (in the electronic basis) [default=NumberOfAtoms]')
    p.add_argument('--FCfirst', dest='FCfirst', type=int, default=1,
                   help='First FC atom index [default=%(default)s]')
    p.add_argument('--FClast', dest='FClast', type=int, default=0,
                   help='Last FC atom index [default=NumberOfAtoms]')
    p.add_argument('--EPHfirst', dest='EPHfirst', type=int, default=1,
                   help='First atom index for which the e-ph. couplings are evaluated [default=FCfirst]')
    p.add_argument('--EPHlast', dest='EPHlast', type=int, default=0,
                   help='Last atom index for which the e-ph. couplings are evaluated [default=FClast]')
    p.add_argument('--PBCFirst', dest='PBCFirst', type=int, default=1,
                   help='For eliminating interactions through periodic boundary conditions in z-direction [default=%(default)s]')
    p.add_argument('--PBCLast', dest='PBCLast', type=int, default=0,
                   help='For eliminating interactions through periodic boundary conditions in z-direction [default=NumberOfAtoms]')
    p.add_argument('--FCwildcard', dest='FCwildcard', type=str, default='./FC*',
                   help='Wildcard for FC directories [default=%(default)s]')
    p.add_argument('--OSdir', dest='onlySdir', type=str, default='./OSrun',
                   help='Location of OnlyS directory [default=%(default)s]')
    p.add_argument('-a', '--AbsoluteEnergyReference', dest='AbsEref', action='store_true', default=False,
                   help='Use an absolute energy reference (Fermi energy of equilibrium structure) for displaced Hamiltonians (e.g., when eF is not well-defined) instead of the instantaneous Fermi energy for the displaced geometries, cf. Eq.(17) in PRB 75, 205413 (2007) [default=%(default)s]')
    p.add_argument('-i', '--Isotopes', dest='Isotopes', default='[]',
                   help='String, formatted as a list [[i1,m1],...], where the mass of atom index i1 (SIESTA numbering) will be set to m1. Alternatively, the argument can be a file with the string [default=%(default)s]')
    p.add_argument('-x', '--k1', dest='k1', default=0.0, type=float,
                   help='k-point along a1 where e-ph couplings are evaluated [%(default)s]')
    p.add_argument('-y', '--k2', dest='k2', default=0.0, type=float,
                   help='k-point along a2 where e-ph couplings are evaluated [%(default)s]')
    p.add_argument('-z', '--k3', dest='k3', default=0.0, type=float,
                   help='k-point along a3 where e-ph couplings are evaluated [%(default)s]')
    p.add_argument('-g', '--WriteGradients', dest='WriteGradients', action='store_true', default=False,
                   help='Write real-space gradients dH/dR to NetCDF [default=%(default)s]')

    options = p.parse_args(argv)

    # Set module name
    options.module = 'Phonons'

    # k-point
    options.kpoint = N.array([options.k1, options.k2, options.k3], N.float)
    del options.k1, options.k2, options.k3

    # Determine array type for H,S,dH,...
    options.GammaPoint = N.dot(options.kpoint, options.kpoint) < 1e-7
    if options.GammaPoint:
        if options.SinglePrec:
            options.atype = N.float32
        else:
            options.atype = N.float64
    else:
        if options.SinglePrec:
            options.atype = N.complex64
        else:
            options.atype = N.complex128

    # Check if we need to set the last atom(s)
    if options.FClast*options.DeviceLast*options.EPHlast*options.PBCLast == 0:
        # We need NumberOfAtoms
        fdf = glob.glob(options.FCwildcard+'/'+options.fdf)
        natoms = SIO.GetFDFlineWithDefault(fdf[0], 'NumberOfAtoms', int, 1000, 'Phonons')
    if options.FClast == 0:
        options.FClast = natoms
    if options.DeviceLast == 0:
        options.DeviceLast = natoms
    if options.EPHlast == 0:
        options.EPHlast = natoms
    if options.PBCLast == 0:
        options.PBCLast = natoms

    # Dynamic atoms
    options.DynamicAtoms = list(range(options.FCfirst, options.FClast+1))

    # EPH atoms - set only different from options.DynamicAtoms if a subset is specified
    if options.EPHfirst >= options.FCfirst and options.EPHlast <= options.FClast:
        options.EPHAtoms = list(range(options.EPHfirst, options.EPHlast+1))
    else:
        options.EPHAtoms = options.DynamicAtoms
    del options.EPHfirst, options.EPHlast
    del options.FCfirst, options.FClast

    # PBCFirst/PBCLast
    if options.PBCFirst < options.DeviceFirst:
        options.PBCFirst = options.DeviceFirst
    if options.PBCLast > options.DeviceLast:
        options.PBCLast = options.DeviceLast

    # Isotopes specified in separate file?
    if os.path.isfile(options.Isotopes):
        f = open(options.Isotopes)
        s = ''
        for line in f.readlines():
            s += line.replace('\n', '')
        options.Isotopes = s
    options.Isotopes = ast.literal_eval(options.Isotopes)

    return options


class FCrun(object):

    def __init__(self, runfdf):
        self.fdf = runfdf
        self.directory, self.tail = os.path.split(runfdf)
        self.systemlabel = SIO.GetFDFlineWithDefault(runfdf, 'SystemLabel', str, 'siesta', 'Phonons')
        FCfirst = SIO.GetFDFlineWithDefault(runfdf, 'MD.FCfirst', int, 0, 'Phonons')
        FClast = SIO.GetFDFlineWithDefault(runfdf, 'MD.FClast', int, 0, 'Phonons')
        # Finite-displacement amplitude
        ampl, unit = SIO.GetFDFline(runfdf, KeyWord='MD.FCDispl')
        if unit.upper() == 'ANG':
            self.Displ = float(ampl)
        elif unit.upper() == 'BOHR':
            self.Displ = float(ampl)*PC.Bohr2Ang
        print('Displacement = %.6f Ang'%self.Displ)
        # Read geometry
        self.geom = MG.Geom(runfdf)
        # Compare with XV file corrected for last displacement
        XV = self.directory+'/%s.XV'%self.systemlabel
        geomXV = MG.Geom(XV)
        geomXV.xyz[FClast-1, 2] -= self.Displ
        if not N.allclose(geomXV.xyz, self.geom.xyz):
            sys.exit('Error: Geometries %s and %s should differ ONLY by displacement of atom %s in z'\
                     %(runfdf, XV, FClast))
        # Set up FC[i,a,j,b]: Force constant (eV/A^2) from moved atom i, axis a to atom j, axis b
        natoms = self.geom.natoms
        self.m = N.zeros((FClast-FCfirst+1, 3, natoms, 3), N.float)
        self.p = N.zeros((FClast-FCfirst+1, 3, natoms, 3), N.float)
        fc = N.array(SIO.ReadFCFile(self.directory+'/%s.FC'%self.systemlabel))
        for i in range(FClast-FCfirst+1):
            for j in range(3):
                self.m[i, j] = fc[2*(3*i+j)*natoms:(2*(3*i+j)+1)*natoms]
                self.p[i, j] = fc[(2*(3*i+j)+1)*natoms:(2*(3*i+j)+2)*natoms]
        # Correct force constants for the moved atom
        # Cf. Eq. (13) in Frederiksen et al. PRB 75, 205413 (2007)
        for i in range(FClast-FCfirst+1):
            for j in range(3):
                self.m[i, j, FCfirst-1+i, :] = 0.0
                self.m[i, j, FCfirst-1+i, :] = -N.sum(self.m[i, j], axis=0)
                self.p[i, j, FCfirst-1+i, :] = 0.0
                self.p[i, j, FCfirst-1+i, :] = -N.sum(self.p[i, j], axis=0)
        self.DynamicAtoms = list(range(FCfirst, FClast+1))

        # Determine TSHS files
        files = glob.glob(self.directory+'/%s*.TSHS'%self.systemlabel)
        files.sort()
        if self.directory+'/%s.TSHS'%self.systemlabel in files:
            # If present, ignore this file which does not origninate from the FCrun
            files.remove(self.directory+'/%s.TSHS'%self.systemlabel)
        if (FClast-FCfirst+1)*6+1 != len(files):
            warnings.warn('Phonons.GetFileLists: WARNING - Inconsistent number of *.TSHS files in %s'%self.directory)
            return

        # Build dictionary over TSHS files and corresponding displacement amplitudes
        self.TSHS = {}
        self.TSHS[0] = files[0] # Equilibrium TSHS
        for i, v in enumerate(self.DynamicAtoms):
            for j in range(3):
                # Shifted TSHS files (atom,axis,direction)
                self.TSHS[v, j, -1] = files[1+6*i+2*j]
                self.TSHS[v, j, 1] = files[1+6*i+2*j+1]

    def GetOrbitalIndices(self):
        # Determine snr (siesta number) for each label
        csl = SIO.GetFDFblock(self.fdf, KeyWord='ChemicalSpeciesLabel')
        csl2snr = {}
        for set in csl:
            csl2snr[set[2]] = set[0]
        # Determine nao (number of orbitals) for each snr
        ionNCfiles = glob.glob(self.directory+'/*.ion.nc*')
        snr2nao = {}
        for ionfile in ionNCfiles:
            if ionfile.endswith('.gz'):
                print('Phonons.GetOrbitalIndices: Unzipping: ' + ionfile)
                os.system('gunzip '+ionfile)
                ionfile = ionfile[:-3]
            file = NC4.Dataset(ionfile, 'r')
            thissnr = int(csl2snr[file.Label])
            snr2nao[thissnr] = file.Number_of_orbitals
            file.close()
        print('Phonons.GetOrbitalIndices: Dictionary snr2nao = ' + str(snr2nao))
        # Determine which orbital indices that belongs to a certain atom
        orbitalIndices = []
        tmpOrb = 0
        for num in self.geom.snr:
            nao = snr2nao[num]
            orbitalIndices.append([tmpOrb, tmpOrb+int(nao)-1])
            tmpOrb += nao
        self.orbitalIndices = N.array(orbitalIndices)
        self.nao = tmpOrb # total number of orbitals
        self.snr2nao = snr2nao # orbitals per species
        return self.orbitalIndices, self.nao


class OTSrun(FCrun): # Only TranSiesta run

    def __init__(self, runfdf):
        self.fdf = runfdf
        self.directory, self.tail = os.path.split(runfdf)
        self.systemlabel = SIO.GetFDFlineWithDefault(runfdf, 'SystemLabel', str, 'siesta', 'Phonons')
        # Read geometry
        self.geom = MG.Geom(runfdf)
        # Compare with XV file corrected for last displacement
        XV = self.directory+'/%s.XV'%self.systemlabel
        # Determine TSHS files
        files = glob.glob(self.directory+'/%s*.TSHS'%self.systemlabel)
        # Build dictionary over TSHS files and corresponding displacement amplitudes
        self.TSHS = {}
        try:
            self.TSHS[0] = files[0] # Equilibrium TSHS
        except:
            warnings.warn('Phonons.GetFileLists: No TSHS file found in %s'%self.directory)


class OSrun(object):

    def __init__(self, onlySdir, kpoint, atype=N.complex):
        print('Phonons.GetOnlyS: Reading from: ' + onlySdir)
        onlySfiles = glob.glob(onlySdir+'/*.onlyS*')
        onlySfiles.sort()
        if len(onlySfiles) < 1:
            sys.exit('Phonons.GetOnlyS: No .onlyS file found!')
        if len(onlySfiles) != 6:
            sys.exit('Phonons.GetOnlyS: Wrong number of onlyS files found!')
        else:
            onlyS = {}
            Displ = {}
            for file in onlySfiles:
                thisHS = SIO.HS(file)
                thisHS.setkpoint(kpoint, atype=atype)
                S = thisHS.S
                del thisHS
                nao = len(S) // 2
                S0 = S[0:nao, 0:nao].copy()
                dmat = S[0:nao, nao:nao*2].copy()
                if file.endswith('_1.onlyS'):
                    onlyS[0, -1] = dmat
                elif file.endswith('_2.onlyS'):
                    onlyS[0, 1] = dmat
                elif file.endswith('_3.onlyS'):
                    onlyS[1, -1] = dmat
                elif file.endswith('_4.onlyS'):
                    onlyS[1, 1] = dmat
                elif file.endswith('_5.onlyS'):
                    onlyS[2, -1] = dmat
                elif file.endswith('_6.onlyS'):
                    onlyS[2, 1] = dmat
            # Loop over the 6 doubled geometries and determine the displacement
            for i in range(1, 7):
                thisd = 1e10
                xyz = N.array(SIO.Getxyz(onlySdir+'/RUN_%i.fdf'%i))
                #for j in range(1, len(xyz)):
                #    thisd = min(thisd, (N.dot(xyz[0]-xyz[j], xyz[0]-xyz[j]))**.5)
                j = len(xyz) // 2
                v = xyz[0]-xyz[j]
                thisd = N.dot(v, v)**.5
                Displ[(i-1) // 2, 1-2*(i % 2)] = thisd
                print('Phonons.GetOnlyS: OnlyS-displacement (min) = %.5f Ang'%thisd)
            # Construct dS array
            self.S0 = S0
            self.dS = N.empty((3,)+dmat.shape, dtype=dmat.dtype)
            for j in range(3):
                self.dS[j] = (onlyS[j, 1]-onlyS[j, -1])/(Displ[j, -1]+Displ[j, 1])
            self.Displ = Displ


class DynamicalMatrix(object):

    def __init__(self, fdfs, DynamicAtoms=None, TSrun=False):
        self.fdfs = fdfs
        if TSrun:
            self.FCRs = [OTSrun(f) for f in fdfs]
            self.TSHS = {}
            self.TSHS[0] = self.FCRs[0].TSHS[0]
        else:
            self.FCRs = [FCrun(f) for f in fdfs]
        self.geom = self.FCRs[0].geom # assume identical geometries
        if DynamicAtoms:
            self.SetDynamicAtoms(DynamicAtoms)

    def SetDynamicAtoms(self, DynamicAtoms):
        self.DynamicAtoms = DynamicAtoms
        NN = len(DynamicAtoms)
        self.m = N.zeros((NN, 3, self.geom.natoms, 3), N.complex)
        self.p = N.zeros((NN, 3, self.geom.natoms, 3), N.complex)
        self.Displ = {}
        self.TSHS = {}
        try:
            self.TSHS[0] = self.FCRs[0].TSHS[0]
            has_TSHS = True
        except:
            has_TSHS = False
        for i, v in enumerate(DynamicAtoms):
            for fcr in self.FCRs:
                if not v in fcr.DynamicAtoms:
                    continue

                j = fcr.DynamicAtoms.index(v)
                print('Reading FC data for dynamic atom %i from %s' %(v, fcr.fdf))
                self.Displ[v] = fcr.Displ
                self.m[i] = fcr.m[j]
                self.p[i] = fcr.p[j]
                if has_TSHS:
                    for k in range(3):
                        self.TSHS[v, k, -1] = fcr.TSHS[v, k, -1]
                        self.TSHS[v, k, 1] = fcr.TSHS[v, k, 1]
                break
            # Check that we really found the required atom
            if len(self.Displ) <= i:
                sys.exit('Error: Did not find FC data for a dynamic atom %i'%v)
        self.mean = (self.m+self.p)/2

    def SetMasses(self, Isotopes=[]):
        self.Masses = []
        # Set default masses
        for i, v in enumerate(self.DynamicAtoms):
            try:
                self.Masses.append(PC.AtomicMass[self.geom.anr[v-1]])
            except:
                print('WARNING: Mass of atom %i unknown, set arbitrarily to 1000.00 amu'%v)
                self.Masses.append(1000.00)
        # Override with specified masses?
        for ii, mass in Isotopes:
            if ii in self.DynamicAtoms:
                j = self.DynamicAtoms.index(ii)
                print('Phonons.Analyse: Setting mass for atom %i (SIESTA numbering) to %f:'%(ii, mass))
                self.Masses[j] = mass

    def ApplySumRule(self, FC):
        FC0 = FC.copy()
        # Correct force constants for the moved atom
        # Cf. Eq. (13) in Frederiksen et al. PRB 75, 205413 (2007)
        for i, v in enumerate(self.DynamicAtoms):
            for j in range(3):
                FC[i, j, v-1, :] = 0.0
                FC[i, j, v-1, :] = -N.sum(FC[i, j], axis=0)
        print('Total sumrule change in FC: %.3e eV/Ang^2' % N.sum(abs(FC0)-abs(FC)))
        return FC

    def ComputePhononModes(self, FC, verbose=True):
        dyn = len(self.DynamicAtoms)
        FCtilde = N.zeros((dyn, 3, dyn, 3), N.complex)
        # Symmetrize and mass-scale
        for i, v in enumerate(self.DynamicAtoms):
            for j, w in enumerate(self.DynamicAtoms):
                FCtilde[i, :, j, :] = 0.5*(FC[i, :, w-1, :]+MM.dagger(FC[j, :, v-1, :]))\
                                      /(self.Masses[i]*self.Masses[j])**0.5
        # Solve eigenvalue problem with symmetric FCtilde
        FCtilde = FCtilde.reshape((3*dyn, 3*dyn), order='C')
        self.FCtilde = FCtilde
        evalue, evec = LA.eigh(FCtilde)
        #evalue,evec = LA.eig(FCtilde)
        evec = N.transpose(evec)
        evalue = N.array(evalue, N.complex)
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
            U[i] = U[i]/(N.dot(N.conjugate(U[i]), U[i])**0.5)
        # Sort in order descending mode energies
        #hw = hw[::-1] # reverse array
        #U = U[::-1] # reverse array
        indx = N.argsort(hw)[::-1] # reverse
        hw = hw[indx]
        U = U[indx]
        # Print mode frequencies
        if verbose:
            print('Phonons.CalcPhonons: Frequencies in meV:')
            for i in range(3*dyn):
                print(('%.3f'%(1000*hw[i])).rjust(9), end='')
                if (i-5)%6 == 0: print()
            if (i-5)%6 != 0: print()
        #print 'Phonons.CalcPhonons: Frequencies in cm^-1:'
        #for i in range(3*dyn):
        #    print ('%.3f'%(hw[i]/PC.invcm2eV)).rjust(9),
        #    if (i-5)%6 == 0: print
        #if (i-5)%6 != 0: print

        # Compute real displacement vectors
        Udisp = U.copy()
        for i in range(3*dyn):
            # Eigenvectors after division by sqrt(mass)
            Udisp[:, i] = U[:, i] / self.Masses[i // 3]**.5

        # Compute displacement vectors scaled for the characteristic length
        Ucl = N.zeros_like(U)
        Ucl[hw > 0, :] = PC.hbar2SI / N.sqrt(
            N.array(self.Masses).repeat(3).reshape(1, -1) * PC.amu2kg
            * hw[hw > 0].reshape(-1, 1) * PC.eV2Joule
            ) * 1e10 * U[hw > 0, :]
        # Note that if we displace by the characteristic length
        # via E=1/2 <u|FC|u>, the energy should change by the
        # characteristic energy (which is hw/2), i.e.,
        # np.diag(Ucl.dot(fcmat).dot(Ucl.T)).real[hw > 0]
        # should be identical to hw

        # Expand vectors to full geometry
        UU = N.zeros((len(hw), self.geom.natoms, 3), N.complex)
        UUdisp = N.zeros((len(hw), self.geom.natoms, 3), N.complex)
        UUcl = N.zeros((len(hw), self.geom.natoms, 3), N.complex)
        for i in range(len(hw)):
            for j, v in enumerate(self.DynamicAtoms):
                UU[i, v-1, :] = [U[i, 3*j], U[i, 3*j+1], U[i, 3*j+2]]
                UUdisp[i, v-1, :] = [Udisp[i, 3*j], Udisp[i, 3*j+1], Udisp[i, 3*j+2]]
                UUcl[i, v-1, :] = [Ucl[i, 3*j], Ucl[i, 3*j+1], Ucl[i, 3*j+2]]
        self.hw = hw
        self.U = U
        self.Udisp = Udisp
        self.Ucl = Ucl
        self.UU = UU
        self.UUdisp = UUdisp
        self.UUcl = UUcl

    def PrepareGradients(self, onlySdir, kpoint, DeviceFirst, DeviceLast, AbsEref, atype, TSrun=False):
        print('\nPhonons.PrepareGradients: Setting up various arrays')
        self.atype = atype
        self.kpoint = kpoint
        self.OrbIndx, nao = self.FCRs[0].GetOrbitalIndices()
        self.TSHS0 = SIO.HS(self.TSHS[0])
        self.TSHS0.setkpoint(kpoint, atype=atype)
        if not TSrun:
            OS = OSrun(onlySdir, kpoint, atype=atype)
            self.dS = OS.dS
            # OS.S0 and TSHS0.S should be identical, but with some versions/compilations
            # of TranSIESTA this is NOT always the case away from k=0 (GammaPoint is OK).
            # It appears to be a bug in TranSIESTA 3.2 and 4.0b affecting runs
            # with TS.onlyS=True, i.e., the quick evaluations in the OSrun folder
            if not N.allclose(OS.S0, self.TSHS0.S):
                sys.exit('Inconsistency detected with your .onlyS files. Perhaps a bug in your TranSIESTA version/compilation or overlaps beyond first neighbor cells.')
        self.invS0H0 = N.empty((2,)+self.TSHS0.H.shape, dtype=self.TSHS0.H.dtype)
        #invS0 = LA.inv(OS.S0) # <--- This choice was used in rev. 324-397
        invS0 = LA.inv(self.TSHS0.S) # Reverting to the matrix used up to rev. 323
        self.nspin = len(self.TSHS0.H)
        for iSpin in range(self.nspin):
            self.invS0H0[0, iSpin, :, :] = MM.mm(invS0, self.TSHS0.H[iSpin, :, :])
            self.invS0H0[1, iSpin, :, :] = MM.mm(self.TSHS0.H[iSpin, :, :], invS0)
        del invS0
        # don't re-create the array every time... too expensive
        self.dSdij = N.zeros((nao, nao), atype)

        # Take Device region
        self.DeviceAtoms = list(range(DeviceFirst, DeviceLast+1))
        first, last = self.OrbIndx[DeviceFirst-1][0], self.OrbIndx[DeviceLast-1][1]
        self.h0 = self.TSHS0.H[:, first:last+1, first:last+1]
        self.s0 = self.TSHS0.S[first:last+1, first:last+1]
        self.DeviceFirst = DeviceFirst
        self.DeviceLast = DeviceLast
        self.AbsEref = AbsEref

    def GetGradient(self, Atom, Axis):
        print('\nPhonons.GetGradient: Computing dH[%i,%i]'%(Atom, Axis))
        # Read TSHS files
        TSHSm = SIO.HS(self.TSHS[Atom, Axis, -1])
        TSHSm.setkpoint(self.kpoint, atype=self.atype)
        TSHSp = SIO.HS(self.TSHS[Atom, Axis, 1])
        TSHSp.setkpoint(self.kpoint, atype=self.atype)
        # Use Fermi energy of equilibrium calculation as energy reference?
        if self.AbsEref:
            print('Computing gradient with absolute energy reference')
            for iSpin in range(self.nspin):
                TSHSm.H[iSpin, :, :] += (TSHSm.ef-self.TSHS0.ef)*TSHSm.S
                TSHSp.H[iSpin, :, :] += (TSHSp.ef-self.TSHS0.ef)*TSHSp.S
        # Compute direct gradient
        dH = (TSHSp.H-TSHSm.H)/(2*self.Displ[Atom])
        del TSHSm, TSHSp
        # Orbital range for the displaced atom:
        f, l = self.OrbIndx[Atom-1]
        self.dSdij[:, f:l+1] = self.dS[Axis, :, f:l+1]
        # Apply Pulay-type corrections
        for iSpin in range(self.nspin):
            dH[iSpin, :, :] -= MM.mm(MM.dagger(self.dSdij), self.invS0H0[0, iSpin, :, :]) \
                             + MM.mm(self.invS0H0[1, iSpin, :, :], self.dSdij)
        self.dSdij[:, f:l+1] = 0. # reset
        return dH

    def ComputeEPHcouplings(self, PBCFirst, PBCLast, EPHAtoms, Restart, CheckPointNetCDF, WriteGradients=False):
        first, last = self.OrbIndx[self.DeviceFirst-1][0], self.OrbIndx[self.DeviceLast-1][1]
        rednao = last+1-first
        const = PC.hbar2SI*(1e20/(PC.eV2Joule*PC.amu2kg))**0.5

        if Restart:
            if not os.path.exists(CheckPointNetCDF):
                print("ERROR!!! You have enforced a restart but you have not provided any checkpoint NetCDF file!")
                sys.exit()
            else:
                print("Restart e-ph. calculation. Partial information read from file: " + CheckPointNetCDF)
                RNCfile = NC4.Dataset(CheckPointNetCDF, 'r')
                Heph = N.array(RNCfile.variables['He_ph'], self.atype)
        else:
            print("Start e-ph. calculation from scratch")
            Heph = N.zeros((len(self.hw), self.nspin, rednao, rednao), self.atype)

        if WriteGradients:
            self.gradients = []
        # Loop over dynamic atoms
        for i, v in enumerate(EPHAtoms):
            # Loop over axes
            for j in range(3):
                # Compute gradient
                dH = self.GetGradient(v, j)
                # Remove Periodic Boundary terms
                nuo = len(dH[0])
                pbcf = self.OrbIndx[PBCFirst-1][0]
                pbcl = self.OrbIndx[PBCLast-1][1]
                if v < PBCFirst:
                    # we have something to remove...
                    print('Warning: Setting certain elements in dH[%i,%i] to zero because %i<PBCFirst'%(v, j, v))
                    #bb = (PBCFirst - FCfirst) * 3
                    dH[:, pbcl+1:nuo, :] = 0.0
                    dH[:, :, pbcl+1:nuo] = 0.0
                if PBCLast < v:
                    print('Warning: Setting certain elements in dH[%i,%i] to zero because PBCLast<%i'%(v, j, v))
                    #aa = (PBCLast - FCfirst) * 3
                    dH[:, :pbcf-1, :] = 0.0
                    dH[:, :, :pbcf-1] = 0.0
                # Device part
                dh = dH[:, first:last+1, first:last+1]
                if WriteGradients:
                    self.gradients.append(dh)
                # Loop over modes and throw away the gradient (to save memory)
                for m in range(len(self.hw)):
                    if self.hw[m] > 0:
                        # Eigenvectors should be real for GammaPoint phonons, hence we always take the real part
                        Heph[m] += const*dh*self.UU[m, v-1, j].real/(2*self.Masses[i]*self.hw[m])**.5
                    elif i == 0:
                        # Print only first time
                        print('Phonons.ComputeEPHcouplings: Mode %i has nonpositive frequency --> Zero-valued coupling matrix'%m)
                        # already zero
        del dH, dh
        self.heph = Heph
        if WriteGradients:
            self.gradients = N.array(self.gradients)

    def WriteOutput(self, label, SinglePrec, GammaPoint):
        print('\nPhonons.WriteOutput')
        ### Write MKL- and xyz-files
        natoms = self.geom.natoms
        hw = self.hw
        # Write only real part of eigenvectors
        UU = self.UU.reshape(len(hw), 3*natoms).real
        UUdisp = self.UUdisp.reshape(len(hw), 3*natoms).real
        UUcl = self.UUcl.reshape(len(hw), 3*natoms).real

        SIO.WriteMKLFile('%s.mkl'%label, self.geom.anr, self.geom.xyz, hw, UU, 1, natoms)
        SIO.WriteMKLFile('%s.real-displ.mkl'%label, self.geom.anr, self.geom.xyz, hw, UUdisp, 1, natoms)
        SIO.WriteXYZFile('%s.xyz'%label, self.geom.anr, self.geom.xyz)
        WriteFreqFile('%s.freq'%label, hw)
        WriteVibDOSFile('%s.Gfdos'%label, hw, type='Gaussian')
        WriteVibDOSFile('%s.Lfdos'%label, hw, type='Lorentzian')
        WriteAXSFFiles('%s.mol.axsf'%label, self.geom.xyz, self.geom.anr, hw, UU, 1, natoms)
        WriteAXSFFiles('%s.mol.real-displ.axsf'%label, self.geom.xyz, self.geom.anr, hw, UUdisp, 1, natoms)
        WriteAXSFFiles('%s.mol.charlength-displ.axsf'%label, self.geom.xyz, self.geom.anr, hw, UUcl, 1, natoms)
        WriteAXSFFilesPer('%s.per.axsf'%label, self.geom.pbc, self.geom.xyz, self.geom.anr, hw, UU, 1, natoms)
        WriteAXSFFilesPer('%s.per.real-displ.axsf'%label, self.geom.pbc, self.geom.xyz, self.geom.anr,
                          hw, UUdisp, 1, natoms)
        WriteAXSFFilesPer('%s.per.charlength-displ.axsf'%label, self.geom.pbc, self.geom.xyz, self.geom.anr,
                          hw, UUcl, 1, natoms)
        # Netcdf format
        ncdffn = '%s.nc'%label
        print('Phonons.WriteOutput: Writing ' + ncdffn)
        ncdf = NC4.Dataset(ncdffn, 'w')
        ncdf.createDimension('one', 1)
        ncdf.createDimension('xyz', 3)
        ncdf.createDimension('modes', len(hw))
        ncdf.createDimension('natoms', self.geom.natoms)
        ncdf.createDimension('dyn_atoms', len(self.DynamicAtoms))
        ncdf.createVariable('hw', 'd', ('modes',))
        ncdf.variables['hw'][:] = hw
        ncdf.variables['hw'].info = 'Phonon frequencies'
        ncdf.variables['hw'].unit = 'eV'
        ncdf.createVariable('U', 'd', ('modes', 'natoms', 'xyz'))
        ncdf.variables['U'][:] = self.UU.real
        ncdf.variables['U'].info = 'Real part of the phonon eigenvectors'
        ncdf.createVariable('Udisp', 'd', ('modes', 'natoms', 'xyz'))
        ncdf.variables['Udisp'][:] = self.UUdisp.real
        ncdf.variables['Udisp'].info = 'Real part of the phonon displacement vectors'
        ncdf.createVariable('Ucl', 'd', ('modes', 'natoms', 'xyz'))
        ncdf.variables['Ucl'][:] = self.UUcl.real
        ncdf.variables['Ucl'].info = 'Real part of displacement vectors scaled for the characteristic length'
        ncdf.createVariable('CellVectors', 'd', ('xyz', 'xyz'))
        ncdf.variables['CellVectors'][:] = self.geom.pbc
        ncdf.variables['CellVectors'].info = 'Unit cell vectors'
        ncdf.variables['CellVectors'].unit = 'Ang'
        ncdf.createVariable('GeometryXYZ', 'd', ('natoms', 'xyz'))
        ncdf.variables['GeometryXYZ'][:] = self.geom.xyz
        ncdf.variables['GeometryXYZ'].info = 'Atomic coordinates of all atoms in cell'
        ncdf.variables['GeometryXYZ'].unit = 'Ang'
        ncdf.createVariable('FC', 'd', ('dyn_atoms', 'xyz', 'natoms', 'xyz'))
        ncdf.variables['FC'][:] = self.mean
        ncdf.variables['FC'].info = 'Force matrix'
        ncdf.variables['FC'].unit = 'eV/Ang^2'
        ncdf.createVariable('AtomNumbers', 'i', ('natoms',))
        ncdf.variables['AtomNumbers'][:] = self.geom.anr
        ncdf.variables['AtomNumbers'].info = 'Element number for each atom (anr)'
        ncdf.createVariable('SpeciesNumbers', 'i', ('natoms',))
        ncdf.variables['SpeciesNumbers'][:] = self.geom.snr
        ncdf.variables['SpeciesNumbers'].info = 'Siesta species number (snr)'
        ncdf.createVariable('Masses', 'd', ('dyn_atoms',))
        ncdf.variables['Masses'][:] = self.Masses
        ncdf.variables['Masses'].info = 'Atomic masses of each dynamic atom'
        ncdf.variables['Masses'].unit = 'Atomic units'
        ncdf.createVariable('DynamicAtoms', 'i', ('dyn_atoms',))
        ncdf.variables['DynamicAtoms'][:] = self.DynamicAtoms
        ncdf.variables['DynamicAtoms'].info = 'Dynamic atoms (SIESTA numbering)'
        try:
            self.h0
            ncdf.createDimension('dev_atoms', len(self.DeviceAtoms))
            ncdf.createDimension('nspin', len(self.h0))
            ncdf.createDimension('norb', len(self.h0[0]))
            ncdf.createVariable('DeviceAtoms', 'i', ('dev_atoms',))
            ncdf.variables['DeviceAtoms'][:] = self.DeviceAtoms
            ncdf.variables['DeviceAtoms'].info = 'Device atoms (SIESTA numbering)'
            ncdf.createVariable('kpoint', 'd', ('xyz',))
            ncdf.variables['kpoint'][:] = self.kpoint
            ncdf.variables['kpoint'].info = 'Vector in orthogonal space (not in reciprocal space)'
            ncdf.createVariable('H0', 'd', ('nspin', 'norb', 'norb'))
            ncdf.variables['H0'][:] = self.h0.real
            ncdf.variables['H0'].info = 'Real part of electronic Hamiltonian'
            ncdf.variables['H0'].unit = 'eV'
            ncdf.createVariable('S0', 'd', ('nspin', 'norb', 'norb'))
            ncdf.variables['S0'][:] = self.s0.real
            ncdf.variables['S0'].info = 'Electronic overlap matrix'
            if not GammaPoint:
                ncdf.createVariable('H0.imag', 'd', ('nspin', 'norb', 'norb'))
                ncdf.variables['H0.imag'][:] = self.h0.imag
                ncdf.variables['H0.imag'].info = 'Imaginary part of Hamiltonian'
                ncdf.createVariable('S0.imag', 'd', ('nspin', 'norb', 'norb'))
                ncdf.variables['S0.imag'][:] = self.s0.imag
                ncdf.variables['S0.imag'].info = 'Imaginary part of overlap'
            print('Phonons.WriteOutput: Wrote H and S to ' + ncdffn)
        except:
            print('Hamiltonian etc not computed')
        # Precision for He_ph (and gradients)
        if SinglePrec:
            atype = 'f'
        else:
            atype = 'd'
        try:
            self.heph
            ncdf.createVariable('He_ph', atype, ('modes', 'nspin', 'norb', 'norb'))
            ncdf.variables['He_ph'][:] = self.heph.real
            ncdf.variables['He_ph'].info = 'Real part of EPH couplings'
            ncdf.variables['He_ph'].unit = 'eV'
            if not GammaPoint:
                ncdf.createVariable('ImHe_ph', atype, ('modes', 'nspin', 'norb', 'norb'))
                ncdf.variables['ImHe_ph'][:] = self.heph.imag
                ncdf.variables['ImHe_ph'].info = 'Imaginary part of EPH couplings'
            print('Phonons.WriteOutput: Wrote He_ph to ' + ncdffn)
        except:
            print('EPH couplings etc not computed')
        try:
            self.gradients
            ncdf.createVariable('grad.re', atype, ('modes', 'nspin', 'norb', 'norb'))
            ncdf.variables['grad.re'][:] = self.gradients.real
            ncdf.variables['grad.re'].info = 'Real part of gradients'
            ncdf.variables['grad.re'].unit = 'eV/Ang'
            if not GammaPoint:
                ncdf.createVariable('grad.im', atype, ('modes', 'nspin', 'norb', 'norb'))
                ncdf.variables['grad.im'][:] = self.gradients.imag
                ncdf.variables['grad.im'].info = 'Imaginary part of gradients'
            print('Phonons.WriteOutput: Wrote gradients to ' + ncdffn)
        except:
            print('Phonons.WriteOutpot: Gradients not computed')
        ncdf.close()
        print('Phonons.WriteOutput: Finished ' + ncdffn)


def WriteFreqFile(filename, hw):
    print('Phonons.WriteFreqFile: Writing ' + filename)
    file = open(filename, 'w')
    file.write('# ')
    for i in range(len(hw)):
        file.write(' %f '%(1000*hw[i]))
    file.write('\n# index   freq/meV \n')
    for i in range(len(hw)):
        file.write('%i  %f \n'%(i, 1000*hw[i]))
    file.close()


def WriteVibDOSFile(filename, hw, type='Gaussian'):
    'Vibrational DOS with Gaussian or Lorentzian broadening'
    fmax = max(hw)
    erng = N.linspace(0, 1.2*fmax, 1001)
    eta = N.linspace(0.001, 0.01, 11)
    ERNG = N.outer(erng, 0*eta+1.)
    ETA = N.outer(0*erng+1, eta)
    spectrum = N.zeros((len(erng), len(eta)), N.float)
    for i in range(len(hw)):
        if type == 'Gaussian':
            spectrum += (2*N.pi)**(-.5)/ETA*N.exp(N.clip(-1.0*(hw[i]-ERNG)**2/(2*ETA**2), -300, 300))
            spectrum -= (2*N.pi)**(-.5)/ETA*N.exp(N.clip(-1.0*(-hw[i]-ERNG)**2/(2*ETA**2), -300, 300))
        elif type == 'Lorentzian':
            spectrum += 1/N.pi*ETA/((hw[i]-ERNG)**2+ETA**2)
            spectrum -= 1/N.pi*ETA/((-hw[i]-ERNG)**2+ETA**2)
    # Write data to file
    print('Phonons.WriteVibDOSFile: Writing ' + filename)
    f = open(filename, 'w')
    f.write('\n# energy/eV  DOS/atom (eta=1,2,3,...,10meV) \n')
    for i in range(len(erng)):
        f.write('%.5e  '%(erng[i]))
        for j in range(len(eta)):
            f.write('%.5e '%(spectrum[i, j]*3/len(hw)))
        f.write('\n')
    f.close()


def WriteAXSFFiles(filename, xyz, anr, hw, U, FCfirst, FClast):
    'Writes the vibrational normal coordinates in xcrysden axsf-format (isolated molecule)'
    print('Phonons.WriteAXSFFile: Writing ' + filename)
    f = open(filename, 'w')
    f.write('ANIMSTEPS %i\n'%len(hw))
    for i in range(len(hw)):
        f.write('ATOMS %i\n'%(i+1))
        for j in range(len(xyz)):
            ln = ' %i'%anr[j]
            for k in range(3):
                ln += ' %.6f'%xyz[j][k]
            if j < FCfirst-1 or j > FClast-1:
                ln += ' %.6f %.6f %.6f'%(0, 0, 0)
            else:
                for k in range(3):
                    try:
                        ln += ' %.6f'%U[i][3*(j+1-FCfirst)+k]
                    except:
                        ln += ' %.6f'%U[i, j, k]
            ln += '\n'
            f.write(ln)
    f.close()


def WriteAXSFFilesPer(filename, vectors, xyz, anr, hw, U, FCfirst, FClast):
    'Writes the vibrational normal coordinates in xcrysden axsf-format (periodic structure)'
    print('Phonons.WriteAXSFFilePer: Writing ' + filename)
    VEC = N.zeros((len(hw), 3*len(xyz)), N.float)
    VEC[:, 3*(FCfirst-1):3*FClast] = U
    f = open(filename, 'w')
    f.write('ANIMSTEPS %i\nCRYSTAL\n'%len(hw))
    for i in range(len(hw)):
        f.write('PRIMVEC %i\n'%(i+1))
        f.write('%.6f %.6f %.6f\n'%(vectors[0][0], vectors[0][1], vectors[0][2]))
        f.write('%.6f %.6f %.6f\n'%(vectors[1][0], vectors[1][1], vectors[1][2]))
        f.write('%.6f %.6f %.6f\n'%(vectors[2][0], vectors[2][1], vectors[2][2]))
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
    """
    Main routine to compute vibrational modes and frequencies (and optionally
    also the corresponding electron-vibration couplings)

    Parameters
    ----------
    options : an ``options`` instance
    """
    Log.CreatePipeOutput(options)
    Log.PrintMainHeader(options)

    # Determine SIESTA input fdf files in FCruns
    fdf = glob.glob(options.FCwildcard+'/'+options.fdf)

    print('Phonons.Analyze: This run uses')
    FCfirst, FClast = min(options.DynamicAtoms), max(options.DynamicAtoms)
    print('  ... FCfirst     = %4i, FClast     = %4i, Dynamic atoms = %4i'
          %(FCfirst, FClast, len(options.DynamicAtoms)))
    print('  ... DeviceFirst = %4i, DeviceLast = %4i, Device atoms  = %4i'
          %(options.DeviceFirst, options.DeviceLast, options.DeviceLast-options.DeviceFirst+1))
    print('  ... PBC First   = %4i, PBC Last   = %4i, Device atoms  = %4i'
          %(options.PBCFirst, options.PBCLast, options.PBCLast-options.PBCFirst+1))

    print('\nSetting array type to %s\n'%options.atype)

    # Build Dynamical Matrix
    DM = DynamicalMatrix(fdf, options.DynamicAtoms)
    DM.SetMasses(options.Isotopes)
    # Compute modes
    DM.ComputePhononModes(DM.mean)
    # Compute e-ph coupling
    if options.CalcCoupl:
        DM.PrepareGradients(options.onlySdir, options.kpoint, options.DeviceFirst, options.DeviceLast, options.AbsEref, options.atype)
        DM.ComputeEPHcouplings(options.PBCFirst, options.PBCLast, options.EPHAtoms, options.Restart, options.CheckPointNetCDF,
                               WriteGradients=options.WriteGradients)
        # Write data to files
        DM.WriteOutput(options.DestDir+'/Output', options.SinglePrec, options.GammaPoint)
        Log.PrintMainFooter(options)
        return DM.h0, DM.s0, DM.hw, DM.heph
    else:
        DM.WriteOutput(options.DestDir+'/Output', options.SinglePrec, options.GammaPoint)
        Log.PrintMainFooter(options)
        return DM.hw
