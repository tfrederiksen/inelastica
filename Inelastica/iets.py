"""

:mod:`Inelastica.iets`
======================

Inelastic transport module.

Scripting
---------

To run inelastic transport simulations as a script one can execute:

.. code-block:: bash

    >>> import Inelatica.iets as IETS
    >>> my_argv = '-F 8 -L 13 -p PHrun/Output.nc -f TSrun/RUN.fdf IETS_dir'
    >>> my_opts = IETS.GetOptions(my_argv)
    >>> IETS.main(my_opts)

Classes
-------

.. autosummary::
   :toctree:

   GetOptions
   main

.. currentmodule:: Inelastica.iets

"""
from __future__ import print_function

import numpy as N
import netCDF4 as NC4
import sys
import Inelastica.physics.constants as PC
import Inelastica.misc.valuecheck as VC
import Inelastica.io.log as Log
import Inelastica.NEGF as NEGF
import Inelastica.io.siesta as SIO
import Inelastica.MakeGeom as MG
import Inelastica.math as MM


def GetOptions(argv):
    """
    Returns an instance of ``options`` for the ``iets`` module

    Parameters
    ----------
    argv : string
        For example `-n 10 test_dir`, which instructs to compute 10 eigenchannels
        and place the results in the output directory `test_dir`.
    """

    # if text string is specified, convert to list
    if isinstance(argv, VC.string_types):
        argv = argv.split()

    import argparse
    p = argparse.ArgumentParser(description='Inelastica script calculates and writes LOE quantities in ascii (Systemlabel.IN) and NetCDF (Systemlabel.IN.nc)')
    p.add_argument('DestDir', help='Destination directory')
    p.add_argument('-n', '--NumChan', dest='numchan', type=int, default=4,
                   help='Number of eigenchannels [default: %(default)s]')
    p.add_argument('-F', '--DeviceFirst', dest='DeviceFirst', default=0, type=int,
                   help='First device atom (SIESTA numbering) [TS.TBT.PDOSFrom]')
    p.add_argument('-L', '--DeviceLast', dest='DeviceLast', default=0, type=int,
                   help='Last device atom (SIESTA numbering) [TS.TBT.PDOSTo]')
    p.add_argument('-e', '--Energy', dest='energy', default=0.0, type=float,
                   help='Energy reference where Greens functions etc are evaluated [default: %(default)s eV]')
    p.add_argument('--eta', dest='eta', default=0.000001, type=float,
                   help='Tiny imag. part in Greens functions etc. [default: %(default)s eV]')
    p.add_argument('-f', '--fdf', dest='fn', default='./RUN.fdf', type=str,
                   help='Input fdf-file for TranSIESTA calculations [default: %(default)s]')
    p.add_argument('-s', '--iSpin', dest='iSpin', default=0, type=int,
                   help='Spin channel [default: %(default)s]')
    p.add_argument('-x', '--k1', dest='k1', default=0.0, type=float,
                   help='k-point along a1 [default: %(default)s]')
    p.add_argument('-y', '--k2', dest='k2', default=0.0, type=float,
                   help='k-point along a2 [default: %(default)s]')
    p.add_argument('-p', '--PhononNetCDF', dest='PhononNetCDF', default='Output.nc', type=str,
                   help='Electron-phonon coupling NetCDF [default: %(default)s]')
    p.add_argument('-t', '--Temp', dest='Temp', default=4.2, type=float,
                   help='Temperature [default: %(default)s K]')
    p.add_argument('-b', '--BiasPoints', dest='biasPoints', default=801, type=int,
                   help='Number of bias points [default: %(default)s]')
    p.add_argument('-v', '--MaxBias', dest='maxBias', default=0.4, type=float,
                   help='Sets the IETS bias range (-MaxBias to MaxBias) [default: %(default)s V]')
    p.add_argument('-c', '--ModeCutoff', dest='modeCutoff', default='0.0025', type=float,
                   help='Ignore phonon modes with lower hw [default: %(default)s eV]')
    p.add_argument('-V', '--Vrms', dest='Vrms', default='0.005', type=float,
                   help='Lock in amplifier broadening [default: %(default)s V]')
    p.add_argument('-H', '--Heating', dest='PhHeating', default=False, action='store_true',
                   help='Include heating of vibrational modes [default: %(default)s]')
    p.add_argument('-d', '--PhExtDamp', dest='PhExtDamp', default=1e-15, type=float,
                   help='External damping [default: %(default)s (?) TODO check unit!]')
    p.add_argument('-u', '--useSigNC', dest='signc', default=False, action='store_true',
                   help='Use SigNCfiles [default: %(default)s]')
    p.add_argument('-l', '--etaLead', dest='etaLead', type=float, default=0.0,
                   help='Additional imaginary part added ONLY in the leads (surface GF) [default: %(default)s eV]')
    p.add_argument('--SpectralCutoff', dest='SpectralCutoff', type=float, default=1e-8,
                   help='Cutoff value for SpectralMatrix functions (for ordinary matrix representation set cutoff<=0.0) [default: %(default)s]')

    # Electrode stuff
    p.add_argument('--bulk', dest='UseBulk', default=-1, action='store_true',
                   help='Use bulk in electrodes. The Hamiltonian from the electrode calculation is inserted into the electrode region in the TranSIESTA cell [TS.UseBulkInElectrodes]')
    p.add_argument('--nobulk', dest='UseBulk', default=-1, action='store_false',
                   help='Use only self-energies in the electrodes. The full Hamiltonian of the TranSIESTA cell is used in combination with self-energies for the electrodes [TS.UseBulkInElectrodes]')

    # Scale (artificially) the coupling to the electrodes
    p.add_argument('--scaleSigL', dest='scaleSigL', type=float, default=1.0,
                   help='Scale factor applied to Sigma_L [default: %(default)s]')
    p.add_argument('--scaleSigR', dest='scaleSigR', type=float, default=1.0,
                   help='Scale factor applied to Sigma_R [default: %(default)s]')

    # Which LOE method?
    p.add_argument('--LOEscale', dest='LOEscale', type=float, default=1.0,
                   help='Scale factor to interpolate between LOE-WBA (0.0) and generalized LOE (1.0), see PRB 89, 081405(R) (2014) [default: %(default)s]')
    p.add_argument('--VfracL', dest='VfracL', type=float, default=0.5,
                   help='Voltage fraction over the left-center interface [default: %(default)s]')

    # Parse the options
    options = p.parse_args(argv)

    # Set module name
    options.module = 'Inelastica'

    # k-point
    options.kpoint = N.array([options.k1, options.k2, 0.0], N.float)
    del options.k1, options.k2

    return options


def main(options):
    """
    Main routine to compute inelastic transport characteristics (dI/dV, d2I/dV2, IETS etc)

    Parameters
    ----------
    options : an ``options`` instance
    """
    Log.CreatePipeOutput(options)
    VC.OptionsCheck(options)
    Log.PrintMainHeader(options)

    options.XV = '%s/%s.XV'%(options.head, options.systemlabel)
    options.geom = MG.Geom(options.XV, BufferAtoms=options.buffer)
    # Voltage fraction over left-center interface
    VfracL = options.VfracL # default is 0.5
    print('Inelastica: Voltage fraction over left-center interface: VfracL =', VfracL)
    # Set up electrodes and device Greens function
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL, options.NA1L, options.NA2L, options.voltage*VfracL)
    elecL.scaling = options.scaleSigL
    elecL.semiinf = options.semiinfL
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR, options.NA1R, options.NA2R, options.voltage*(VfracL-1.))
    elecR.scaling = options.scaleSigR
    elecR.semiinf = options.semiinfR
    # Read phonons
    NCfile = NC4.Dataset(options.PhononNetCDF, 'r')
    print('Inelastica: Reading ', options.PhononNetCDF)
    hw = NCfile.variables['hw'][:]
    # Work with GFs etc for positive (V>0: \mu_L>\mu_R) and negative (V<0: \mu_L<\mu_R) bias voltages
    GFp = NEGF.GF(options.TSHS, elecL, elecR,
                  Bulk=options.UseBulk, DeviceAtoms=options.DeviceAtoms,
                  BufferAtoms=options.buffer)
    # Prepare lists for various trace factors
    #GF.dGnout = []
    #GF.dGnin = []
    GFp.P1T = N.zeros(len(hw), N.float)     # M.A.M.A (total e-h damping)
    GFp.P2T = N.zeros(len(hw), N.float)     # M.AL.M.AR (emission)
    GFp.ehDampL = N.zeros(len(hw), N.float) # M.AL.M.AL (L e-h damping)
    GFp.ehDampR = N.zeros(len(hw), N.float) # M.AR.M.AR (R e-h damping)
    GFp.nHT = N.zeros(len(hw), N.float)     # non-Hilbert/Isym factor
    GFp.HT = N.zeros(len(hw), N.float)      # Hilbert/Iasym factor
    GFp.dIel = N.zeros(len(hw), N.float)
    GFp.dIinel = N.zeros(len(hw), N.float)
    GFp.dSel = N.zeros(len(hw), N.float)
    GFp.dSinel = N.zeros(len(hw), N.float)
    #
    GFm = NEGF.GF(options.TSHS, elecL, elecR,
                  Bulk=options.UseBulk, DeviceAtoms=options.DeviceAtoms,
                  BufferAtoms=options.buffer)
    GFm.P1T = N.zeros(len(hw), N.float)     # M.A.M.A (total e-h damping)
    GFm.P2T = N.zeros(len(hw), N.float)     # M.AL.M.AR (emission)
    GFm.ehDampL = N.zeros(len(hw), N.float) # M.AL.M.AL (L e-h damping)
    GFm.ehDampR = N.zeros(len(hw), N.float) # M.AR.M.AR (R e-h damping)
    GFm.nHT = N.zeros(len(hw), N.float)     # non-Hilbert/Isym factor
    GFm.HT = N.zeros(len(hw), N.float)      # Hilbert/Iasym factor
    GFm.dIel = N.zeros(len(hw), N.float)
    GFm.dIinel = N.zeros(len(hw), N.float)
    GFm.dSel = N.zeros(len(hw), N.float)
    GFm.dSinel = N.zeros(len(hw), N.float)
    # Calculate transmission at Fermi level
    GFp.calcGF(options.energy+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
               etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
    L = options.bufferL
    # Pad lasto with zeroes to enable basis generation...
    lasto = N.zeros((GFp.HS.nua+L+1,), N.int)
    lasto[L:] = GFp.HS.lasto
    basis = SIO.BuildBasis(options.fn, options.DeviceAtoms[0]+L, options.DeviceAtoms[1]+L, lasto)
    basis.ii -= L
    TeF = MM.trace(GFp.TT).real
    GFp.TeF = TeF
    GFm.TeF = TeF
    # Check consistency of PHrun vs TSrun inputs
    IntegrityCheck(options, GFp, NCfile)
    # Calculate trace factors one mode at a time
    print('Inelastica: LOEscale =', options.LOEscale)
    if options.LOEscale == 0.0:
        # LOEscale=0.0 => Original LOE-WBA method, PRB 72, 201101(R) (2005) [cond-mat/0505473].
        GFp.calcGF(options.energy+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                   etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
        GFm.calcGF(options.energy+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                   etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
        for ihw in (hw > options.modeCutoff).nonzero()[0]:
            calcTraces(options, GFp, GFm, basis, NCfile, ihw)
            calcTraces(options, GFm, GFp, basis, NCfile, ihw)
        writeFGRrates(options, GFp, hw, NCfile)
    else:
        # LOEscale=1.0 => Generalized LOE, PRB 89, 081405(R) (2014) [arXiv:1312.7625]
        for ihw in (hw > options.modeCutoff).nonzero()[0]:
            GFp.calcGF(options.energy+hw[ihw]*options.LOEscale*VfracL+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                       etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
            GFm.calcGF(options.energy+hw[ihw]*options.LOEscale*(VfracL-1.)+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                       etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
            calcTraces(options, GFp, GFm, basis, NCfile, ihw)
            if VfracL != 0.5:
                GFp.calcGF(options.energy-hw[ihw]*options.LOEscale*(VfracL-1.)+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                           etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
                GFm.calcGF(options.energy-hw[ihw]*options.LOEscale*VfracL+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                           etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
            calcTraces(options, GFm, GFp, basis, NCfile, ihw)

    # Multiply traces with voltage-dependent functions
    data = calcIETS(options, GFp, GFm, basis, hw)
    NCfile.close()
    NEGF.SavedSig.close()
    Log.PrintMainFooter(options)
    return data

########################################################


def IntegrityCheck(options, GF, NCfile):
    # Perform consistency checks for device region in
    # PH and TS calculations by comparing coordinates
    # and atom numbers
    PH_dev = NCfile.variables['DeviceAtoms'][:]
    PH_dyn = NCfile.variables['DynamicAtoms'][:]
    PH_xyz = NCfile.variables['GeometryXYZ'][:]
    PH_anr = NCfile.variables['AtomNumbers'][:]
    TS_dev = list(range(options.DeviceAtoms[0], options.DeviceAtoms[1]+1))
    TS_anr = options.geom.anr
    TS_xyz = options.geom.xyz
    print('\nInelastica.IntegrityCheck:')
    print('A = %s'%options.PhononNetCDF)
    print('B = %s'%options.XV)
    print(' idxA    xA       yA       zA    anrA   ', end=' ')
    print(' idxB    xB       yB       zB    anrB')
    for i in range(max(len(PH_dev), len(TS_dev))):
        # Geom A
        if PH_dev[0]+i in PH_dev:
            s = ('%i'%(PH_dev[0]+i)).rjust(5)
            for j in range(3):
                s += ('%.4f'%(PH_xyz[PH_dev[0]-1+i, j])).rjust(9)
            s += ('%i'%PH_anr[PH_dev[0]-1+i]).rjust(5)
        else:
            s = ('---').center(36)
        s += '  vs'
        # Geom B
        if options.DeviceAtoms[0]+i in TS_dev:
            s += ('%i'%(options.DeviceAtoms[0]+i)).rjust(5)
            for j in range(3):
                s += ('%.4f'%(TS_xyz[options.DeviceAtoms[0]-1+i, j])).rjust(9)
            s += ('%i'%TS_anr[options.DeviceAtoms[0]-1+i]).rjust(5)
        else:
            s += ('---').center(36)
        print(s)
        # Dynamic region
        if PH_dev[0]+i+1 == PH_dyn[0]:
            print('        -------------------- Dynamic region begins ---------------------')
        if PH_dev[0]+i == PH_dyn[-1]:
            print('        --------------------- Dynamic region ends ----------------------')
    # - check 1: Matrix sizes
    PH_H0 = NCfile.variables['H0'][:]
    check1 = N.shape(PH_H0[0]) == N.shape(GF.Gr)
    if check1:
        print('... Check 1 passed: Device orb. space matches')
    else:
        print('... Check 1 failed: Device orb. space do not match!!!')
    # - check 2&3: Geometry and atom number
    dist_xyz = 0.0
    dist_anr = 0.0
    for i in range(len(PH_dev)):
        # Geometric distance between atoms
        if i == 0:
            # Allow for a global offset of coordinates R
            R = PH_xyz[PH_dev[0]-1+i]-TS_xyz[options.DeviceAtoms[0]-1+i]
            print('Global offset R = [%.3f %.3f %.3f]'%(R[0], R[1], R[2]))
        d = PH_xyz[PH_dev[0]-1+i]-TS_xyz[options.DeviceAtoms[0]-1+i] - R
        dist_xyz += N.dot(d, d)**.5
        # Difference between atom numbers
        # Modulo 200
        a = (PH_anr[PH_dev[0]-1+i])%200 - (TS_anr[options.DeviceAtoms[0]-1+i])%200
        dist_anr += abs(a)
    if dist_xyz < 1e-3:
        print('... Check 2 passed: Atomic coordinates consistent')
        check2 = True
    elif dist_xyz < 0.1:
        print('... Check 2 WARNING: Atomic coordinates deviate by %.3f Ang!!!'%dist_xyz)
        check2 = True
    else:
        print('... Check 2 failed: Atomic coordinates deviate by %.3f Ang!!!'%dist_xyz)
        check2 = False
    if dist_anr == 0:
        print('... Check 3 passed: Atomic numbers consistent')
        check3 = True
    else:
        print('... Check 3 failed: Atomic numbers inconsistent!!!')
        check3 = False
    if (not check1) or (not check2) or (not check3):
        sys.exit('Inelastica: Error - inconsistency detected for device region.\n')


def get_coupling(ncfile, ihw, ispin):
    ihw = int(ihw)
    dims = ncfile.dimensions
    vars = ncfile.variables
    norb = len(dims["norb"])
    eph = N.arange(norb).reshape(-1, 1)
    if "norb_eph" in dims:
        eph = vars["eph_orbs_in_dev"][:].data.reshape(-1, 1)

    dtype = complex if "ImHe_ph" in vars else float
    M = N.zeros((norb, norb), dtype=dtype)
    M.real[eph, eph.T] = vars['He_ph'][ihw, ispin, :, :]
    if "ImHe_ph" not in vars:
        print('Warning: Variable ImHe_ph not found')
    else:
        M.imag[eph, eph.T] = vars['ImHe_ph'][ihw, ispin, :, :]
    return M


def calcTraces(options, GF1, GF2, basis, NCfile, ihw):
    # Calculate various traces over the electronic structure
    # Electron-phonon couplings
    ihw = int(ihw)
    M = get_coupling(NCfile, ihw, options.iSpin)
    # Calculation of intermediate quantity
    MARGLGM = MM.mm(M, GF1.ARGLG, M)
    MARGLGM2 = MM.mm(M, GF2.ARGLG, M)
    # LOE expressions in compact form
    t1 = MM.mm(MARGLGM, GF2.AR)
    t2 = MM.mm(MARGLGM2, GF1.AL)
    # Note that compared with Eq. (10) of PRB89, 081405 (2014) we here use
    # the definition B_lambda = MM.trace(t1-dagger(t2)), which in turn gives
    # ReB = MM.trace(t1).real-MM.trace(t2).real
    # ImB = MM.trace(t1).imag+MM.trace(t2).imag
    K23 = MM.trace(t1).imag+MM.trace(t2).imag
    K4 = MM.trace(MM.mm(M, GF1.ALT, M, GF2.AR))
    aK23 = 2*(MM.trace(t1).real-MM.trace(t2).real) # asymmetric part
    # Non-Hilbert term defined here with a minus sign
    GF1.nHT[ihw] = NEGF.AssertReal(K23+K4, 'nHT[%i]'%ihw)
    GF1.HT[ihw] = NEGF.AssertReal(aK23, 'HT[%i]'%ihw)
    # Power, damping and current rates
    GF1.P1T[ihw] = NEGF.AssertReal(MM.trace(MM.mm(M, GF1.A, M, GF2.A)), 'P1T[%i]'%ihw)
    GF1.P2T[ihw] = NEGF.AssertReal(MM.trace(MM.mm(M, GF1.AL, M, GF2.AR)), 'P2T[%i]'%ihw)
    GF1.ehDampL[ihw] = NEGF.AssertReal(MM.trace(MM.mm(M, GF1.AL, M, GF2.AL)), 'ehDampL[%i]'%ihw)
    GF1.ehDampR[ihw] = NEGF.AssertReal(MM.trace(MM.mm(M, GF1.AR, M, GF2.AR)), 'ehDampR[%i]'%ihw)
    # Remains from older version (see before rev. 219):
    #GF.dGnout.append(EC.calcCurrent(options,basis,GF.HNO,mm(Us,-0.5j*(tmp1-dagger(tmp1)),Us)))
    #GF.dGnin.append(EC.calcCurrent(options,basis,GF.HNO,mm(Us,mm(G,MA1M,Gd)-0.5j*(tmp2-dagger(tmp2)),Us)))
    # NB: TF Should one use GF.HNO (nonorthogonal) or GF.H (orthogonalized) above?

    if options.LOEscale == 0.0:
        # Check against original LOE-WBA formulation
        isym1 = MM.mm(GF1.ALT, M, GF2.AR, M)
        isym2 = MM.mm(MM.dagger(GF1.ARGLG), M, GF2.A, M)
        isym3 = MM.mm(GF1.ARGLG, M, GF2.A, M)
        isym = MM.trace(isym1)+1j/2.*(MM.trace(isym2)-MM.trace(isym3))
        print('LOE-WBA check: Isym diff', K23+K4-isym)
        iasym1 = MM.mm(MM.dagger(GF1.ARGLG), M, GF2.AR-GF2.AL, M)
        iasym2 = MM.mm(GF1.ARGLG, M, GF2.AR-GF2.AL, M)
        iasym = MM.trace(iasym1)+MM.trace(iasym2)
        print('LOE-WBA check: Iasym diff', aK23-iasym)

        # Compute inelastic shot noise terms according to the papers
        # Haupt, Novotny & Belzig, PRB 82, 165441 (2010) and
        # Avriller & Frederiksen, PRB 86, 155411 (2012)
        # Zero-temperature limit
        TT = MM.mm(GF1.GammaL, GF1.AR) # this matrix has the correct shape for MM
        ReGr = (GF1.Gr+GF1.Ga)/2.
        tmp = MM.mm(GF1.Gr, M, ReGr, M, GF1.AR)
        tmp = tmp+MM.dagger(tmp)
        Tlambda0 = MM.mm(GF1.GammaL, tmp)
        tmp1 = MM.mm(M, GF1.AR, M)
        tmp2 = MM.mm(M, GF1.A, M, GF1.Gr, GF1.GammaR)
        tmp = tmp1+1j/2.*(MM.dagger(tmp2)-tmp2)
        Tlambda1 = MM.mm(GF1.GammaL, GF1.Gr, tmp, GF1.Ga)
        MARGL = MM.mm(M, GF1.AR, GF1.GammaL)
        tmp1 = MM.mm(MARGL, GF1.AR, M)
        tmp2 = MM.mm(MARGL, GF1.Gr, M, GF1.Gr, GF1.GammaR)
        tmp = tmp1+tmp2
        tmp = tmp + MM.dagger(tmp)
        Qlambda = MM.mm(-GF1.Ga, GF1.GammaL, GF1.Gr, tmp)
        tmp = -2*TT
        OneMinusTwoT = tmp+N.identity(len(GF1.GammaL))
        # Store relevant traces
        GF1.dIel[ihw] = NEGF.AssertReal(MM.trace(Tlambda0), 'dIel[%i]'%ihw)
        GF1.dIinel[ihw] = NEGF.AssertReal(MM.trace(Tlambda1), 'dIinel[%i]'%ihw)
        GF1.dSel[ihw] = NEGF.AssertReal(MM.trace(MM.mm(OneMinusTwoT, Tlambda0)), 'dSel[%i]'%ihw)
        GF1.dSinel[ihw] = NEGF.AssertReal(MM.trace(Qlambda+MM.mm(OneMinusTwoT, Tlambda1)), 'dSinel[%i]'%ihw)


def calcIETS(options, GFp, GFm, basis, hw):
    # Calculate product of electronic traces and voltage functions
    print('Inelastica.calcIETS: nHTp =', GFp.nHT*PC.unitConv) # OK
    print('Inelastica.calcIETS: nHTm =', GFm.nHT*PC.unitConv) # OK
    print('Inelastica.calcIETS: HTp  =', GFp.HT) # OK
    print('Inelastica.calcIETS: HTm  =', GFm.HT) # OK

    # Set up grid and Hilbert term
    kT = options.Temp/11604.0 # (eV)

    # Generate grid for numerical integration of Hilbert term
    max_hw = max(hw)
    max_win = max(-options.minBias, max_hw)+20*kT+4*options.Vrms
    min_win = min(-options.maxBias, -max_hw)-20*kT-4*options.Vrms
    pts = int(N.floor((max_win-min_win)/kT*3))
    Egrid = N.linspace(min_win, max_win, pts)
    print("Inelastica.calcIETS: Hilbert integration grid : %i pts [%f,%f]" % (pts, min(Egrid), max(Egrid)))

    NN = options.biasPoints
    print('Inelastica.calcIETS: Biaspoints =', NN)

    # Add some points for the Lock in broadening
    approxdV = (options.maxBias-options.minBias)/(NN-1)
    NN += int(((8*options.Vrms)/approxdV)+.5)
    Vl = N.linspace(options.minBias-4*options.Vrms, options.maxBias+4*options.Vrms, NN)

    # Vector implementation on Vgrid:
    wp = (1+N.sign(Vl))/2. # weights for positive V
    wm = (1-N.sign(Vl))/2. # weights for negative V

    # Mode occupation and power dissipation
    Pow = N.zeros((len(hw), NN), N.float) # (modes,Vgrid)
    nPh = N.zeros((len(hw), NN), N.float)
    t0 = N.clip(Vl/kT, -700, 700)
    cosh0 = N.cosh(t0) # Vgrid
    sinh0 = N.sinh(t0)
    for i in (hw > options.modeCutoff).nonzero()[0]:
        P1T = wm*GFm.P1T[i]+wp*GFp.P1T[i]
        P2T = wm*GFm.P2T[i]+wp*GFp.P2T[i]
        # Bose distribution
        nB = 1/(N.exp(N.clip(hw[i]/kT, -300, 300))-1) # number
        t1 = N.clip(hw[i]/(2*kT), -700, 700) # number
        coth1 = N.cosh(t1)/N.sinh(t1)
        # Emission rate and e-h damping
        damp = P1T*hw[i]/N.pi # Vgrid
        emis = P2T*(hw[i]*(cosh0-1)*coth1-Vl*sinh0)/(N.cosh(2*t1)-cosh0)/N.pi
        # Determine mode occupation
        if options.PhHeating:
            nPh[i, :] = emis/(hw[i]*P1T/N.pi+options.PhExtDamp)+nB
        else:
            nPh[i, :] = nB
        # Mode-resolved power dissipation
        Pow[i, :] = hw[i]*((nB-nPh[i])*damp+emis)

    # Current: non-Hilbert part (InH)
    InH = N.zeros((NN,), N.float) # Vgrid
    IsymF = N.zeros((NN,), N.float)
    for i in (hw > options.modeCutoff).nonzero()[0]:
        nHT = wm*GFm.nHT[i]+wp*GFp.nHT[i] # Vgrid
        t1 = hw[i]/(2*kT) # number
        t1 = N.clip(t1, -700, 700)
        coth1 = N.cosh(t1)/N.sinh(t1)
        t2 = (hw[i]+Vl)/(2*kT) # Vgrid
        t2 = N.clip(t2, -700, 700)
        coth2 = N.cosh(t2)/N.sinh(t2)
        t3 = (hw[i]-Vl)/(2*kT) # Vgrid
        t3 = N.clip(t3, -700, 700)
        coth3 = N.cosh(t3)/N.sinh(t3) # Vgrid
        # Isym function
        Isym = 0.5*(hw[i]+Vl)*(coth1-coth2) # Vgrid
        Isym -= 0.5*(hw[i]-Vl)*(coth1-coth3)
        # non-Hilbert part
        InH += (Isym+2*Vl*nPh[i])*nHT # Vgrid
        IsymF += Isym

    # Current: Add Landauer part, GFm.TeF = GFp.TeF
    InH += GFp.TeF*Vl # Vgrid

    # Current: Asymmetric/Hilbert part (IH)
    try:
        import scipy.special as SS
        print("Inelastica: Computing asymmetric term using digamma function,")
        print("... see G. Bevilacqua et al., Eur. Phys. J. B (2016) 89: 3")
        IH = N.zeros((NN,), N.float)
        IasymF = N.zeros((NN,), N.float)
        for i in (hw > options.modeCutoff).nonzero()[0]:
            v0 = hw[i]/(2*N.pi*kT)
            vp = (hw[i]+Vl)/(2*N.pi*kT)
            vm = (hw[i]-Vl)/(2*N.pi*kT)
            Iasym = kT*(2*v0*SS.psi(1.j*v0)-vp*SS.psi(1.j*vp)-vm*SS.psi(1.j*vm)).real
            IasymF += Iasym
            IH += GFp.HT[i]*N.array(Vl > 0.0, dtype=int)*Iasym
            IH += GFm.HT[i]*N.array(Vl < 0.0, dtype=int)*Iasym
    except:
        print("Computing using explit Hilbert transformation")
        IH = N.zeros((NN,), N.float)
        IasymF = N.zeros((NN,), N.float)
        # Prepare box/window function on array
        Vl2 = N.outer(Vl, N.ones(Egrid.shape))
        Egrid2 = N.outer(N.ones(Vl.shape), Egrid)
        # Box/window function nF(E-Vl2)-nF(E-0):
        kasse = MM.box(0, -Vl2, Egrid2, kT) # (Vgrid,Egrid)
        ker = None
        for i in (hw > options.modeCutoff).nonzero()[0]:
            # Box/window function nF(E-hw)-nF(E+hw)
            tmp = MM.box(-hw[i], hw[i], Egrid, kT)
            hilb, ker = MM.Hilbert(tmp, ker) # Egrid
            # Calculate Iasym for each bias point
            for j in range(len(Vl)):
                Iasym = MM.trapez(Egrid, kasse[j]*hilb, equidistant=True).real/2
                IasymF[j] += Iasym
                if Vl[j] > 0:
                    IH[j] += GFp.HT[i]*Iasym
                else:
                    IH[j] += GFm.HT[i]*Iasym

    # Compute inelastic shot noise terms here:
    absVl = N.absolute(Vl)
    Inew = N.zeros(len(Vl), N.float)
    Snew = N.zeros(len(Vl), N.float)
    print('Noise factors:')
    print(GFp.dIel)
    print(GFp.dIinel)
    print(GFp.dSel)
    print(GFp.dSinel)
    for i in (hw > options.modeCutoff).nonzero()[0]:
        # Elastic part
        Inew += GFp.dIel[i]*Vl
        Snew += GFp.dSel[i]*absVl
        # Inelastic part
        indx = (absVl-hw[i] < 0).nonzero()[0]
        fct = absVl-hw[i]
        fct[indx] = 0.0 # set elements to zero
        Inew += GFp.dIinel[i]*fct*N.sign(Vl)
        Snew += GFp.dSinel[i]*fct

    # Get the right units for gamma_eh, gamma_heat
    gamma_eh_p = N.zeros((len(hw),), N.float)
    gamma_eh_m = N.zeros((len(hw),), N.float)
    gamma_heat_p = N.zeros((len(hw),), N.float)
    gamma_heat_m = N.zeros((len(hw),), N.float)
    for i in (hw > options.modeCutoff).nonzero()[0]:
        # Units [Phonons per Second per dN where dN is number extra phonons]
        gamma_eh_p[i] = GFp.P1T[i]*hw[i]*PC.unitConv
        gamma_eh_m[i] = GFm.P1T[i]*hw[i]*PC.unitConv
        # Units [Phonons per second per eV [eV-ihw]
        gamma_heat_p[i] = GFp.P2T[i]*PC.unitConv
        gamma_heat_m[i] = GFm.P2T[i]*PC.unitConv

    print('Inelastica.calcIETS: gamma_eh_p =', gamma_eh_p) # OK
    print('Inelastica.calcIETS: gamma_eh_m =', gamma_eh_m) # OK
    print('Inelastica.calcIETS: gamma_heat_p =', gamma_heat_p) # OK
    print('Inelastica.calcIETS: gamma_heat_m =', gamma_heat_m) # OK

    V, I, dI, ddI, BdI, BddI = Broaden(options, Vl, InH+IH)
    V, Is, dIs, ddIs, BdIs, BddIs = Broaden(options, Vl, IsymF)
    V, Ia, dIa, ddIa, BdIa, BddIa = Broaden(options, Vl, IasymF)

    # Interpolate quantities to new V-grid
    NPow = N.zeros((len(Pow), len(V)), N.float)
    NnPh = N.zeros((len(Pow), len(V)), N.float)
    for ii in range(len(Pow)):
        NPow[ii] = MM.interpolate(V, Vl, Pow[ii])
        NnPh[ii] = MM.interpolate(V, Vl, nPh[ii])

    # Interpolate inelastic noise
    NV, NI, NdI, NddI, NBdI, NBddI = Broaden(options, Vl, GFp.TeF*Vl+Inew)
    NV, NS, NdS, NddS, NBdS, NBddS = Broaden(options, Vl, Snew)

    print('Inelastica.calcIETS: V[:5]        =', V[:5]) # OK
    print('Inelastica.calcIETS: V[-5:][::-1] =', V[-5:][::-1]) # OK
    print('Inelastica.calcIETS: I[:5]        =', I[:5]) # OK
    print('Inelastica.calcIETS: I[-5:][::-1] =', I[-5:][::-1]) # OK
    print('Inelastica.calcIETS: BdI[:5]        =', BdI[:5]) # OK
    print('Inelastica.calcIETS: BdI[-5:][::-1] =', BdI[-5:][::-1]) # OK
    print('Inelastica.calcIETS: BddI[:5]        =', BddI[:5]) # OK
    print('Inelastica.calcIETS: BddI[-5:][::-1] =', BddI[-5:][::-1]) # OK

    datafile = '%s/%s.IN'%(options.DestDir, options.systemlabel)
    # ascii format
    writeLOEData2Datafile(datafile+'p', hw, GFp.TeF, GFp.nHT, GFp.HT)
    writeLOEData2Datafile(datafile+'m', hw, GFm.TeF, GFm.nHT, GFm.HT)
    # netcdf format
    outNC = initNCfile(datafile, hw, V)
    write2NCfile(outNC, BddI/BdI, 'IETS', 'Broadened BddI/BdI [1/V]')
    write2NCfile(outNC, ddI/dI, 'IETS_0', 'Intrinsic ddI/dI [1/V]')
    write2NCfile(outNC, BdI, 'BdI', 'Broadened BdI, G0')
    write2NCfile(outNC, BddI, 'BddI', 'Broadened BddI, G0')
    write2NCfile(outNC, I, 'I', 'Intrinsic I, G0 V')
    write2NCfile(outNC, dI, 'dI', 'Intrinsic dI, G0')
    write2NCfile(outNC, ddI, 'ddI', 'Intrinsic ddI, G0/V')
    if options.LOEscale == 0.0:
        write2NCfile(outNC, GFp.nHT, 'ISymTr', 'Trace giving Symmetric current contribution (prefactor to universal function)')
        write2NCfile(outNC, GFp.HT, 'IAsymTr', 'Trace giving Asymmetric current contribution (prefactor to universal function)')
        write2NCfile(outNC, gamma_eh_p, 'gamma_eh', 'e-h damping [*deltaN=1/Second]')
        write2NCfile(outNC, gamma_heat_p, 'gamma_heat', 'Phonon heating [*(bias-hw) (eV) = 1/Second]')
        # New stuff related to the noise implementation
        write2NCfile(outNC, NI, 'Inew', 'Intrinsic Inew (new implementation incl. elastic renormalization, T=0)')
        write2NCfile(outNC, NdI, 'dInew', 'Intrinsic dInew (new implementation incl. elastic renormalization, T=0)')
        write2NCfile(outNC, NddI, 'ddInew', 'Intrinsic ddInew (new implementation incl. elastic renormalization, T=0)')
        write2NCfile(outNC, NddI/NdI, 'IETSnew_0', 'Intrinsic ddInew/dInew (new implementation incl. elastic renormalization, T=0) [1/V]')
        write2NCfile(outNC, NBdI, 'BdInew', 'Broadened BdInew (new implementation incl. elastic renormalization, T=0)')
        write2NCfile(outNC, NBddI, 'BddInew', 'Broadened BddInew (new implementation incl. elastic renormalization, T=0)')
        write2NCfile(outNC, NBddI/NBdI, 'IETSnew', 'Broadened BddInew/BdInew (new implementation incl. elastic renormalization, T=0) [1/V]')
        write2NCfile(outNC, NdS, 'dSnew', 'Inelastic first-derivative of the shot noise dSnew (T=0)')
    else:
        write2NCfile(outNC, GFp.nHT, 'ISymTr_p', 'Trace giving Symmetric current contribution (prefactor to universal function)')
        write2NCfile(outNC, GFp.HT, 'IAsymTr_p', 'Trace giving Asymmetric current contribution (prefactor to universal function)')
        write2NCfile(outNC, GFm.nHT, 'ISymTr_m', 'Trace giving Symmetric current contribution (prefactor to universal function)')
        write2NCfile(outNC, GFm.HT, 'IAsymTr_m', 'Trace giving Asymmetric current contribution (prefactor to universal function)')
        write2NCfile(outNC, gamma_eh_p, 'gamma_eh_p', 'e-h damping [*deltaN=1/Second]')
        write2NCfile(outNC, gamma_heat_p, 'gamma_heat_p', 'Phonon heating [*(bias-hw) (eV) = 1/Second]')
        write2NCfile(outNC, gamma_eh_m, 'gamma_eh_m', 'e-h damping [*deltaN=1/Second]')
        write2NCfile(outNC, gamma_heat_m, 'gamma_heat_m', 'Phonon heating [*(bias-hw) (eV) = 1/Second]')
    # Phonon occupations and power balance
    write2NCfile(outNC, NnPh, 'nPh', 'Number of phonons')
    write2NCfile(outNC, N.sum(NnPh, axis=0), 'nPh_tot', 'Total number of phonons')
    write2NCfile(outNC, NPow, 'Pow', 'Mode-resolved power balance')
    write2NCfile(outNC, N.sum(NPow, axis=0), 'Pow_tot', 'Total power balance')
    # Write "universal functions"
    write2NCfile(outNC, dIs, 'dIs', 'dIsym function')
    write2NCfile(outNC, dIa, 'dIa', 'dIasym function')
    write2NCfile(outNC, ddIs, 'ddIs', 'ddIasym function')
    write2NCfile(outNC, ddIa, 'ddIa', 'ddIasym function')
    # Write energy reference where Greens functions are evaluated
    outNC.createDimension('number', 1)
    tmp = outNC.createVariable('EnergyRef', 'd', ('number',))
    tmp[:] = N.array(options.energy)
    # Write LOEscale
    tmp = outNC.createVariable('LOEscale', 'd', ('number',))
    tmp[:] = N.array(options.LOEscale)
    # Write k-point
    outNC.createDimension('vector', 3)
    tmp = outNC.createVariable('kpoint', 'd', ('vector',))
    tmp[:] = N.array(options.kpoint)
    outNC.close()

    return V, I, dI, ddI, BdI, BddI


########################################################
def writeFGRrates(options, GF, hw, NCfile):
    print('Inelastica.writeFGRrates: Computing FGR rates')
    # Eigenchannels
    GF.calcEigChan(channels=options.numchan)
    NCfile = NC4.Dataset(options.PhononNetCDF, 'r')
    print('Reading ', options.PhononNetCDF)

    outFile = open('%s/%s.IN.FGR'%(options.DestDir, options.systemlabel), 'w')
    outFile.write('Total transmission [in units of (1/s/eV)] : %e\n' % (PC.unitConv*GF.TeF,))

    for ihw in range(len(hw)):
        SIO.printDone(ihw, len(hw), 'Golden Rate')
        M = get_coupling(NCfile, ihw, options.iSpin)
        rate = N.zeros((len(GF.ECleft), len(GF.ECright)), N.float)
        totrate = 0.0
        inter, intra = 0.0, 0.0 # splitting total rate in two
        for iL in range(len(GF.ECleft)):
            for iR in range(len(GF.ECright)):
                tmp = N.dot(N.conjugate(GF.ECleft[iL]), MM.mm(M, GF.ECright[iR]))
                rate[iL, iR] = (2*N.pi)**2*abs(tmp)**2
                totrate += rate[iL, iR]
                if iL == iR: intra += rate[iL, iR]
                else: inter += rate[iL, iR]

        outFile.write('\nPhonon mode %i : %f eV [Rates in units of (1/s/eV)]\n' % (ihw, hw[ihw]))
        outFile.write('eh-damp : %e (1/s) , heating %e (1/(sV)))\n' % (GF.P1T[ihw]*PC.unitConv*hw[ihw], GF.P2T[ihw]*PC.unitConv))
        outFile.write('eh-damp 1, 2 (MALMAL, MARMAR): %e (1/s) , %e (1/(s)))\n' % (GF.ehDampL[ihw]*PC.unitConv*hw[ihw], GF.ehDampR[ihw]*PC.unitConv*hw[ihw]))
        outFile.write('SymI : %e (1/(sV)) , AsymI %e (?))\n' % (GF.nHT[ihw]*PC.unitConv, GF.HT[ihw]*PC.unitConv))
        #outFile.write('Elast : %e (1/(sV)) , Inelast %e (1/(sV)))\n' % (GF.nHTel[ihw]*PC.unitConv,GF.nHTin[ihw]*PC.unitConv))
        outFile.write('down=left EC, right=right EC\n')
        if GF.P2T[ihw] > 0.0:
            if abs(totrate/(GF.P2T[ihw])-1) < 0.05:
                outFile.write('Sum/Tr[MALMAR] , Tr: %1.3f  %e\n'%(totrate/(GF.P2T[ihw]), PC.unitConv*GF.P2T[ihw]))
            else:
                outFile.write('WARNING: !!!! Sum/Tr[MALMAR] , Tr: %2.2e  %e\n'%(totrate/(GF.P2T[ihw]), PC.unitConv*GF.P2T[ihw]))
        else:
            outFile.write(' Tr:  %e\n'%(PC.unitConv*GF.P2T[ihw]))
        inter = inter/GF.P2T[ihw]
        intra = intra/GF.P2T[ihw]
        outFile.write('Interchannel ratio: Sum(inter)/Tr[MALMAR]      = %.4f \n'%inter)
        outFile.write('Intrachannel ratio: Sum(intra)/Tr[MALMAR]      = %.4f \n'%intra)
        outFile.write('Inter+intra ratio: Sum(inter+intra)/Tr[MALMAR] = %.4f \n'%(inter+intra))
        for iL in range(len(GF.ECleft)):
            for iR in range(len(GF.ECright)):
                outFile.write('%e ' % (PC.unitConv*rate[iL, iR],))
            outFile.write('\n')
    outFile.close()


################################################################
# Broadening due to Vrms
################################################################

def Broaden(options, VV, II):
    """
    Broadening corresponding to Lock in measurements for the
    conductance and IETS spectra. Also resample II, Pow, and nPh
    to match a common voltage list
    """

    II = II.copy()
    II = II.real

    # First derivative dI and bias list dV
    dI = (II[1:len(II)]-II[:-1])/(VV[1]-VV[0])
    dV = (VV[1:len(VV)]+VV[:-1])/2
    # Second derivative and bias ddV
    ddI = (dI[1:len(dI)]-dI[:-1])/(VV[1]-VV[0])
    ddV = (dV[1:len(dV)]+dV[:-1])/2

    # Modulation amplitude
    VA = N.sqrt(2.0)*options.Vrms

    # New bias ranges for broadening
    tmp = int(N.floor(VA/(dV[1]-dV[0]))+1)
    BdV = dV[tmp:-tmp]
    BddV = ddV[tmp:-tmp]

    # Initiate derivatives
    BdI = 0*BdV
    BddI = 0*BddV

    # Calculate first derivative with Vrms broadening
    for iV, V in enumerate(BdV):
        SIO.printDone(iV, len(BdV), 'Inelastica.Broaden: First-derivative Vrms broadening')
        wt = (N.array(list(range(200)))/200.0-0.5)*N.pi
        VL = V+VA*N.sin(wt)
        dIL = MM.interpolate(VL, dV, dI)
        BdI[iV] = 2/N.pi*N.sum(dIL*(N.cos(wt)**2))*(wt[1]-wt[0])

    # Calculate second derivative with Vrms broadening
    for iV, V in enumerate(BddV):
        SIO.printDone(iV, len(BddV), 'Inelastica.Broaden: Second-derivative Vrms broadening')
        wt = (N.array(list(range(200)))/200.0-0.5)*N.pi
        VL = V+VA*N.sin(wt)
        ddIL = MM.interpolate(VL, ddV, ddI)
        BddI[iV] = 8.0/3.0/N.pi*N.sum(ddIL*(N.cos(wt)**4))*(wt[1]-wt[0])

    # Reduce to one voltage grid
    NN = options.biasPoints
    V = N.linspace(options.minBias, options.maxBias, NN)

    NI = MM.interpolate(V, VV, II)
    NdI = MM.interpolate(V, dV, dI)
    NddI = MM.interpolate(V, ddV, ddI)
    NBdI = MM.interpolate(V, BdV, BdI)
    NBddI = MM.interpolate(V, BddV, BddI)

    return V, NI, NdI, NddI, NBdI, NBddI

################################################################
# Output to NetCDF file
################################################################


def initNCfile(filename, hw, V):
    'Initiate netCDF file'
    print('Inelastica: Initializing nc-file')
    ncfile = NC4.Dataset(filename+'.nc', 'w')
    ncfile.title = 'Inelastica Output'
    ncfile.version = 1
    ncfile.createDimension('Nph', len(hw))
    # Phonon mode frequencies
    tmp = ncfile.createVariable('hw', 'd', ('Nph',))
    tmp[:] = N.array(hw)
    tmp.units = 'eV'
    # Voltage
    ncfile.createDimension('V', len(V))
    tmp = ncfile.createVariable('V', 'd', ('V',))
    tmp[:] = N.array(V)
    tmp.units = 'V'
    return ncfile


def write2NCfile(NCfile, var, name, unit):
    var = N.array(var)
    dim1 = len(NCfile.variables['hw'][:])
    dim2 = len(NCfile.variables['V'][:])
    print('Inelastica.write2NCfile: Writing', name, var.shape)
    if var.shape == (dim1,):
        tmp = NCfile.createVariable(name, 'd', ('Nph',))
    elif var.shape == (dim2,):
        tmp = NCfile.createVariable(name, 'd', ('V',))
    elif var.shape == (dim1, dim2):
        tmp = NCfile.createVariable(name, 'd', ('Nph', 'V'))
    tmp[:] = var.real
    tmp.units = unit


def writeLOEData2Datafile(file, hw, T, nHT, HT):
    f = open(file, 'w')
    f.write("## Almost First Order Born Approximation (kbT=0)\n")
    f.write("## hw(eV)      Trans        non-H-term (e/pi/hbar)   H-term\n")
    for ii, val in enumerate(hw):
        f.write("## %e %e %e %e\n" % (val, T, nHT[ii], HT[ii]))
    f.write("\n")
    f.close()
