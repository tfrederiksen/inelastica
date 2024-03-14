"""

:mod:`Inelastica.EigenChannels`
===============================

1. Eigenchannels, method from Paulsson and Brandbyge PRB 2007
2. Calculate "bond" currents

Scripting
---------

To run ``EigenChannels`` as a script one can execute:

.. code-block:: bash

    >>> import Inelatica.EigenChannels as EC
    >>> my_argv = '-F 8 -L 13 -f TSrun/RUN.fdf EC_dir'
    >>> my_opts = EC.GetOptions(my_argv)
    >>> EC.main(my_opts)

Classes
-------

.. autosummary::
   :toctree:

   GetOptions
   main

.. currentmodule:: Inelastica.EigenChannels

"""
from __future__ import print_function

import numpy as N
import netCDF4 as NC4
import struct
import Inelastica.physics.constants as PC
import Inelastica.misc.valuecheck as VC
import Inelastica.io.log as Log
import Inelastica.NEGF as NEGF
import Inelastica.io.siesta as SIO
import Inelastica.MakeGeom as MG
import Inelastica.math as MM


def GetOptions(argv):
    """
    Returns an instance of ``options`` for the ``EigenChannels`` module

    Parameters
    ----------
    argv : string
        For example `-n 2 test_dir`, which instructs to compute only the two most transmitting
        eigenchannel scattering states and place the results in the output directory `test_dir`.
    """

    # if text string is specified, convert to list
    if isinstance(argv, VC.string_types):
        argv = argv.split()

    import argparse

    p = argparse.ArgumentParser(description='Eigenchannels, see Paulsson et al. PRB 76, 115117 (2007)')
    p.add_argument('DestDir', help='Destination directory')
    p.add_argument('-F', '--DeviceFirst', dest='DeviceFirst', default=0, type=int,
                   help='First device atom (SIESTA numbering) [default: TS.TBT.PDOSFrom]')
    p.add_argument('-L', '--DeviceLast', dest='DeviceLast', default=0, type=int,
                   help='Last device atom (SIESTA numbering) [default: TS.TBT.PDOSTo]')
    p.add_argument('-n', '--NumChan', dest='numchan', type=int, default=4,
                   help='Number of eigenchannels [default: %(default)s]')
    p.add_argument('-B', '--BothSides', dest='bothsides', default=False, action='store_true',
                   help='Calculate eigenchannels from both sides [default: %(default)s]')
    p.add_argument('-M', '--MPSH', dest='MolStates', default=0.0, type=float,
                   help='Calculate eigenstates of the device region Hamiltonian (Molecular Projected Selfconsistent Hamiltonian, MPSH) within [default: +/- %(default)s] eV from Ef')
    p.add_argument('-r', '--Res', dest='res', default=0.4, type=float,
                   help='Resolution [default: %(default)s Ang]')
    p.add_argument('-w', '--format', dest='format', default='XSF', type=str,
                   help='Wavefunction format (macu, cube, XSF, or nc) [default: %(default)s]')
    p.add_argument('-e', '--Energy', dest='energy', default=0.0, type=float,
                   help='Energy where eigenchannel scattering states are evaluated [default: %(default)s eV]')
    p.add_argument('--eta', dest='eta', type=float, default=0.000001,
                   help='Imaginary part added to all energies (device and leads) [default: %(default)s eV]')
    p.add_argument('-l', '--etaLead', dest='etaLead', type=float, default=0.0,
                   help='Additional imaginary part added ONLY in the leads (surface GF) [default: %(default)s eV]')
    p.add_argument('-f', '--fdf', dest='fn', default='./RUN.fdf', type=str,
                   help='Input fdf-file for TranSIESTA calculations [default: %(default)s]')
    p.add_argument('-s', '--iSpin', dest='iSpin', default=0, type=int,
                   help='Spin channel [default: %(default)s]')
    p.add_argument('-x', '--k1', dest='k1', default=0.0, type=float,
                   help='k-point along a1 [default: %(default)s]')
    p.add_argument('-y', '--k2', dest='k2', default=0.0, type=float,
                   help='k-point along a2 [default: %(default)s]')
    p.add_argument('-u', '--useSigNC', dest='signc', default=False, action='store_true',
                   help='Use SigNCfiles [default: %(default)s]')

    # Electrode stuff
    p.add_argument('--bulk', dest='UseBulk', default=-1, action='store_true',
                   help='Use bulk in electrodes. The Hamiltonian from the electrode calculation is inserted into the electrode region in the TranSIESTA cell [default: TS.UseBulkInElectrodes]')
    p.add_argument('--nobulk', dest='UseBulk', default=-1, action='store_false',
                 help='Use only self-energies in the electrodes. The full Hamiltonian of the TranSIESTA cell is used in combination with self-energies for the electrodes [default: TS.UseBulkInElectrodes]')

    # Scale (artificially) the coupling to the electrodes
    p.add_argument('--scaleSigL', dest='scaleSigL', type=float, default=1.0,
                   help='Scale factor applied to Sigma_L [default=%(default)s]')
    p.add_argument('--scaleSigR', dest='scaleSigR', type=float, default=1.0,
                   help='Scale factor applied to Sigma_R [default=%(default)s]')

    # Use spectral matrices?
    p.add_argument('--SpectralCutoff', dest='SpectralCutoff', type=float, default=0.0,
                   help='Cutoff value for SpectralMatrix functions (for ordinary matrix representation set cutoff<=0.0) [default=%(default)s]')

    options = p.parse_args(argv)

    # Set module name
    options.module = 'EigenChannels'

    # k-point
    options.kpoint = N.array([options.k1, options.k2, 0.0], N.float)
    del options.k1, options.k2

    return options


def main(options):
    """
    Main routine to compute eigenchannel scattering states

    Parameters
    ----------
    options : an ``options`` instance
    """

    Log.CreatePipeOutput(options)
    VC.OptionsCheck(options)
    Log.PrintMainHeader(options)

    # Read geometry
    XV = '%s/%s.XV'%(options.head, options.systemlabel)
    geom = MG.Geom(XV, BufferAtoms=options.buffer)

    # Set up device Greens function
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL, options.NA1L, options.NA2L, options.voltage/2.)
    elecL.scaling = options.scaleSigL
    elecL.semiinf = options.semiinfL
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR, options.NA1R, options.NA2R, -options.voltage/2.)
    elecR.scaling = options.scaleSigR
    elecR.semiinf = options.semiinfR
    DevGF = NEGF.GF(options.TSHS, elecL, elecR, Bulk=options.UseBulk,
                    DeviceAtoms=options.DeviceAtoms,
                    BufferAtoms=options.buffer)
    DevGF.calcGF(options.energy+options.eta*1.0j, options.kpoint[0:2], ispin=options.iSpin,
                 etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
    NEGF.SavedSig.close() # Make sure saved Sigma is written to file

    # Transmission
    print('Transmission Ttot(%.4feV) = %.16f'%(options.energy, N.trace(DevGF.TT).real))

    # Build basis
    options.nspin = DevGF.HS.nspin
    L = options.bufferL
    # Pad lasto with zeroes to enable basis generation...
    lasto = N.zeros((DevGF.HS.nua+L+1,), N.int)
    lasto[L:] = DevGF.HS.lasto
    basis = SIO.BuildBasis(options.fn,
                           options.DeviceAtoms[0]+L,
                           options.DeviceAtoms[1]+L, lasto)
    basis.ii -= L

    # Calculate Eigenchannels
    DevGF.calcEigChan(options.numchan)

    # Compute bond currents?
    if options.kpoint[0] != 0.0 or options.kpoint[1] != 0.0:
        print('Warning: The current implementation of bond currents is only valid for the Gamma point (should be easy to fix)')
        BC = False
    else:
        BC = True

    # Eigenchannels from left
    ECleft = DevGF.ECleft
    for jj in range(options.numchan):
        options.iSide, options.iChan = 0, jj+1
        writeWavefunction(options, geom, basis, ECleft[jj])
        if BC:
            Curr = calcCurrent(options, basis, DevGF.H, ECleft[jj])
            writeCurrent(options, geom, Curr)

    # Calculate eigenchannels from right
    if options.bothsides:
        ECright = DevGF.ECright
        for jj in range(options.numchan):
            options.iSide, options.iChan = 1, jj+1
            writeWavefunction(options, geom, basis, ECright[jj])
            if BC:
                Curr = calcCurrent(options, basis, DevGF.H, ECright[jj])
                writeCurrent(options, geom, Curr)

    # Calculate total "bond currents"
    if BC:
        Curr = -calcCurrent(options, basis, DevGF.H, DevGF.AL)
        options.iChan, options.iSide = 0, 0
        writeCurrent(options, geom, Curr)
        Curr = -calcCurrent(options, basis, DevGF.H, DevGF.AR)
        options.iSide = 1
        writeCurrent(options, geom, Curr)

    # Calculate eigenstates of device Hamiltonian (MPSH)
    if options.MolStates > 0.0:
        try:
            import scipy.linalg as SLA
            ev, es = SLA.eigh(DevGF.H, DevGF.S)
            print('EigenChannels: Eigenvalues (in eV) of computed molecular eigenstates:')
            print(ev)
            # Write eigenvalues to file
            fn = options.DestDir+'/'+options.systemlabel+'.EIGVAL'
            print('EigenChannels: Writing', fn)
            fnfile = open(fn, 'w')
            fnfile.write('# Device region = [%i,%i], units in eV\n'%(options.DeviceFirst, options.DeviceLast))
            for i, val in enumerate(ev):
                fnfile.write('%i %.8f\n'%(i, val))
            fnfile.close()
            # Compute selected eigenstates
            for ii, val in enumerate(ev):
                if N.abs(val) < options.MolStates:
                    fn = options.DestDir+'/'+options.systemlabel+'.S%.3i.E%.3f'%(ii, val)
                    writeWavefunction(options, geom, basis, es[:, ii], fn=fn)
        except:
            print('You need to install scipy to solve the generalized eigenvalue problem')
            print('for the molecular eigenstates in the nonorthogonal basis')

    Log.PrintMainFooter(options)
    return DevGF

def calcWF(options, geom, basis, Y):
    """
    Calculate wavefunction, returns:
    YY : complex wavefunction on regular grid
    dstep : stepsize
    origo : vector
    nx, ny, nz : number of grid points
    """

    xyz = N.array(geom.xyz[options.DeviceAtoms[0]-1:options.DeviceAtoms[1]])

    # Size of cube
    xmin, xmax = min(xyz[:, 0])-5.0, max(xyz[:, 0])+5.0
    ymin, ymax = min(xyz[:, 1])-5.0, max(xyz[:, 1])+5.0
    zmin, zmax = min(xyz[:, 2])-5.0, max(xyz[:, 2])+5.0
    xl, yl, zl = xmax-xmin, ymax-ymin, zmax-zmin
    dx, dy, dz = options.res, options.res, options.res
    nx, ny, nz = int(xl/dx)+1, int(yl/dy)+1, int(zl/dz)+1

    origo = N.array([xmin, ymin, zmin], N.float)

    # Def cube
    YY = N.zeros((nx, ny, nz), N.complex)
    rx = N.array(list(range(nx)), N.float)*dx+origo[0]
    ry = N.array(list(range(ny)), N.float)*dy+origo[1]
    rz = N.array(list(range(nz)), N.float)*dz+origo[2]

    for ii, Yval in enumerate(Y):
        if ii > 0:# and ii%(int(len(Y)/10)) == 0:
            SIO.printDone(ii, len(Y), 'Wavefunction')

        rax, ray, raz = basis.xyz[ii, 0], basis.xyz[ii, 1], basis.xyz[ii, 2]
        # Only calulate in subset
        ixmin, ixmax = int((rax-origo[0]-basis.coff[ii])/dx), \
                       int((rax-origo[0]+basis.coff[ii])/dx)
        iymin, iymax = int((ray-origo[1]-basis.coff[ii])/dy), \
                       int((ray-origo[1]+basis.coff[ii])/dy)
        izmin, izmax = int((raz-origo[2]-basis.coff[ii])/dz), \
                       int((raz-origo[2]+basis.coff[ii])/dz)

        ddx, ddy, ddz = rx[ixmin:ixmax]-rax, ry[iymin:iymax]-ray, rz[izmin:izmax]-raz

        dr = N.sqrt(MM.outerAdd(ddx*ddx, ddy*ddy, ddz*ddz))
        drho = N.sqrt(MM.outerAdd(ddx*ddx, ddy*ddy, 0*ddz))

        imax = (basis.coff[ii]-2*basis.delta[ii])/basis.delta[ii]
        ri = dr/basis.delta[ii]
        ri = N.where(ri < imax, ri, imax)
        ri = ri.astype(N.int)
        costh = MM.outerAdd(0*ddx, 0*ddy, ddz)/dr
        cosfi, sinfi = MM.outerAdd(ddx, 0*ddy, 0*ddz)/drho, MM.outerAdd(0*ddx, ddy, 0*ddz)/drho

        # Numpy has changed the choose function to crap!
        RR = N.take(basis.orb[ii], ri)

        # Calculate spherical harmonics
        l = basis.L[ii]
        m = basis.M[ii]
        if l == 3:
            print('f-shell : l=%i, m=%i (NOT TESTED!!)'%(l, m))
        thisSphHar = MM.sphericalHarmonics(l, m, costh, sinfi, cosfi)

        YY[ixmin:ixmax, iymin:iymax, izmin:izmax] = YY[ixmin:ixmax, iymin:iymax, izmin:izmax]+\
                                                    RR*thisSphHar*Yval

    print("Wave function norm on real space grid:", N.sum(YY.conjugate()*YY)*dx*dy*dz)

    return YY, options.res, origo, nx, ny, nz

###################### Current densities #########################


def calcCurrent(options, basis, H, Y):
    """
    Calculate current density in atomic bonds
    Y : complex scattering state or
    Y : A_l or A_r! (for total current)
    """

    if isinstance(Y, MM.SpectralMatrix):
        Y = MM.mm(Y.L, Y.R)
    NN = len(H)
    NN2 = options.DeviceAtoms[1]-options.DeviceAtoms[0]+1
    Curr = N.zeros((NN2, NN2), N.float)

    if len(Y.shape) == 2:
        for ii in range(NN):
            a1 = basis.ii[ii]-options.DeviceAtoms[0]
            for jj in range(NN):
                a2 = basis.ii[jj]-options.DeviceAtoms[0]
                tmp = H[jj, ii]*Y[ii, jj]/2/N.pi
                # Note that taking the imaginary part is only the valid
                # expression for Gamma point calculations
                Curr[a1, a2] = Curr[a1, a2]+4*N.pi*tmp.imag
    else:
        for ii in range(NN):
            a1 = basis.ii[ii]-options.DeviceAtoms[0]
            for jj in range(NN):
                a2 = basis.ii[jj]-options.DeviceAtoms[0]
                tmp = H[ii, jj]*N.conjugate(Y[ii])*Y[jj]
                Curr[a1, a2] = Curr[a1, a2]+4*N.pi*tmp.imag

    return Curr


def writeCurrent(options, geom, Curr):
    """
    Write bond currents to file
    fn : Filename
    Curr : Bond current matrix
    """
    fn = fileName(options)
    xyz = N.array(geom.xyz[options.DeviceAtoms[0]-1:options.DeviceAtoms[1]])
    atomnum = geom.anr[options.DeviceAtoms[0]-1:options.DeviceAtoms[1]]

    foC = open(fn+'.curr', 'w')
    foC.write('%i\n'%(options.DeviceAtoms[1]-options.DeviceAtoms[0]+1))

    for ii in range(len(xyz)):
        foC.write('%i %e %e %e\n'%(atomnum[ii], xyz[ii, 0], xyz[ii, 1], xyz[ii, 2]))

    for ii in range(len(Curr)):
        for jj in range(len(Curr)):
            foC.write('%e '%(Curr[ii, jj]))
        foC.write('\n')
    foC.close()
    # Current vector for atom i
    Curr2 = N.zeros(geom.xyz.shape)
    for i in range(len(Curr)):
        for j in range(len(Curr)):
            if i != j:
                R = N.zeros(3, N.float)
                for k in range(-1, 2): # Loop over neighbors
                    for l in range(-1, 2):
                        r = xyz[i]-xyz[j]+k*geom.pbc[0]+l*geom.pbc[1]
                        dr = N.dot(r, r)**.5 # Norm
                        R += r*N.exp(-dr) # Exponential weighting
                R = R/N.dot(R, R)**.5
                Curr2[options.DeviceAtoms[0]-1+i] += Curr[i, j]*R/2
    SIO.WriteAXSFFiles(fn+'.curr.AXSF', [geom], [Curr2])

########################################################


def writenetcdf(geom, fn, YY, nx, ny, nz, origo, dstep):
    """
    THF: Write eigenchannels to netcdf format
    """
    ncfile = NC4.Dataset(fn, 'w')
    ncfile.createDimension('nx', nx)
    ncfile.createDimension('ny', ny)
    ncfile.createDimension('nz', nz)
    ncfile.createDimension('natoms', len(geom.xyz))
    ncfile.createDimension('naxes', 3)
    ncfile.createDimension('number', 1)
    #ncfile.createDimension('pair',2)

    # Grid
    #grid = ncfile.createVariable('grid','d',('naxes','pair'))
    #tmp  = [origo,[dstep,dstep,dstep]]
    #grid[:] = N.transpose(N.array(tmp))

    # Fields
    varRe = ncfile.createVariable('Re-Psi', 'd', ('nx', 'ny', 'nz'))
    varRe[:] = YY.real

    varIm = ncfile.createVariable('Im-Psi', 'd', ('nx', 'ny', 'nz'))
    varIm[:] = YY.imag

    varAbsSq = ncfile.createVariable('Abs-sqr-Psi', 'd', ('nx', 'ny', 'nz'))
    varAbsSq[:] = N.absolute(N.square(YY))

    vardstep = ncfile.createVariable('dstep', 'd', ('number', ))
    vardstep[:] = dstep
    vardstep.units = 'Ang'

    varorigo = ncfile.createVariable('origo', 'd', ('naxes', ))
    varorigo[:] = origo
    varorigo.units = 'Ang'

    # Old stuff for OpenDX visualization
    vargeom = ncfile.createVariable('xyz-dx', 'f', ('natoms', 'naxes'))
    tmp = []
    for i in range(len(geom.xyz)):
        tmp2 = (geom.xyz[i]-origo)/dstep
        tmp.append([float(tmp2[0]),
                    float(tmp2[1]),
                    float(tmp2[2])])
    vargeom[:] = tmp
    vargeom.units = '(xyz-origo)/dstep'

    # Absolute coordiates of atoms and cell
    vargeom = ncfile.createVariable('xyz', 'f', ('natoms', 'naxes'))
    vargeom[:] = geom.xyz
    vargeom.units = 'Ang'

    varcell = ncfile.createVariable('cell', 'f', ('naxes', 'naxes'))
    varcell[:] = geom.pbc
    varcell.units = 'Ang'

    varanr = ncfile.createVariable('anr', 'i', ('natoms',))
    varanr[:] = N.array(geom.anr, N.int32)

    # Set attributes
    #setattr(varanr, 'field', 'anr')
    #setattr(varanr, 'positions', 'xyz')

    #setattr(varRe,'field','Re-Psi')
    #setattr(varRe,'positions','grid, compact')

    #setattr(varIm,'field','Im-Psi')
    #setattr(varIm,'positions','grid, compact')

    ncfile.close()

################# Write cube file ######################


def writecube(geom, fn, YY, nx, ny, nz, origo, dstep):
    """
    Write wavefunction to cube file.
    """
    xyz = N.array(geom.xyz)
    anr = geom.anr

    foR = open(fn, 'w')
    foR.write('Eigenchannel wavefunction\n%s\n'%fn)
    foR.write('%i %f %f %f\n'% (len(xyz), origo[0]/PC.Bohr2Ang, origo[1]/PC.Bohr2Ang, origo[2]/PC.Bohr2Ang))
    foR.write('%i %f %f %f\n'% (nx, dstep/PC.Bohr2Ang, 0.0, 0.0))
    foR.write('%i %f %f %f\n'% (ny, 0.0, dstep/PC.Bohr2Ang, 0.0))
    foR.write('%i %f %f %f\n'% (nz, 0.0, 0.0, dstep/PC.Bohr2Ang))
    # Write atom coordinates
    for ii in range(len(xyz)):
        foR.write('%i %f '% (anr[ii], 0.0))
        tmp = xyz[ii, :]/PC.Bohr2Ang
        foR.write('%f %f %f\n'% (tmp[0], tmp[1], tmp[2]))
    # Write wavefunction
    YYY = YY.real*(PC.Bohr2Ang**(3.0/2.0))
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                foR.write('%1.3e\n' % YYY[ix, iy, iz])
    foR.close()

################# Write wave function in macu format ######################


def writemacubin(fn, YY, nx, ny, nz, origo, dstep):
    """
    Write molekel binary format
    """

    fo = open(fn, 'w')
    fo.write(struct.pack('i', 36))
    xmin, xmax = origo[0], origo[0]+dstep*(nx-1)
    ymin, ymax = origo[1], origo[1]+dstep*(ny-1)
    zmin, zmax = origo[2], origo[2]+dstep*(nz-1)
    fo.write(struct.pack('6f4i',
                         xmin, xmax, ymin, ymax, zmin, zmax, nx, ny, nz, nx*ny*nz))

    for ii in range(nz):
        mlklbin = struct.pack('i', nx*ny*4)
        for kk in range(ny):
            for jj in range(nx):
                mlklbin += struct.pack('f', YY[jj, kk, ii])
        mlklbin += struct.pack('i', nx*ny*4)
        fo.write(mlklbin)

    fo.close()

################ Write wave function in XSF format ##################################


def writeXSF(geom, fn, YY, nx, ny, nz, origo, dstep):
    """
    Write XSF datagrid for XCrysden
    """
    fo = open(fn, 'w')
    speciesnumber = geom.snr
    atomnumber = geom.anr
    xyz = geom.xyz
    xmin, xmax = origo[0], origo[0]+dstep*(nx-1)
    ymin, ymax = origo[1], origo[1]+dstep*(ny-1)
    zmin, zmax = origo[2], origo[2]+dstep*(nz-1)

    fo.write(' ATOMS\n')
    numberOfAtoms = len(speciesnumber)

    # Write out the position for the atoms
    for i in range(numberOfAtoms):
        line = '   %i  ' %atomnumber[i]
        for j in range(3):
            line += ('%.9f'%xyz[i][j]).rjust(16)
        fo.write(line+'\n')

    #Write the datagrid
    fo.write('BEGIN_BLOCK_DATAGRID_3D\n')
    fo.write(' DATA_from:Inelastica\n')
    # Real part
    fo.write(' BEGIN_DATAGRID_3D_REAL\n')
    fo.write('    %3.0i    %3.0i    %3.0i\n'%(nx, ny, nz))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (origo[0], origo[1], origo[2]))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (xmax-xmin, 0.0000, 0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000, ymax-ymin, 0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000, 0.0000, zmax-zmin))
    data = []
    for ii in range(nz):
        for kk in range(ny):
            for jj in range(nx):
                data.append(YY.real[jj, kk, ii])
    for iii in range((nx*ny*nz)):
        if (iii+1)%6 == 0:
            fo.write('  %1.5E\n'% (data[iii]))
        else:
            fo.write('  %1.5E'% (data[iii]))
    fo.write('\n END_DATAGRID_3D\n')
    # Imaginary part
    fo.write(' BEGIN_DATAGRID_3D_IMAG\n')
    fo.write('    %3.0i    %3.0i    %3.0i\n'%(nx, ny, nz))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (origo[0], origo[1], origo[2]))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (xmax-xmin, 0.0000, 0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000, ymax-ymin, 0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000, 0.0000, zmax-zmin))
    data = []
    for ii in range(nz):
        for kk in range(ny):
            for jj in range(nx):
                data.append(YY.imag[jj, kk, ii])
    for iii in range((nx*ny*nz)):
        if (iii+1)%6 == 0:
            fo.write('  %1.5E\n'% (data[iii]))
        else:
            fo.write('  %1.5E'% (data[iii]))
    fo.write('\n END_DATAGRID_3D\n')
    # Absolute square
    fo.write(' BEGIN_DATAGRID_3D_ABS\n')
    fo.write('    %3.0i    %3.0i    %3.0i\n'%(nx, ny, nz))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (origo[0], origo[1], origo[2]))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (xmax-xmin, 0.0000, 0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000, ymax-ymin, 0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000, 0.0000, zmax-zmin))
    data = []
    YYA2 = N.absolute(N.square(YY))
    for ii in range(nz):
        for kk in range(ny):
            for jj in range(nx):
                data.append(YYA2[jj, kk, ii])
    for iii in range((nx*ny*nz)):
        if (iii+1)%6 == 0:
            fo.write('  %1.5E\n'% (data[iii]))
        else:
            fo.write('  %1.5E'% (data[iii]))
    fo.write('\n END_DATAGRID_3D\n')
    fo.write('END_BLOCK_DATAGRID_3D')
    fo.close()


########################################################
def writeWavefunction(options, geom, basis, Y, fn=None):
    """
    Writes wave function to the specified output type
    Y: vector for wavefunction

    """
    if fn == None:
        fn = fileName(options)
    print('Eigenchannels.writeWavefunction: Writing', fn)
    # Rotate in complex space
    max_amp = -1.0
    phase = 1.0+0.0j

    for Ykk in Y:
        if abs(Ykk) > max_amp:
            max_amp = abs(Ykk)
            phase = Ykk/max_amp
    Y = Y/phase

    foT = open(fn+'.abs.txt', 'w')
    foT.write('index Atom nr   M   L abs(Y)         re(Y)        im(Y)\n')
    for ii, Yval in enumerate(Y):
        foT.write('%5.1i %3.0i %3.0i %3.1i %3.1i %1.8f %15.5e %12.5e\n'%
        (ii, basis.ii[ii], basis.atomnum[ii], basis.M[ii], basis.L[ii], abs(Yval), Yval.real, Yval.imag))
    foT.close()

    YY, dstep, origo, nx, ny, nz = calcWF(options, geom, basis, Y)

    # Write wave function in specified file format
    if options.format.lower() == 'macu':
        writemacubin(fn+'.Re.macu', YY.real, nx, ny, nz, origo, dstep)
        writemacubin(fn+'.Im.macu', YY.imag, nx, ny, nz, origo, dstep)
        writemacubin(fn+'.Abs.macu', N.absolute(N.square(YY)), nx, ny, nz, origo, dstep)

    if options.format.lower() == 'cube':
        writecube(geom, fn+'.Re.cube', YY.real, nx, ny, nz, origo, dstep)
        writecube(geom, fn+'.Im.cube', YY.imag, nx, ny, nz, origo, dstep)
        writecube(geom, fn+'.Abs.cube', N.absolute(N.square(YY)), nx, ny, nz, origo, dstep)

    if options.format.lower() == 'xsf':
        writeXSF(geom, fn+'.XSF', YY, nx, ny, nz, origo, dstep)

    if options.format.lower() == 'nc':
        writenetcdf(geom, fn+'.nc', YY, nx, ny, nz, origo, dstep)


################# File names ##################################
def fileName(options):
    systemlabel = options.systemlabel
    # Generate filename '.EC.{1,Tot}{L,R,,In,Out}[UP][Ef=1.0].'
    if options.iChan == 0:
        fn = systemlabel+'.EC.Tot%s'%(['L', 'R', '', 'In', 'Out'][options.iSide])
    else:
        fn = systemlabel+'.EC.%i%s'%(options.iChan, ['L', 'R', '', 'In', 'Out'][options.iSide])
    if options.nspin == 2:
        fn += ['UP', 'DOWN'][options.iSpin]
    if options.energy != 0.0:
        fn += '_E%.3f'%options.energy
    if options.kpoint[0] != 0.0 or options.kpoint[1] != 0.0:
        fn += '_kx%.3f_ky%.3f'%(options.kpoint[0], options.kpoint[1])
    return options.DestDir+'/'+fn
