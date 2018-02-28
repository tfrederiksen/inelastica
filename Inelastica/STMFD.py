"""

STMFD (:mod:`Inelastica.STMFD`)
===============================

.. currentmodule:: Inelastica.STMFD

"""

import numpy               as N
import scipy.linalg        as SLA
import netCDF4             as NC
from   scipy.sparse.linalg import isolve
from   scipy               import sparse
from   scipy               import interpolate
import time
import glob
import Inelastica.io.netcdf as wNC
import Inelastica.physics.constants as PC
import Inelastica.io.siesta as SIO


def main(options, kpoint, ikpoint):
    Ef = SIO.HS(options.systemlabel+'.TSHS').ef/PC.Rydberg2eV
    kpt = ikpoint
    path = './'+options.DestDir+'/'
    pathkpt = './'+options.DestDir+'/'+str(kpt)+'/'
    print 'k-point: '+str(ikpoint)+'/'
    print '(k1,k2): ('+str(kpoint[0])+','+str(kpoint[1])+')'
    print 'Fermi energy ('+str(options.systemlabel)+'.TSHS):', N.round(Ef, 4), 'Ry =', N.round(Ef*PC.Rydberg2eV, 4), 'eV'

    #Determine substrate layers and tip height
    posZMol, posZTip = LayersAndTipheight(options, kpoint, ikpoint)

    #Read total potential, loc-basis states, lattice constants etc.
    tmp    = readDFT(options, kpoint, pathkpt, posZMol, posZTip)
    Subwfs, Tipwfs, Vsub, Vtip = tmp[0], tmp[1], tmp[2], tmp[3]
    scSize, ucSize, Max, dS, theta = tmp[4], tmp[5], tmp[6], tmp[7], tmp[8]
    SubChans, TipChans, SubPot, TipPot, SubRho, TipRho, MeshCutoff = tmp[9], tmp[10], tmp[11], tmp[12], tmp[13], tmp[14], tmp[15]

    Nx, Ny, Nz = scSize[0], scSize[1], scSize[2]
    NN = Nx*Ny*Nz
    a1, a2, a3 = ucSize[0], ucSize[1], ucSize[2]
    if options.samplingscale != 1:
        print '\nReal-space sampling in xy plane coarser by factor', options.samplingscale, '...'
        tmp = sampling(options, Nx, Ny, Nz, Subwfs, Tipwfs, SubChans, TipChans, SubPot, TipPot, SubRho, TipRho)
        Nx, Ny, Nz  = tmp[0], tmp[1], tmp[2]
        NN        = Nx*Ny*Nz
        a1, a2     = a1*options.samplingscale, a2*options.samplingscale
        ucSize    = [a1, a2, a3]
        Subwfs, Tipwfs, Vsub, Vtip, SubRho, TipRho = tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8]
        print 'New lattice in xy plane obtained (corresponding to '+str(MeshCutoff/options.samplingscale**2)+' Ry):'
        print '          [Nx Ny Nz] = ['+str(Nx), str(Ny), str(Nz)+']'
        print '          [a1 a2 a3] = ['+str(N.round(a1*PC.Bohr2Ang, 4)), str(N.round(a2*PC.Bohr2Ang, 4)), str(N.round(a3*PC.Bohr2Ang, 4))+'] Ang'
    else:
        print '\nReal-space sampling omitted (same grid as TranSiesta).\n'

    #Substrate part
    print '\nCalculation from substrate side:'
    WFs, Chans, rho = Subwfs, SubChans, SubRho
    inda, indb = reindexing(options, Nx, Ny, Nz, NN, WFs, Chans, rho, ucSize, theta)
    Pot       = Vsub
    Ham       = Hamiltonian(options, a1, a2, a3, Nx, Ny, Nz, NN, Pot, theta, kpoint)
    Tau, Hb    = SplitHam(Ham, inda, indb)
    ChansL    = Chans
    propSubModes = LinearSolve(options, WFs, inda, indb, Ef, Tau, Hb, NN, Nx, Ny, Nz, Chans)

    #Tip part
    print '\nCalculation from tip side:'
    WFs, Chans, rho = Tipwfs, TipChans, TipRho
    inda, indb = reindexing(options, Nx, Ny, Nz, NN, WFs, Chans, rho, ucSize, theta)
    Pot       = Vtip
    Ham       = Hamiltonian(options, a1, a2, a3, Nx, Ny, Nz, NN, Pot, theta, kpoint)
    Tau, Hb    = SplitHam(Ham, inda, indb)
    ChansR    = Chans
    propTipModes = LinearSolve(options, WFs, inda, indb, Ef, Tau, Hb, NN, Nx, Ny, Nz, Chans)

    #Combine sub and tip to find conductance
    STMimage, Currmat = Current(options, propSubModes, propTipModes, Nx, Ny, Nz, Max,\
                               ucSize, ChansL, ChansR, dS, kpoint, pathkpt, ikpoint)
    STMimage = STMimage[:Nx, :Ny]


def LayersAndTipheight(options, kpoint, ikpoint):
    tmp  = NC.Dataset('./'+options.DestDir+'/'+str(ikpoint)+'/'+options.systemlabel+'.AL0.nc', 'r')
    xyzSupercell = N.array(tmp.variables['xyz'][:], N.float)
    xyz = xyzSupercell.copy()
    noAtoms = len(xyzSupercell[:, 0])
    for ii in range(noAtoms):
        xyz[ii, 2] = xyzSupercell[ii, 2]-xyzSupercell[0, 2]
    atomlayer1 = N.zeros(len(xyz))
    #First substrate layer
    for ii in range(len(xyz)-1):
        atomlayer1[ii] += xyz[ii][2]
        if N.abs(N.sum(atomlayer1)/(ii+1)-atomlayer1[ii])>.5:
            break
    Layer1 = ii
    print '\nThe atomic layers seem to consist of '+str(ii)+' atoms.'
    zSurfLayer1 = N.sum(atomlayer1[0:ii])/ii
    #Second substrate layer?
    atomlayer2 = N.zeros(len(xyz))
    for jj in range(ii, len(xyz)-1):
        atomlayer2[jj] += xyz[jj][2]
        if N.abs(N.sum(atomlayer2[:])/(jj-ii+1)-atomlayer2[jj])>.5:
            break
    Layer2 = jj
    #Third!?
    atomlayer3 = N.zeros(len(xyz))
    for kk in range(jj, len(xyz)-1):
        atomlayer3[kk] += xyz[kk][2]
        if N.abs(N.sum(atomlayer3[:])/(kk-jj+1)-atomlayer3[kk])>.5:
            break
    Layer3 = kk
    if Layer2-Layer1<Layer1:
        print 'The substrate seems to consist a single atomic layer.'
    elif Layer3-Layer2<Layer1:
        print 'The substrate seems to consist of two atomic layers.'
    else:
        print 'The substrate seems to consist of at three (or more) atomic layers.'
        print '(One or two layers are enough. Check your -F and -L flags!)'
    for ii in range(len(xyz)-1):
        if xyz[ii+1][2]-xyz[ii][2]>5/PC.Bohr2Ang: #Minimum gap: 5 Ang!
            TipHeightMol = N.round(N.abs(xyz[ii+1][2]-xyz[ii][2])*PC.Bohr2Ang, 3)
            print 'Vacuum gap along z appears to be '+str(TipHeightMol)+' Ang'
            break
    Molidx = ii; Tipidx = ii+1
    posZMol = xyz[Molidx][2]; posZTip = xyz[Tipidx][2]
    TipHeightLayer1 = N.round(N.abs(posZTip-zSurfLayer1)*PC.Bohr2Ang, 3)
    if Layer2-Layer1<Layer1:
        print 'Substrate surface  --- tip-apex distance: '+str(TipHeightLayer1)+' Ang'
    if Layer2-Layer1 == Layer1:
        zSurfLayer2 = N.sum(atomlayer2[Layer1:Layer2])/(Layer2-Layer1)
        TipHeightLayer2 = N.round(N.abs(posZTip-zSurfLayer2)*PC.Bohr2Ang, 3)
        print 'Substrate surface (1st) --- tip-apex distance: '+str(TipHeightLayer1)+' Ang'
        print 'Substrate surface (2nd) --- tip-apex distance: '+str(TipHeightLayer2)+' Ang'
    if Layer3-Layer2 == Layer1:
        zSurfLayer3 = N.sum(atomlayer3[Layer2:Layer3])/(Layer3-Layer2)
        TipHeightLayer3 = N.round(N.abs(posZTip-zSurfLayer3)*PC.Bohr2Ang, 3)
        print 'Substrate surface (3rd) --- tip-apex distance: '+str(TipHeightLayer3)+' Ang'
    return posZMol, posZTip


def readDFT(options, kpt, pathkpt, posZMol, posZTip):
    file = NC.Dataset('TotalPotential.grid.nc', 'r')
    print '\nReading potential from:    TotalPotential.grid.nc'
    pot = N.array(file.variables['gridfunc'][:], N.float)[0]

    print 'Reading DFT density from:  Rho.grid.nc'
    rho = N.array(file.variables['gridfunc'][:], N.float)[0]

    #Import supercell grid from one localized-basis state
    tmp  = NC.Dataset(pathkpt+options.systemlabel+'.AL0.nc', 'r')
    tmp2 = N.array(tmp.variables['Re-Psi'][:], N.float)
    dim  = N.shape(tmp2)
    Nx, Ny, Nz = dim[0], dim[1], dim[2];

    SubChans = N.shape(glob.glob(pathkpt+options.systemlabel+'.AL*.nc'))[0]
    TipChans = N.shape(glob.glob(pathkpt+options.systemlabel+'.AR*.nc'))[0]

    Subwfs = N.zeros((SubChans, dim[0], dim[1], dim[2]), N.complex)
    Tipwfs = N.zeros((TipChans, dim[0], dim[1], dim[2]), N.complex)

    print '\nReading localized-basis wave functions from substrate side...'
    for ii in range(SubChans):
        file = NC.Dataset(pathkpt+options.systemlabel+'.AL'+N.str(ii)+'.nc', 'r')
        print options.systemlabel+'.AL'+N.str(ii)+'.nc'
        Re_AL = N.array(file.variables['Re-Psi'][:], N.float)
        Im_AL = N.array(file.variables['Im-Psi'][:], N.float)
        Subwfs[ii, :, :, :] = Re_AL+1j*Im_AL
        file.close()

    print 'Reading localized-basis wave functions from tip side ...'
    for ii in range(TipChans):
        file = NC.Dataset(pathkpt+options.systemlabel+'.AR'+N.str(ii)+'.nc', 'r')
        print options.systemlabel+'.AR'+N.str(ii)+'.nc'
        Re_AR = N.array(file.variables['Re-Psi'][:], N.float)
        Im_AR = N.array(file.variables['Im-Psi'][:], N.float)
        Tipwfs[ii, :, :, ::-1] = Re_AR+1j*Im_AR
        file.close()
    file = NC.Dataset(pathkpt+options.systemlabel+'.AL0.nc', 'r')
    steps = N.array(file.variables['steps'][:], N.float)
    dS = SLA.norm(N.cross(steps[0], steps[1]))
    theta = N.arccos(N.dot(steps[0], steps[1])/(SLA.norm(steps[0])*SLA.norm(steps[1])))
    a1, a2, a3 = SLA.norm(steps[0]), SLA.norm(steps[1]), SLA.norm(steps[2])
    avec = [a1, a2, a3]
    orig = N.array(file.variables['origin'][:], N.float)[2]
    Nzi  = N.int(orig/a3)
    Nzf  = Nzi+Nz
    PotDevice = pot[Nzi:Nzf, :, :]
    SubPot = N.transpose(PotDevice, [2, 1, 0])
    TipPot = SubPot[:, :, ::-1].copy()
    PotDev = SubPot.copy()
    SubRho = N.transpose(rho[Nzi:Nzf, :, :], [2, 1, 0])
    TipRho = SubRho[:, :, ::-1].copy()
    print '\nPotential is cut out to fit the real-space projected localized-basis grid'

    try:
        for line in open('RUN.fdf').readlines():
            if 'MeshCutoff' in line:
                break
        eval(line.split()[1])
    except:
        for line in open('Default.fdf').readlines():
            if 'MeshCutoff' in line:
                break
        eval(line.split()[1])

    MeshCutoff = eval(line.split()[1])

    print 'Lattice obtained by using DFT energy cutoff '+str(MeshCutoff)+' Ry:'
    print '          [Nx Ny Nz] = ['+str(Nx), str(Ny), str(Nz)+']'
    print '          [a1 a2 a3] = ['+str(N.round(a1*PC.Bohr2Ang, 4)), str(N.round(a2*PC.Bohr2Ang, 4)), str(N.round(a3*PC.Bohr2Ang, 4))+'] Ang'
    print 'Angle spanned by the lateral unit vectors a1 and a2: '+str(N.round(theta, 5))+' rad (Pi/'+str(N.round(N.pi/theta, 3))+')'
    print '\nPositioning the separation plane on which the wave functions are evaluated:'
    iShiftSep = N.int(options.ShiftSeparationPlane/PC.Bohr2Ang/a3)
    usedPlane = iShiftSep*a3*PC.Bohr2Ang
    if options.ShiftSeparationPlane != 0:
        print 'Requested separation-plane shift: '+str(options.ShiftSeparationPlane)+' Ang. Used: '+str(N.round(usedPlane, 4))+' Ang.'
    tmp = [N.sum(PotDev[:, :, ii]) for ii in range(Nz)]
    MaxPot = N.max(tmp)
    MaxIdx = tmp.index(MaxPot)+iShiftSep
    print 'Maximum potential at slice', MaxIdx, 'from left contact ('+str(Nz-MaxIdx)+' from right)'
    posZSepSurf = MaxIdx*a3
    print 'Position, separation plane: '+str(N.round(N.abs(posZSepSurf-posZMol)*PC.Bohr2Ang, 3))+' Ang above uppermost substrate atom'
    print '                            '+str(N.round(N.abs(posZTip-posZSepSurf)*PC.Bohr2Ang, 3))+' Ang below the tip apex'

    if N.abs(posZMol-posZSepSurf) < 1./PC.Bohr2Ang:
        print 'Warning: The separation plane seems to be too close to the substrate.'
        print '(Should be at least 2 Ang. Change with flag --ssp)'
    if N.abs(posZSepSurf-posZTip) < 1./PC.Bohr2Ang:
        print 'Warning: The separation plane seems to be too close to the tip.'
        print '(Should be at least 2 Ang. Change with flag --ssp)'

    VacPot = N.mean(SubPot[:, :, MaxIdx])
    print '\nVacuum potential at (and away from) the separation plane:', N.round(VacPot, 4), 'Ry =', N.round(VacPot*PC.Rydberg2eV, 4), 'eV'
    SubPot[:, :, MaxIdx:] = VacPot
    TipPot[:, :, Nz-1-MaxIdx:] = VacPot
    #Flatten potential to fit FD Hamiltonian
    Vsub = [SubPot[ii, jj, kk] for kk in range(Nz) \
               for jj in range(Ny) for ii in range(Nx)]
    Vtip = [TipPot[ii, jj, kk] for kk in range(Nz) \
               for jj in range(Ny) for ii in range(Nx)]
    SubRho[:, :, MaxIdx:] = 0
    TipRho[:, :, Nz-1-MaxIdx:] = 0
    return Subwfs, Tipwfs, Vsub, Vtip, dim, avec, MaxIdx, dS, theta, SubChans, TipChans, SubPot, TipPot, SubRho, TipRho, MeshCutoff


def sampling(options, Nx, Ny, Nz, Subwfs, Tipwfs, SubChans, TipChans, SubPot, TipPot, SubRho, TipRho):
    ssc = options.samplingscale
    interMethod = 'linear'
    NX, NY, NZ = N.int(N.floor(Nx/ssc)), N.int(N.floor(Ny/ssc)), Nz
    x, y = N.arange(0, Nx), N.arange(0, Ny)
    xnew, ynew = N.arange(0, Nx, ssc), N.arange(0, Ny, ssc)

    ReNewSubwfs, ImNewSubwfs  = N.zeros((SubChans, NX, NY, NZ)), N.zeros((SubChans, NX, NY, NZ))
    for imode in range(SubChans):
        for iz in range(NZ):
            tmp1 = N.real(Subwfs[imode, :, :, iz]).T
            f1 = interpolate.interp2d(x, y, tmp1, kind=interMethod)
            ReNewSubwfs[imode, :, :, iz] = f1(xnew, ynew).T
            tmp2 = N.imag(Subwfs[imode, :, :, iz]).T
            f2 = interpolate.interp2d(x, y, tmp2, kind=interMethod)
            ImNewSubwfs[imode, :, :, iz] = f2(xnew, ynew).T
    NewSubwfs = ReNewSubwfs+1j*ImNewSubwfs

    ReNewTipwfs, ImNewTipwfs  = N.zeros((TipChans, NX, NY, NZ)), N.zeros((TipChans, NX, NY, NZ))
    for imode in range(TipChans):
        for iz in range(NZ):
            tmp1 = N.real(Tipwfs[imode, :, :, iz]).T
            f1 = interpolate.interp2d(x, y, tmp1, kind=interMethod)
            ReNewTipwfs[imode, :, :, iz] = f1(xnew, ynew).T
            tmp2 = N.imag(Tipwfs[imode, :, :, iz]).T
            f2 = interpolate.interp2d(x, y, tmp2, kind=interMethod)
            ImNewTipwfs[imode, :, :, iz] = f2(xnew, ynew).T
    NewTipwfs = ReNewTipwfs+1j*ImNewTipwfs
    NewSubPot, NewTipPot = N.zeros((NX, NY, NZ)), N.zeros((NX, NY, NZ))
    for iz in range(NZ):
        tmp1 = SubPot[:, :, iz].T
        f1 = interpolate.interp2d(x, y, tmp1, kind=interMethod)
        NewSubPot[:, :, iz] = f1(xnew, ynew).T
        tmp2 = TipPot[:, :, iz].T
        f2 = interpolate.interp2d(x, y, tmp2, kind=interMethod)
        NewTipPot[:, :, iz] = f2(xnew, ynew).T
    NewSubRho, NewTipRho = N.zeros((NX, NY, NZ)), N.zeros((NX, NY, NZ))
    for iz in range(NZ):
        tmp1 = SubRho[:, :, iz].T
        f1 = interpolate.interp2d(x, y, tmp1, kind=interMethod)
        NewSubRho[:, :, iz] = f1(xnew, ynew).T
        tmp2 = TipRho[:, :, iz].T
        f2 = interpolate.interp2d(x, y, tmp2, kind=interMethod)
        NewTipRho[:, :, iz] = f2(xnew, ynew).T
    #Flatten new potential to fit discrete Hamiltonian
    NewVSub = [NewSubPot[ii, jj, kk] for kk in range(NZ) \
               for jj in range(NY) for ii in range(NX)]
    NewVTip = [NewTipPot[ii, jj, kk] for kk in range(NZ) \
               for jj in range(NY) for ii in range(NX)]
    return NX, NY, NZ, NewSubwfs, NewTipwfs, NewVSub, NewVTip, NewSubRho, NewTipRho


def reindexing(options, Nx, Ny, Nz, NN, WFs, Chans, rho, ucSize, theta):
    rhoindx = N.array([rho[ii, jj, kk] for kk in range(Nz) for jj in range(Ny) for ii in range(Nx)])
    inda = N.where(rhoindx>options.rhoiso)[0]
    tmp  = set(inda)
    indb = [ii for ii in N.arange(NN) if ii not in tmp]
    return inda, indb


def Hamiltonian(options, a1, a2, a3, Nx, Ny, Nz, NN, Pot, theta, kpoint):
    print 'Creating the finite-difference Laplacian including the total DFT potential'
    kx = -kpoint[0]*2.0*N.pi
    ky = -kpoint[1]*2.0*N.pi
    #Determine the coefficients of nabla^2, A dx^2 + B dy^2 + C dxdy: (C!=0 in skew system)
    mat = SLA.inv([[1, 0], [N.cos(theta), N.sin(theta)]])
    mixX  = mat[0, 0]**2 + mat[1, 0]**2
    mixY  = mat[0, 1]**2 + mat[1, 1]**2
    mixXY = 2*mat[0, 0]*mat[0, 1] + 2*mat[1, 0]*mat[1, 1]
    mixX, mixY, mixXY = N.round(mixX, 8), N.round(mixY, 8), N.round(mixXY, 8)
    print 'nabla^2_{xy} = '+str(mixX)+'d^2/dx^2+'+str(mixY)+'d^2/dy^2+'+str(mixXY)+'d^2/dxdy'
    #print 'd^2x:',mixX,',d2^y:',mixY,',dxdy:',mixXY
    bx = mixX/a1**2; by = mixY/a2**2; bz = 1./a3**2;
    bands = 23
    D2    = N.zeros((bands, NN), N.complex)

    #tau z
    D2[0, Nx*Ny:NN] = [-bz for ii in range(Nx*Ny, NN)]
    #mix pbc xy, upper1
    for ii in range(Nx*Ny-1, NN, Nx*Ny):
        D2[1, ii] = -mixXY/(4*a1*a2)*N.exp(1j*(kx+ky))
    #mix pbc xy, upper2
    for ii in range(Nx*Ny-Nx+1, NN, Nx*Ny):
        D2[2, ii:(ii+Nx-1)] = mixXY/(4*a1*a2)*N.exp(1j*ky)
    #pbc y
    for ii in range(Nx*Ny-Nx, NN, Nx*Ny):
        D2[3, ii:(ii+Nx)] = -by*N.exp(1j*ky)
    #mix pbc xy, lower2
    for ii in range(Nx*Ny-Nx, NN, Nx*Ny):
        D2[4, ii:(ii+Nx-1)] = -mixXY/(4*a1*a2)*N.exp(1j*ky)
    #mix pbc xy, lower1
    for ii in range(Nx*Ny-Nx, NN, Nx*Ny):
        D2[5, ii] = mixXY/(4*a1*a2)*N.exp(1j*(ky-kx))
    #mix pbc x, upper
    for ii in range(2*Nx-1, NN, Nx):
        D2[6, ii] = mixXY/(4*a1*a2)*N.exp(1j*kx)
    for ii in range(1, Nz):
        D2[6, Nx-1+Nx*Ny*ii] = 0
    #mix xy, upper
    for ii in range(Nx+1, NN, Nx):
        D2[7, ii:(ii+Nx-1)] = -mixXY/(4*a1*a2)
    for ii in range(1, Nz):
        D2[7, 1+Nx*Ny*ii:Nx*Ny*ii+Nx] = 0
    #tau y
    for ii in range(Nx, NN):
        D2[8, ii] = -by
    for ii in range(1, Nz):
        D2[8, Nx*Ny*ii:Nx+Nx*Ny*ii] = 0
    #mix xy, lower
    for ii in range(Nx, NN, Nx):
        D2[9, ii:(ii+Nx-1)] = mixXY/(4*a1*a2)
    for ii in range(1, Nz):
        D2[9, Nx*Ny*ii:Nx*Ny*ii+Nx-1] = 0
    #pbc x, upper
    for ii in range(Nx-1, NN, Nx):
        D2[9, ii] = -bx*N.exp(1j*kx)
    #mix pbc x, lower
    for ii in range(Nx, NN, Nx):
        D2[10, ii] = -mixXY/(4*a1*a2)*N.exp(-1j*kx) #<- Minus sign!
    for ii in range(1, Nz):
        D2[10, Nx*Ny*ii] = 0
    #tau x, upper
    for ii in range(1, NN, Nx):
        D2[10, ii:(ii+Nx-1)] = -bx
    #diagonal
    D2[11, :] = 2*(bx+by+bz)
    D2[11] = D2[11]+Pot
    #Bands below diagonal:
    for ii in range(N.int((bands-1)/2)):
        D2[12+ii, :] = D2[10-ii, ::-1]
    #Define Hamiltonian as sparse array
    mtx = sparse.spdiags([D2[0],               D2[1],               D2[2],
                          D2[3],               D2[4],               D2[5],
                          D2[6],               D2[7],               D2[8],
                          D2[9],               D2[10],              D2[11],
                          N.conjugate(D2[12]), N.conjugate(D2[13]), D2[14],
                          D2[15],              N.conjugate(D2[16]), N.conjugate(D2[17]),
                          N.conjugate(D2[18]), N.conjugate(D2[19]), N.conjugate(D2[20]),
                          N.conjugate(D2[21]), D2[22]],

                         [Nx*Ny,               Nx*Ny-1,             Nx*Ny-Nx+1,
                          Nx*Ny-Nx,            Nx*Ny-Nx-1,          Nx*Ny-2*Nx+1,
                          2*Nx-1,              Nx+1,                Nx,
                          Nx-1,                1,                   0,
                          -1,                  -(Nx-1),             -Nx,
                          -(Nx+1),             -(2*Nx-1),           -(Nx*Ny-2*Nx+1),
                          -(Nx*Ny-Nx-1),       -(Nx*Ny-Nx),         -(Nx*Ny-Nx+1),
                          -(Nx*Ny-1),          -Nx*Ny], NN, NN)
    #Convert Hamiltonian to sparse.csr format
    Ham = sparse.csr_matrix(mtx)
    return Ham


def SplitHam(Ham, inda, indb):
    print 'Splitting FD Hamiltonian in regions w.r.t. the iso surface:'
    print '      -                 -'
    print '     |  Ha           tau |'
    print 'H -> |                   |'
    print '     | tau^\dagger    Hb |'
    print '      -                 -'
    Na, Nb = len(inda), len(indb)
    Tau   = Ham[indb[0:Nb]][:, inda[:Na]]
    Hb    = Ham[indb[0:Nb]][:, indb[0:Nb]]
    print 'Linear system of equations: (E-Hb).phi = tau^{\dagger}.psi, where'
    print 'dim(Hb)  = ('+str(Nb)+'x'+str(Nb)+') lattice points'
    print 'dim(psi) = ('+str(Nb)+'x1) lattice points'
    return Tau, Hb


def LinearSolve(options, WFs, inda, indb, Ef, Tau, Hb, NN, Nx, Ny, Nz, Chans):
    Na, Nb  = len(inda), len(indb)
    pmodes = N.zeros((Chans, NN), N.complex)
    timemodeprop = time.clock()
    print '\nMode propagation starts'
    for mode in range(Chans):
        #Mode at iso surface to propagate (flattened and indexed)
        psiEC = [WFs[mode, ii, jj, kk] for kk in range(Nz) \
                     for jj in range(Ny) for ii in range(Nx)]
        psiA = [psiEC[ii] for ii in inda]
        rhs = Tau.dot(psiA) #rhs of linear system to solve: lhs x = rhs
        EfI = sparse.spdiags([Ef*N.ones(Nb)], [0], Nb, Nb)
        lhs = EfI-Hb
        timesolve = time.clock()
        sol = isolve.gmres(lhs, rhs, tol=options.tolsolve)
        print 'Mode %s found in:'%mode, N.round(time.clock()-timesolve, 2), 's'
        index = N.zeros(NN)
        index[0:Na] = inda
        index[Na:NN]= indb
        invindex = N.argsort(index)
        phi = N.zeros(NN, N.complex)
        phi[0:Na] = psiA
        phi[Na:NN]= sol[0]
        phi = phi[invindex]
        pmodes[mode, :] = phi
        times = N.round(time.clock()-timemodeprop, 2)
        timem = N.round(times/60, 2)
    print 'Finished in', times, 's =', timem, 'min'
    return pmodes


def Current(options, propSubModes, propTipModes, Nx, Ny, Nz, Max, ucSize, ChansL, ChansR, dS, kpoint, pathkpt, ikpoint):
    Ltmp, Rtmp = propSubModes, propTipModes
    a3 = ucSize[2]
    print '\nComputing wave functions and gradients at separation surface ...'
    L  = N.zeros((ChansL, Nx, Ny), N.complex)
    R  = N.zeros((ChansR, Nx, Ny), N.complex)
    dL = N.zeros((ChansL, Nx, Ny), N.complex)
    dR = N.zeros((ChansR, Nx, Ny), N.complex)
    L2, dL2 = N.zeros((ChansL, 2*Nx, 2*Ny), N.complex), N.zeros((ChansL, 2*Nx, 2*Ny), N.complex)
    R2, dR2 = N.zeros((ChansR, 2*Nx, 2*Ny), N.complex), N.zeros((ChansR, 2*Nx, 2*Ny), N.complex)
    for ii in range(ChansL):
        L2[ii, :Nx, :Ny]  = Ltmp[ii, Max*Nx*Ny:(Max+1)*Nx*Ny].reshape(Ny, Nx).transpose()
        L2[ii, Nx:, Ny:]  = L2[ii, :Nx, :Ny]*N.exp(2.0j*N.pi*(kpoint[0]+kpoint[1]))
        L2[ii, Nx:, :Ny]  = L2[ii, :Nx, :Ny]*N.exp(2.0j*N.pi*kpoint[0])
        L2[ii, :Nx, Ny:]  = L2[ii, :Nx, :Ny]*N.exp(2.0j*N.pi*kpoint[1])
        dL2[ii, :Nx, :Ny] = ((Ltmp[ii, (Max+1)*Nx*Ny:(Max+2)*Nx*Ny]\
                        -Ltmp[ii, (Max-1)*Nx*Ny:Max*Nx*Ny])/2/a3).reshape(Ny, Nx).transpose()
        dL2[ii, Nx:, Ny:] = dL2[ii, :Nx, :Ny]*N.exp(2.0j*N.pi*(kpoint[0]+kpoint[1]))
        dL2[ii, Nx:, :Ny] = dL2[ii, :Nx, :Ny]*N.exp(2.0j*N.pi*kpoint[0])
        dL2[ii, :Nx, Ny:] = dL2[ii, :Nx, :Ny]*N.exp(2.0j*N.pi*kpoint[1])
    for ii in range(ChansR):
        R2[ii, :Nx, :Ny]  = Rtmp[ii, (Nz-Max-1)*Nx*Ny:(Nz-Max)*Nx*Ny].reshape(Ny, Nx).transpose()[::-1, ::-1]
        dR2[ii, :Nx, :Ny] = -((Rtmp[ii, (Nz-Max)*Nx*Ny:(Nz-Max+1)*Nx*Ny]\
                             -Rtmp[ii, (Nz-Max-2)*Nx*Ny:(Nz-Max-1)*Nx*Ny])/2/a3).reshape(Ny, Nx).transpose()[::-1, ::-1]
    n=wNC.NCfile(pathkpt+'FD'+N.str(ikpoint)+'.nc')
    n.write(L2.real, 'reL')
    n.write(L2.imag, 'imL')
    n.write(R2.real, 'reR')
    n.write(R2.imag, 'imR')
    n.write(dL2.real, 'redL')
    n.write(dL2.imag, 'imdL')
    n.write(dR2.real, 'redR')
    n.write(dR2.imag, 'imdR')
    n.write(kpoint, 'kpnt')
    n.close()
    L, R, dL, dR = L2, R2, dL2, dR2
    Nx, Ny = 2*Nx, 2*Ny
    scale    = 4*N.pi**2*7.748e-5*1e9*dS**2*options.samplingscale**4
    print '\nSimulating tip scanning by fast Fourier transform ...'
    tot = N.zeros((Nx, Ny))
    currmat = N.zeros((ChansL, ChansR))
    for ii in range(ChansL):
        for jj in range(ChansR):
            tmp1 = N.fft.ifft2(N.fft.fft2(L[ii, ::-1, ::-1])*N.conjugate(N.fft.fft2(dR[jj, :, :])))
            tmp2 = N.fft.ifft2(N.fft.fft2(dL[ii, ::-1, ::-1])*N.conjugate(N.fft.fft2(R[jj, :, :])))
            tot  = tot+N.abs(tmp1-tmp2)**2
            currmat[ii, jj] = N.sum(N.abs(tmp1-tmp2)**2)
    STMcurrent = scale*tot
    STMcurrent = STMcurrent[:Nx/2, :Ny/2]
    n=wNC.NCfile(pathkpt+'FDcurr'+N.str(ikpoint)+'.nc')
    n.write(STMcurrent, 'Curr')
    n.write(kpoint, 'kpnt')
    n.close()
    MaxCurr, MinCurr = N.max(STMcurrent[:Nx]), N.min(STMcurrent)
    print '\nMax conductance:', N.round(MaxCurr, 10), 'nA/V'
    print 'Min conductance:', N.round(MinCurr, 10), 'nA/V'
    print 'Min/Max:', N.round(MinCurr/MaxCurr*100, 3), '%\n'
    return STMcurrent, currmat
