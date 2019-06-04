"""

:mod:`Inelastica.STM`
=====================

Script that calculates STM images using the Bardeen approximation
outlined in PRB 93 115434 (2016) and PRB 96 085415 (2017). The
script is divided into 3 parts:

1. Calculation of the scattering states at the Fermi-energy on the
same real space grid as `TranSIESTA`_ (real-space cutoff). These are
saved in `DestDir/SystemLabel.A[LR][0-99].nc` files and are reused if
found. **NEEDS:** `TranSIESTA`_ calculation.

2. Propagation of the scattering states from a surface (defined by a
constant charge density) out into the vacuum region. After the x-y plane,
where the average potential of the slice is maximum (the separation
plane), is found, the potential is ascribed a constant value at this
average. Saves the propagated wavefunctions at the separation plane in
`DestDir/[kpoint]/FD[kpoint].nc`. **NEEDS:** *TotalPotential.grid.nc* and
*Rho.grid.nc*.

3. Conductance calculation where the tip/substrate wavefunctions are
displaced to simulate the conductance at different tip-positions. The k
averaged STM image and the STM images of individual k points are saved in
`DestDir/STMimage.nc`.

.. currentmodule:: Inelastica.STM

"""
from __future__ import print_function

import numpy as N
import numpy.linalg as LA
import netCDF4 as NC
import sys
import glob
import os
import ast
import time
import Inelastica.misc.valuecheck as VC
import Inelastica.io.log as Log
import Inelastica.io.netcdf as writeNC
import Inelastica.STMFD as STMFD
import Inelastica.NEGF as NEGF
import Inelastica.io.siesta as SIO
import Inelastica.MakeGeom as MG
import Inelastica.math as MM
import Inelastica.physics.constants as PC
import Inelastica.physics.mesh as Kmesh
import Inelastica.misc.multiproc as MP

#Units: Bohr and Rydberg!

# ID for boundaries
zeroSide, wfSide, pbcSide = 1, 2, 3 # ID
DEBUG = False
DEBUGPDE = False


def GetOptions(argv):
    # if text string is specified, convert to list
    if isinstance(argv, str): argv = argv.split()

    import optparse as o

    d = """Script that calculates STM images using the Bardeen approximation outlined in PRB 93 115434 (2016) and PRB 96 085415 (2017). The script is divided into 3 parts:
1) Calculation of the scattering states at the Fermi-energy on the same real space grid as TranSiesta (real-space cutoff). These are saved in DestDir/SystemLabel.A[LR][0-99].nc files and are reused if found. NEEDS: TranSiesta calculation.
2) Propagation of the scattering states from a surface (defined by a constant charge density) out into the vacuum region. After the x-y plane, where the average potential of the slice is maximum (the separation plane), is found, the potential is ascribed a constant value at this average. Saves the propagated wavefunctions at the separation plane in DestDir/[kpoint]/FD[kpoint].nc. NEEDS: TotalPotential.grid.nc and Rho.grid.nc.
3) Conductance calculation where the tip/substrate wavefunctions are displaced to simulate the conductance at different tip-positions. The k averaged STM image and the STM images of individual k points are saved in DestDir/STMimage.nc.
"""

    p = o.OptionParser("usage: %prog [options] DestinationDirectory", description=d)
    p.add_option("-F", "--DeviceFirst", dest='DeviceFirst', default=0, type='int',
                 help="First device atom (SIESTA numbering) [TS.TBT.PDOSFrom]")
    p.add_option("-L", "--DeviceLast", dest='DeviceLast', default=0, type='int',
                 help="Last device atom (SIESTA numbering) [TS.TBT.PDOSTo]")
    p.add_option("-e", "--Energy", dest='energy', default=0.0, type='float',
                 help="Energy where scattering states are evaluated [%default eV]")
    p.add_option("--eta", dest="eta", help="Imaginary part added to all energies (device and leads) [%default eV]",
                 type='float', default=0.000001)
    p.add_option("-l", "--etaLead", dest="etaLead", help="Additional imaginary part added ONLY in the leads (surface GF) [%default eV]",
                 type='float', default=0.0)
    p.add_option("-f", "--fdf", dest='fn', default='./RUN.fdf', type='string',
                 help="Input fdf-file for TranSIESTA calculations [%default]")
    p.add_option("-s", "--iSpin", dest='iSpin', default=0, type='int',
                 help="Spin channel [%default]")
    p.add_option("-p", "--savePOS", dest='savePOS', default=False, action='store_true',
                 help="Save the individual solutions as .pos files")
    p.add_option("--shift", dest='shift', default=False, action='store_true',
                 help="Shift current 1/2 cell in x, y directions")

    # Electrode stuff
    p.add_option("--bulk", dest='UseBulk', default=-1, action='store_true',
                 help="Use bulk in electrodes. The Hamiltonian from the electrode calculation is inserted into the electrode region in the TranSIESTA cell [TS.UseBulkInElectrodes]")
    p.add_option("--nobulk", dest='UseBulk', default=-1, action='store_false',
                 help="Use only self-energies in the electrodes. The full Hamiltonian of the TranSIESTA cell is used in combination with self-energies for the electrodes [TS.UseBulkInElectrodes]")

    # Scale (artificially) the coupling to the electrodes
    p.add_option("--scaleSigL", dest="scaleSigL", help="Scale factor applied to Sigma_L [default=%default]",
                 type='float', default=1.0)
    p.add_option("--scaleSigR", dest="scaleSigR", help="Scale factor applied to Sigma_R [default=%default]",
                 type='float', default=1.0)

    p.add_option("-u", "--useSigNC", dest='signc', default=False, action='store_true',
                 help="Use SigNCfiles [%default]")

    # Use spectral matrices?
    p.add_option("--SpectralCutoff", dest="SpectralCutoff", help="Cutoff value for SpectralMatrix functions (for ordinary matrix representation set cutoff<=0.0) [default=%default]",
                 type='float', default=0.0)
    p.add_option("-n", "--nCPU", dest='nCPU', default=1, type='int',
                 help="Number of processors [%default]")

    # k-grid
    p.add_option("-x", "--Nk1", dest='Nk1', default=1, type='int',
                  help="k-points Nk1 along a1 [%default]")
    p.add_option("-y", "--Nk2", dest='Nk2', default=1, type='int',
                  help="k-points Nk2 along a2 [%default]")

    # FD calculation
    p.add_option("-r", "--rhoiso", dest='rhoiso', default=1e-3, type='float',
                 help="Density at the isosurface from which the localized-basis wave functions are propagated [default=%default Bohr^-3Ry^-1]")
    p.add_option("--ssp", dest='ShiftSeparationPlane', default=0, type='float',
                 help="Manually shift the separation plane (>0 means away from the substrate) [default=%default Ang]")
    p.add_option("--sc", dest='samplingscale', default=2, type='int',
                 help="Sampling scale of wave functions in lateral plane. 1 means same real-space resolution as used in TranSiesta, while 2 means doubling of the lateral lattice constants, etc. [default=%default]")
    p.add_option("-t", "--tolsolve", dest='tolsolve', default=1e-6, type='float',
                 help="Tolerance for the iterative linear solver (scipy.linalg.isolve.gmres) [default=%default]")
    p.add_option("--savelocwfs", dest='savelocwfs', default=False, action='store_true',
                 help="Save localized-basis wave functions.")

    (options, args) = p.parse_args(argv)

    # Get the last positional argument
    options.DestDir = VC.GetPositional(args, "You need to specify a destination directory!")

    # Set module name
    options.module = 'STM'

    options.kpoints = Kmesh.kmesh(options.Nk1, options.Nk2, 1, meshtype=['LIN', 'LIN', 'LIN'], invsymmetry=False)
    return options


def main(options):
    dd = len(glob.glob(options.DestDir))
    if len(glob.glob('TotalPotential.grid.nc')) == 0 and len(glob.glob('Rho.grid.nc')) == 0:
        sys.exit('TotalPotential.nc and Rho.grid.nc not found! Add "SaveTotalPotential true" and "SaveRho true" to RUN.fdf')
    if len(glob.glob('TotalPotential.grid.nc')) == 0:
        sys.exit('TotalPotential.nc not found! Add "SaveTotalPotential true" to RUN.fdf')
    if len(glob.glob('Rho.grid.nc')) == 0:
        sys.exit('Rho.grid.nc not found! Add "SaveRho true" to RUN.fdf')

    Log.CreatePipeOutput(options)
    VC.OptionsCheck(options)
    Log.PrintMainHeader(options)

    ## Step 1: Calculate scattering states from L/R on TranSiesta real space grid.
    if glob.glob(options.DestDir+'/kpoints') != []: # Check previous k-points
        f, oldk = open(options.DestDir+'/kpoints', 'r'), []
        f.readline()
        for ii in f.readlines():
            oldk += [N.array(ii.split(), N.float)]
        oldk = N.array(oldk)

    options.kpoints.mesh2file(options.DestDir+'/kpoints')
    doK = []
    for ikpoint in range(len(options.kpoints.k)):
        ikdir = options.DestDir+'/%i'%(ikpoint)
        if not os.path.isdir(ikdir):
            os.mkdir(ikdir)
        doK += [ikpoint]

    #Check if some k points are already done
    doK = []
    nokpts = len(options.kpoints.k)
    noFDcalcs = len(glob.glob('./'+options.DestDir+'/*/FDcurr*.nc'))
    for ii in range(nokpts):
        if len(glob.glob('./'+options.DestDir+'/'+str(ii)+'/FDcurr'+str(ii)+'.nc')) == 1:
            pass
        else:
            doK += [ii]

    if noFDcalcs > nokpts-1:
        print('\nAll ('+str(nokpts)+') k points are already finished!\n')
    elif noFDcalcs > 0 and noFDcalcs < nokpts:
        print(str(nokpts-len(doK))+' ('+str(N.round(100.*((1.*nokpts-len(doK))/nokpts), 1))+'%) k points already done. Will proceed with the rest.')
        print('You should perhaps remove the loc-basis states in the most recent k folder ')
        print('since some of these may not have been calculated before interuption.\n')
    else:
        if dd == 1:
            print('No STM calculations found in existing directory '+str(options.DestDir)+'/. Starting from scratch.')
            print('(...by possibly using saved localized-basis states)\n')
        else:
            print('STM calculation starts.')
    args = [(options, ik) for ik in doK]
    tmp = MP.runParallel(calcTSWFPar, args, nCPU=options.nCPU)

    print('Calculating k-point averaged STM image')

    def ShiftOrigin(mat, x, y):
        Nx = N.shape(mat)[0]
        Ny = N.shape(mat)[1]
        NewMat = N.zeros((Nx, Ny))
        for ii in range(Nx):
            for jj in range(Ny):
                NewMat[N.mod(ii+x, Nx), N.mod(jj+y, Ny)] = mat[ii, jj]
        return NewMat

    file = NC.Dataset('TotalPotential.grid.nc', 'r')
    steps = N.array(file.variables['cell'][:], N.float)
    theta = N.arccos(N.dot(steps[0], steps[1])/(LA.norm(steps[0])*LA.norm(steps[1])))

    currtmp = NC.Dataset('./'+options.DestDir+'/0/FDcurr0.nc', 'r')
    dimSTMimage = N.shape(currtmp.variables['Curr'][:, :])
    dim1 = dimSTMimage[0]
    dim2 = dimSTMimage[1]
    STMimage = N.zeros((dim1, dim2))
    tmpSTM = N.zeros((dim1*nokpts, dim2))
    for ii in range(nokpts):
        file = NC.Dataset('./'+options.DestDir+'/'+str(ii)+'/FDcurr'+str(ii)+'.nc', 'r')
        ikSTMimage = file.variables['Curr'][:, :]/(nokpts)
        STMimage += ShiftOrigin(ikSTMimage, dim1 // 2, dim2 // 2)
        tmpSTM[ii*dim1:(ii+1)*dim1, :] = ShiftOrigin(ikSTMimage, dim1 // 2, dim2 // 2)
    STMimagekpt = N.zeros((dim1*options.Nk1, dim2*options.Nk2))
    for ii in range(options.Nk1):
        for jj in range(options.Nk2):
            N1 = ii+1
            N2 = options.Nk2-jj
            kk = (ii+1)*options.Nk2-jj-1
            STMimagekpt[(N1-1)*dim1:N1*dim1, (N2-1)*dim2:N2*dim2] = tmpSTM[kk*dim1:(kk+1)*dim1, ::-1]

    tmp = open(options.systemlabel+'.XV').readlines()[4+options.DeviceFirst-1:4+options.DeviceLast]
    xyz = N.zeros((len(tmp), 4))
    for ii in range(len(tmp)):
        for jj in range(1, 5):
            tmp2 = tmp[ii].split()
            xyz[ii, jj-1] = ast.literal_eval(tmp2[jj])
            if jj > 1:
                xyz[ii, jj-1] = xyz[ii, jj-1]*PC.Bohr2Ang

    drAtoms = N.zeros(len(xyz))
    for atomIdx in range(len(xyz)-1):
        drAtoms[atomIdx] = N.round((xyz[atomIdx+1][3]-xyz[atomIdx][3]), 3)
    TipHeight = N.max(drAtoms)

    n = writeNC.NCfile('./'+options.DestDir+'/STMimage.nc')
    n.write(STMimage, 'STMimage')
    n.write(theta, 'theta')
    n.write(options.kpoints.k, 'kpoints')
    n.write(options.Nk1, 'Nk1')
    n.write(options.Nk2, 'Nk2')
    n.write(STMimagekpt, 'STMkpoint')
    n.write(xyz, 'Geometry')
    n.write(TipHeight, 'TipHeight')
    n.close()

    Log.PrintMainFooter(options)

########################################################


def calcTSWF(options, ikpoint):
    kpoint = options.kpoints.k[ikpoint]

    def calcandwrite(A, txt, kpoint, ikpoint):
        if isinstance(A, MM.SpectralMatrix):
            A = A.full()
        A = A*PC.Rydberg2eV # Change to 1/Ryd
        ev, U = LA.eigh(A)
        Utilde = N.empty(U.shape, U.dtype)
        for jj, val in enumerate(ev): # Problems with negative numbers
            if val < 0:
                val = 0
            Utilde[:, jj] = N.sqrt(val/(2*N.pi))*U[:, jj]
        indx2 = N.where(abs(abs(ev) > 1e-4))[0] # Pick non-zero states

        ev = ev[indx2]
        Utilde = Utilde[:, indx2]
        indx = ev.real.argsort()[::-1]
        fn = options.DestDir+'/%i/'%(ikpoint)+options.systemlabel+'.%s'%(txt)

        path = './'+options.DestDir+'/'
        #FermiEnergy = SIO.HS(options.systemlabel+'.TSHS').ef

        def noWfs(side):
            tmp = len(glob.glob(path+str(ikpoint)+'/'+options.systemlabel+'.A'+side+'*'))
            return tmp
        if noWfs('L') == 0 or noWfs('R') == 0 and len(glob.glob(path+str(ikpoint)+'/FD*')) == 0:
            print('Calculating localized-basis states from spectral function %s ...'%(txt))
            tlb = time.clock()
            calcWF2(options, geom, options.DeviceAtoms, basis, Utilde[:, indx], [N1, N2, N3, minN3, maxN3], Fold=True, k=kpoint, fn=fn)
            times = N.round(time.clock()-tlb, 2)
            timem = N.round(times/60, 2)
            print('Finished in '+str(times)+' s = '+str(timem)+' min')

        if noWfs('L') > 0 and noWfs('R') > 0 and len(glob.glob(path+str(ikpoint)+'/FD*')) == 0 and str('%s'%(txt)) == str('AR'):
            print('\nLocalized-basis states are calculated in k point '+str(ikpoint)+'/.')
            print('------------------------------------------------------')
            print('Finite-difference calculation of vacuum states starts!')
            timeFD = time.clock()
            STMFD.main(options, kpoint, ikpoint)
            print('FD calculation in k-point folder '+str(ikpoint)+'/ done in '+str(N.round((time.clock()-timeFD)/60, 2))+' min.')
            print('------------------------------------------------------')
            if options.savelocwfs == False:
                os.system('rm -f '+path+str(ikpoint)+'/'+options.systemlabel+'*')

    #Read geometry
    XV = '%s/%s.XV'%(options.head, options.systemlabel)
    geom = MG.Geom(XV, BufferAtoms=options.buffer)

    #Set up device Greens function
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL, options.NA1L, options.NA2L, options.voltage/2.)
    elecL.scaling = options.scaleSigL
    elecL.semiinf = options.semiinfL
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR, options.NA1R, options.NA2R, -options.voltage/2.)
    elecR.scaling = options.scaleSigR
    elecR.semiinf = options.semiinfR
    DevGF = NEGF.GF(options.TSHS, elecL, elecR, Bulk=options.UseBulk,
                    DeviceAtoms=options.DeviceAtoms,
                    BufferAtoms=options.buffer)

    DevGF.calcGF(options.energy+options.eta*1.0j, kpoint[0:2], ispin=options.iSpin,
                 etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
    NEGF.SavedSig.close() #Make sure saved Sigma is written to file
    #Transmission
    print('Transmission Ttot(%.4feV) = %.16f'%(options.energy, N.trace(DevGF.TT).real))

    #Build basis
    options.nspin = DevGF.HS.nspin
    L = options.bufferL
    #Pad lasto with zeroes to enable basis generation...
    lasto = N.zeros((DevGF.HS.nua+L+1,), N.int)
    lasto[L:] = DevGF.HS.lasto
    basis = SIO.BuildBasis(options.fn,
                           options.DeviceAtoms[0]+L,
                           options.DeviceAtoms[1]+L, lasto)
    basis.ii -= L

    file = NC.Dataset('TotalPotential.grid.nc', 'r')
    N1, N2, N3 = len(file.dimensions['n1']), len(file.dimensions['n2']), len(file.dimensions['n3'])
    cell = file.variables['cell'][:]
    file.close()

    #Find device region in a3 axis
    U = LA.inv(N.array([cell[0]/N1, cell[1]/N2, cell[2]/N3]).transpose())
    gridindx = N.dot(geom.xyz[options.DeviceAtoms[0]-1:options.DeviceAtoms[1]]/PC.Bohr2Ang, U)
    minN3, maxN3 = N.floor(N.min(gridindx[:, 2])).astype(N.int), N.ceil(N.max(gridindx[:, 2])).astype(N.int)
    if not N.allclose(geom.pbc, cell*PC.Bohr2Ang):
        print('Error: TotalPotential.grid.nc has different cell compared to geometry')
        sys.exit(1)

    print('\nk-point folder '+str(ikpoint)+'/')
    print('Checking/calculating localized-basis states ...')
    calcandwrite(DevGF.AL, 'AL', kpoint, ikpoint)
    calcandwrite(DevGF.AR, 'AR', kpoint, ikpoint)


def calcTSWFPar(resQue, ii, options, ik):
    calcTSWF(options, ik)
    resQue.put((ii, True))


def calcWF2(options, geom, DeviceAtoms, basis, Y, NN, Fold=True, k=[0, 0, 0], a=None, fn=''):
    """
    Calculate wavefunction on real space mesh with optional folding
    of the periodic boundary conditions. If Fold==True the vectors 'a'
    will be choosen to be the periodic boundary condition vectors and
    the real space wavefunction folded into the cell using the k-vector.
    If Folded==False, you can choose any 'a' but folding by the PBC will
    not be done.
    Note: Folding assumes coupling only between N.N. cells

    INPUT:
      geom        : MakeGeom structure for full geometry
      DeviceAtoms : [first, last] numbering from 1 to N
      basis       : Basis struct for ONLY the device region
      Y           : Wavefunction for the device region
      NN          : [N1,N2,N3,minN3,maxN3] number of points along [a1, a2, a3]
                    only calculate from minN3 to maxN3
      Fold        : fold periodic boundary conditions
      k           : k-vector to use for folding [-0.5,0.5]
      a           : if Fold==True, uses geom.pbc, if false: [a1, a2, a3]
                    along these directions, i.e., da1=a1/N1 ...

    RETURNS:
    YY : complex wavefunction as N1, N2, N3 matrix
    """
    xyz = N.array(geom.xyz[DeviceAtoms[0]-1:DeviceAtoms[1]])
    N1, N2, N3, minN3, maxN3 = NN[0], NN[1], NN[2], NN[3], NN[4]
    NN3 = maxN3-minN3+1
    if Fold:
        da1, da2, da3 = geom.pbc[0]/N1, geom.pbc[1]/N2, geom.pbc[2]/N3
    else:
        da1, da2, da3 = a[0]/N1, a[1]/N2, a[2]/N3
    try:
        NNY = Y.shape[1]
    except:
        NNY = 1
    # Def cube
    YY = [N.zeros((N1, N2, NN3), N.complex) for ii in range(NNY)]
    # Help indices
    i1 = MM.outerAdd(N.arange(N1), N.zeros(N2), N.zeros(NN3))
    i2 = MM.outerAdd(N.zeros(N1), N.arange(N2), N.zeros(NN3))
    i3 = MM.outerAdd(N.zeros(N1), N.zeros(N2), N.arange(NN3)+minN3)

    # Positions x,y,z for points in matrix form
    rx = i1*da1[0]+i2*da2[0]+i3*da3[0]
    ry = i1*da1[1]+i2*da2[1]+i3*da3[1]
    rz = i1*da1[2]+i2*da2[2]+i3*da3[2]

    if Fold:
        pbc = N.array([da1*N1, da2*N2, da3*N3])
        orig = da3*minN3
        b = N.transpose(LA.inv(pbc))
        olist = [orig, orig+pbc[0]+pbc[1]+pbc[2]]
        pairs = N.array([[b[1], b[2]], [b[0], b[2]], [b[0], b[1]]])
        pbcpairs = N.array([[pbc[1], pbc[2]], [pbc[0], pbc[2]], [pbc[0], pbc[1]]])

    for iiatom in range(DeviceAtoms[1]-DeviceAtoms[0]+1):
        if iiatom > 0:
            #Wavefunction percent done...
            SIO.printDone(iiatom, DeviceAtoms[1]-DeviceAtoms[0]+1, 'Wave function')

        if not Fold:
            shiftvec = [0, 0, 0]
        else: # include shift vectors if sphere cuts planes of PBC
            basisindx = N.where(basis.ii == iiatom+DeviceAtoms[0])[0]
            basisradius = N.max(basis.coff[basisindx])
            minshift, maxshift = [0, 0, 0], [0, 0, 0]
            for iithiso, thiso in enumerate(olist):
                c = xyz[iiatom, :]-thiso
                for ip in range(3):
                    n = N.cross(pbcpairs[ip, 0], pbcpairs[ip, 1])
                    n = n / LA.norm(n)
                    dist = N.abs(N.dot(n, N.dot(N.dot(pairs[ip], c), pbcpairs[ip])-c))
                    if dist < basisradius:
                        if iithiso == 0:
                            minshift[ip] = -1
                        else:
                            maxshift[ip] = 1
            shiftvec = []
            for ix in range(minshift[0], maxshift[0]+1):
                for iy in range(minshift[1], maxshift[1]+1):
                    for iz in range(minshift[2], maxshift[2]+1):
                        shiftvec += [[ix, iy, iz]]

        #print(shiftvec)
        for shift in shiftvec:
            # Atom position shifted
            vec = -N.dot(N.array(shift), geom.pbc)
            phase = N.exp(-1.0j*(2*N.pi*N.dot(N.array(k), N.array(shift))))

            # Difference and polar, cylinder distances
            dx, dy, dz = rx-xyz[iiatom, 0]-vec[0], ry-xyz[iiatom, 1]-vec[1], rz-xyz[iiatom, 2]-vec[2]
            dr2 = dx*dx+dy*dy+dz*dz
            drho2 = dx*dx+dy*dy

            basisindx = N.where(basis.ii == iiatom+DeviceAtoms[0])[0]
            old_coff, old_delta = 0.0, 0.0
            for basisorb in basisindx:
                if not (basis.coff[basisorb] == old_coff and basis.delta[basisorb] == old_delta):
                    old_coff, old_delta = basis.coff[basisorb], basis.delta[basisorb]
                    indx = N.where(dr2 < basis.coff[basisorb]**2) # Find points close to atom

                    idr, idrho = N.sqrt(dr2[indx]), N.sqrt(drho2[indx])
                    iri = (idr/basis.delta[basisorb]).astype(N.int)
                    idx, idy, idz = dx[indx], dy[indx], dz[indx]

                    costh = idz/idr
                    cosfi, sinfi = idx/idrho, idy/idrho

                    # Fix divide by zeros
                    indxRzero, indxRhozero = N.where(idr == 0.0), N.where(idrho == 0.0)
                    costh[indxRzero] = 1.0
                    cosfi[indxRhozero], sinfi[indxRhozero] = 1.0, 0.0

                    # Numpy has changed the choose function to crap!
                    RR = N.take(basis.orb[basisorb], iri)
                # Calculate spherical harmonics
                if len(idr) > 0:
                    l = basis.L[basisorb]
                    m = basis.M[basisorb]
                    if l == 3:
                        print('f-shell : l=%i, m=%i (NOT TESTED!!)'%(l, m))
                    thisSphHar = MM.sphericalHarmonics(l, m, costh, sinfi, cosfi)
                    for iy in range(NNY):
                        YY[iy][indx] = YY[iy][indx]+RR*thisSphHar*Y[basisorb, iy]*phase
    N1, N2, N3, minN3, maxN3 = NN[0], NN[1], NN[2], NN[3], NN[4]
    for ii in range(NNY):
        writenetcdf2(geom, fn+'%i.nc'%ii, YY[ii], N1, N2, N3, minN3, maxN3, geom.pbc, options.DeviceAtoms)
    return YY

########################################################


def writenetcdf2(geom, fn, YY, nx, ny, nz, minnz, maxnz, pbc, DeviceAtoms):
    """
    MP: Write STM files to netcdf
    """
    file = NC.Dataset(fn, 'w', 'Created '+time.ctime(time.time()))
    file.createDimension('nx', nx)
    file.createDimension('ny', ny)
    file.createDimension('nz', maxnz-minnz+1)
    file.createDimension('natoms', DeviceAtoms[1]-DeviceAtoms[0]+1)
    file.createDimension('naxes', 3)
    file.createDimension('number', 1)
    # Fields
    varRe = file.createVariable('Re-Psi', 'd', ('nx', 'ny', 'nz'), zlib=True)
    varRe[:] = YY.real*(PC.Bohr2Ang**(3.0/2.0))
    varRe.units = '1/[Ryd^(1/2) Bohr^(3/2)]'
    varIm = file.createVariable('Im-Psi', 'd', ('nx', 'ny', 'nz'), zlib=True)
    varIm[:] = YY.imag*(PC.Bohr2Ang**(3.0/2.0))
    varIm.units = '1/[Ryd^(1/2) Bohr^(3/2)]'
    varcell = file.createVariable('cell', 'd', ('naxes', 'naxes'))
    varcell[:] = pbc/PC.Bohr2Ang
    varcell.units = 'Bohr'
    vardsteps = file.createVariable('steps', 'd', ('naxes', 'naxes'))
    steps = N.array([pbc[0]/nx, pbc[1]/ny, pbc[2]/nz])/PC.Bohr2Ang
    vardsteps[:] = steps
    vardsteps.units = 'Bohr'
    varorig = file.createVariable('origin', 'd', ('naxes',))
    varorig[:] = pbc[2]/nz*minnz/PC.Bohr2Ang
    varorig.units = 'Bohr'
    vargeom = file.createVariable('xyz', 'd', ('natoms', 'naxes'))
    vargeom[:] = geom.xyz[list(range(DeviceAtoms[0]-1, DeviceAtoms[1]))]/PC.Bohr2Ang
    vargeom.units = 'Bohr'
    varanr = file.createVariable('anr', 'i', ('natoms',))
    varanr[:] = N.array(geom.anr[DeviceAtoms[0]-1:DeviceAtoms[1]], N.int32)
    file.close()

##################### Start main routine #####################
if __name__ == '__main__':
    from datetime import datetime
    start = datetime.now()
    options = GetOptions(sys.argv[1:])
    print(dir(options))
    #profile.run('main(options)')
    main(options)

    dT = datetime.now()-start
    Log.PrintScriptSummary(sys.argv, dT)
