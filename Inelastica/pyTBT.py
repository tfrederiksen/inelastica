"""
################################################################

 python TBTrans 
 Magnus Paulsson magnus.paulsson@hik.se

 Requires: numpy (compile it linked with mkl, acml or atlas!)
           For speed compile the fortran subroutines in F90 
           (cd F90;source compile.bat)

  UNITS! Always eV and Angstrom!
         k-values always given in range [0,1.0] (or [-0.5,0.5])
         They are not in reciprocal space. Instead they corresponds
         to the mathematical orthogonal space that is fourier 
         transformed.

################################################################
"""

import Inelastica.io.siesta as SIO
import Inelastica.MiscMath as MM
import Inelastica.NEGF
import numpy as N
import Inelastica.physics.mesh as Kmesh
import Inelastica.ValueCheck as VC
import Inelastica.CommonFunctions as CF
import Inelastica.io.xmgrace as XMGR

# For doing loops with pyTBT we encourage the usage of this function
# By creating the parser locally we can actually pass down these informations easily.
# DIRECTLY in python


def GetOptions(argv, **kwargs):
    # if text string is specified, convert to list
    if isinstance(argv, VC.string_types):
        argv = argv.split()

    import argparse

    p = argparse.ArgumentParser(description='pyTBT is the Python version of TBtrans originally developed by Mads Brandbyge.')
    p.add_argument('DestDir', help='Destination directory')
    # keywords with defaults from fdf-file:
    p.add_argument('-f', '--fdf', dest='fn', default='./RUN.fdf', type=str,
                 help='Input fdf-file for TranSIESTA calculation [%(default)s]')
    p.add_argument('-F', '--DeviceFirst', dest='DeviceFirst', default=0, type=int,
                 help='First device atom (SIESTA numbering) [TS.TBT.PDOSFrom]')
    p.add_argument('-L', '--DeviceLast', dest='DeviceLast', default=0, type=int,
                 help='Last device atom (SIESTA numbering) [TS.TBT.PDOSTo]')
    p.add_argument('-N', '--NPoints', dest='NPoints', default=0, type=int,
                 help='Energy points [TS.TBT.NPoints]')
    p.add_argument('--Emin', dest='Emin', default=1e10, type=float,
                 help='First energy point [TS.TBT.Emin]')
    p.add_argument('--Emax', dest='Emax', default=1e10, type=float,
                 help='Last energy point [TS.TBT.Emax]')

    # k-points related
    p.add_argument('-x', '--Nk1', dest='Nk1', default=1, type=int,
                 help='k-points Nk1 along a1 [%(default)s]')
    p.add_argument('-y', '--Nk2', dest='Nk2', default=1, type=int,
                 help='k-points Nk2 along a2 [%(default)s]')
    p.add_argument('-a', '--Gk1', dest='Gk1', default=0, type=int,
                 help='Gaussian quadrature k-point sampling for a1 direction (2*GK1+1 points) [%(default)s]')
    p.add_argument('-b', '--Gk2', dest='Gk2', default=0, type=int,
                 help='Gaussian quadrature k-point sampling for a2 direction (2*GK2+1 points) [%(default)s]')
    p.add_argument('-s', '--skipsym', dest='skipsymmetry', default=False, action='store_true',
                 help='Skip inversion (time-reversal) symmetry (i.e., k=-k) that reduces the number of k-point evaluations')
    p.add_argument('-j', '--singlejunction', dest='singlejunction', default=False, action='store_true',
                 help='k-point sample only electrode self-energies')

    # Imaginary part to Greens functions
    p.add_argument('-e', '--eta', dest='eta', type=float, default=1e-6,
                 help='Imaginary part added to all energies (device and leads) [%(default)s eV]')
    p.add_argument('-l', '--etaLead', dest='etaLead', type=float, default=0.0,
                 help='Additional imaginary part added ONLY in the leads (surface GF) [%(default)s eV]')

    # Other options
    p.add_argument('-d', '--skipDOS', dest='dos', default=True, action='store_false',
                 help='Skip calculation of PDOS')
    p.add_argument('--useSigNC', dest='signc', default=False, action='store_true',
                 help='Use SigNCfiles')
    p.add_argument('--NumChan', dest='numchan', type=int, default=10,
                 help='Number of eigenchannels [%(default)s]')

    # Electrode stuff
    p.add_argument('--bulk', dest='UseBulk', default=-1, action='store_true',
                 help='Use bulk in electrodes. The Hamiltonian from the electrode calculation is inserted into the electrode region in the TranSIESTA cell [TS.UseBulkInElectrodes]')
    p.add_argument('--nobulk', dest='UseBulk', default=-1, action='store_false',
                 help='Use only self-energies in the electrodes. The full Hamiltonian of the TranSIESTA cell is used in combination with self-energies for the electrodes [TS.UseBulkInElectrodes]')

    # Scale (artificially) coupling to electrodes
    p.add_argument('--scaleSigL', dest='scaleSigL', type=float, default=1.0,
                 help='Scale factor applied to Sigma_L [default=%(default)s]')
    p.add_argument('--scaleSigR', dest='scaleSigR', type=float, default=1.0,
                 help='Scale factor applied to Sigma_R [default=%(default)s]')

    # Use spectral matrices?
    p.add_argument('--SpectralCutoff', dest='SpectralCutoff', type=float, default=0.0,
                 help='Cutoff value for SpectralMatrix functions (for ordinary matrix representation set cutoff<=0.0) [default=%(default)s]')

    # Parse the options
    options = p.parse_args(argv)

    # With this one can overwrite the logging information
    if "log" in kwargs:
        options.Logfile = kwargs['log']
    else:
        options.Logfile = 'pyTBT.log'

    return options


def main(options):
    CF.CreatePipeOutput(options.DestDir+'/'+options.Logfile)
    VC.OptionsCheck(options, 'pyTBT')
    CF.PrintMainHeader('pyTBT', options)

    # K-points
    if options.Gk1>1:
        Nk1, t1 = options.Gk1, 'GK'
    else:
        Nk1, t1 = options.Nk1, 'LIN'
    if options.Gk2>1:
        Nk2, t2 = options.Gk2, 'GK'
    else:
        Nk2, t2 = options.Nk2, 'LIN'
    # Generate full k-mesh:
    mesh = Kmesh.kmesh(Nk1, Nk2, Nk3=1, meshtype=[t1, t2, 'LIN'], invsymmetry=not options.skipsymmetry)
    mesh.mesh2file('%s/%s.%ix%i.mesh'%(options.DestDir, options.systemlabel, mesh.Nk[0], mesh.Nk[1]))
    # Setup self-energies and device GF
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL, options.NA1L, options.NA2L, options.voltage/2.)
    elecL.scaling = options.scaleSigL
    elecL.semiinf = options.semiinfL
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR, options.NA1R, options.NA2R, -options.voltage/2.)
    elecR.scaling = options.scaleSigR
    elecR.semiinf = options.semiinfR
    DevGF = NEGF.GF(options.TSHS, elecL, elecR, Bulk=options.UseBulk,
                    DeviceAtoms=options.DeviceAtoms,
                    BufferAtoms=options.buffer)
    nspin = DevGF.HS.nspin

    # k-sample only self-energies?
    if options.singlejunction:
        elecL.mesh = mesh
        mesh = Kmesh.kmesh(3, 3, 1)

    if options.dos:
        DOSL = N.zeros((nspin, len(options.Elist), DevGF.nuo), N.float)
        DOSR = N.zeros((nspin, len(options.Elist), DevGF.nuo), N.float)

        # MPSH projections?
        MPSHL = N.zeros((nspin, len(options.Elist), DevGF.nuo), N.float)
        MPSHR = N.zeros((nspin, len(options.Elist), DevGF.nuo), N.float)
        # evaluate eigenstates at Gamma
        import scipy.linalg as SLA
        DevGF.setkpoint(N.zeros(2))
        ev0, es0 = SLA.eigh(DevGF.H, DevGF.S)
        print 'MPSH eigenvalues:', ev0
        #print 'MPSH eigenvector normalizations:',N.diag(MM.mm(MM.dagger(es0),DevGF.S,es0)).real # right

    # Loop over spin
    for iSpin in range(nspin):
        # initialize transmission and shot noise arrays
        Tkpt=N.zeros((len(options.Elist), mesh.NNk, options.numchan+1), N.float)
        SNkpt=N.zeros((len(options.Elist), mesh.NNk, options.numchan+1), N.float)
        # prepare output files
        outFile = options.DestDir+'/%s.%ix%i'%(options.systemlabel, mesh.Nk[0], mesh.Nk[1])
        if nspin<2: thisspinlabel = outFile
        else: thisspinlabel = outFile+['.UP', '.DOWN'][iSpin]
        fo=open(thisspinlabel+'.AVTRANS', 'write')
        fo.write('# Nk1(%s)=%i Nk2(%s)=%i eta=%.2e etaLead=%.2e\n'%(mesh.type[0], mesh.Nk[0], mesh.type[1], mesh.Nk[1], options.eta, options.etaLead))
        fo.write('# E   Ttot(E)   Ti(E)(i=1-%i)   RelErrorTtot(E)\n'%options.numchan)
        foSN=open(thisspinlabel+'.AVNOISE', 'write')
        foSN.write('# Nk1(%s)=%i Nk2(%s)=%i eta=%.2e etaLead=%.2e\n'%(mesh.type[0], mesh.Nk[0], mesh.type[1], mesh.Nk[1], options.eta, options.etaLead))
        foSN.write('# E   SNtot(E)   SNi(E)(i=1-%i)\n'%options.numchan)
        foFF=open(thisspinlabel+'.FANO', 'write')
        foFF.write('# Nk1(%s)=%i Nk2(%s)=%i eta=%.2e etaLead=%.2e\n'%(mesh.type[0], mesh.Nk[0], mesh.type[1], mesh.Nk[1], options.eta, options.etaLead))
        foFF.write('# E   Fano factor \n')
        # Loop over energy
        for ie, ee in enumerate(options.Elist):
            Tavg = N.zeros((options.numchan+1, len(mesh.w)), N.float)
            SNavg = N.zeros((options.numchan+1, len(mesh.w)), N.float)
            AavL = N.zeros((DevGF.nuo, DevGF.nuo), N.complex)
            AavR = N.zeros((DevGF.nuo, DevGF.nuo), N.complex)
            # Loops over k-points
            for ik in range(mesh.NNk):
                DevGF.calcGF(ee+options.eta*1.0j, mesh.k[ik, :2], ispin=iSpin,
                             etaLead=options.etaLead, useSigNCfiles=options.signc, SpectralCutoff=options.SpectralCutoff)
                # Transmission and shot noise
                T, SN = DevGF.calcTEIG(options.numchan)
                for iw in range(len(mesh.w)):
                    Tavg[:, iw] += T*mesh.w[iw, ik]
                    SNavg[:, iw] += SN*mesh.w[iw, ik]
                Tkpt[ie, ik] = T
                SNkpt[ie, ik] = SN
                # DOS calculation:
                if options.dos:
                    if options.SpectralCutoff>0.0:
                        AavL += mesh.w[0, ik]*MM.mm(DevGF.AL.L, DevGF.AL.R, DevGF.S)
                        AavR += mesh.w[0, ik]*MM.mm(DevGF.AR.L, DevGF.AR.R, DevGF.S)
                    else:
                        AavL += mesh.w[0, ik]*MM.mm(DevGF.AL, DevGF.S)
                        AavR += mesh.w[0, ik]*MM.mm(DevGF.AR, DevGF.S)
            # Print calculated quantities
            err = (N.abs(Tavg[0, 0]-Tavg[0, 1])+N.abs(Tavg[0, 0]-Tavg[0, 2]))/2
            relerr = err/Tavg[0, 0]
            print 'ispin= %i, e= %.4f, Tavg= %.8f, RelErr= %.1e'%(iSpin, ee, Tavg[0, 0], relerr)
            transline = '\n%.10f '%ee
            noiseline = '\n%.10f '%ee
            for ichan in range(options.numchan+1):
                if ichan==0:
                    transline += '%.8e '%Tavg[ichan, 0]
                    noiseline += '%.8e '%SNavg[ichan, 0]
                else:
                    transline += '%.4e '%Tavg[ichan, 0]
                    noiseline += '%.4e '%SNavg[ichan, 0]
            transline += '%.2e '%relerr
            fo.write(transline)
            foSN.write(noiseline)
            foFF.write('\n%.10f %.8e'%(ee, SNavg[0, 0]/Tavg[0, 0]))
            # Partial density of states:
            if options.dos:
                DOSL[iSpin, ie, :] += N.diag(AavL).real/(2*N.pi)
                DOSR[iSpin, ie, :] += N.diag(AavR).real/(2*N.pi)
                MPSHL[iSpin, ie, :] += N.diag(MM.mm(MM.dagger(es0), AavL, es0)).real/(2*N.pi)
                MPSHR[iSpin, ie, :] += N.diag(MM.mm(MM.dagger(es0), AavR, es0)).real/(2*N.pi)
                print 'ispin= %i, e= %.4f, DOSL= %.4f, DOSR= %.4f'%(iSpin, ee, N.sum(DOSL[iSpin, ie, :]), N.sum(DOSR[iSpin, ie, :]))
        fo.write('\n')
        fo.close()
        foSN.write('\n')
        foSN.close()
        foFF.write('\n')
        foFF.close()

        # Write k-point-resolved transmission
        fo=open(thisspinlabel+'.TRANS', 'write')
        for ik in range(mesh.NNk):
            w = mesh.w[:, ik]
            fo.write('\n\n# k = %f, %f    w = %f %f %f %f'%(mesh.k[ik, 0], mesh.k[ik, 1], w[0], w[1], w[2], w[3]))
            for ie, ee in enumerate(options.Elist):
                transline = '\n%.10f '%ee
                for ichan in range(options.numchan+1):
                    if ichan==0:
                        transline += '%.8e '%Tkpt[ie, ik, ichan]
                    else:
                        transline += '%.4e '%Tkpt[ie, ik, ichan]
                fo.write(transline)
        fo.write('\n')
        fo.close()

        # Write k-point-resolved shot noise
        fo=open(thisspinlabel+'.NOISE', 'write')
        for ik in range(mesh.NNk):
            w = mesh.w[:, ik]
            fo.write('\n\n# k = %f, %f    w = %f %f %f %f'%(mesh.k[ik, 0], mesh.k[ik, 1], w[0], w[1], w[2], w[3]))
            for ie, ee in enumerate(options.Elist):
                noiseline = '\n%.10f '%ee
                for ichan in range(options.numchan+1):
                    if ichan==0:
                        noiseline += '%.8e '%SNkpt[ie, ik, ichan]
                    else:
                        noiseline += '%.4e '%SNkpt[ie, ik, ichan]
                fo.write(noiseline)
        fo.write('\n')
        fo.close()

    # End loop over spin
    NEGF.SavedSig.close() # Make sure saved Sigma is written to file

    if options.dos:
        # Read basis
        L = options.bufferL
        # Pad lasto with zeroes to enable basis generation...
        lasto = N.zeros((DevGF.HS.nua+L+1,), N.int)
        lasto[L:] = DevGF.HS.lasto
        basis = SIO.BuildBasis(options.fn, 1+L, DevGF.HS.nua+L,
                               lasto)
        basis.ii -= L
        WritePDOS(outFile+'.PDOS.gz', options, DevGF, DOSL+DOSR, basis)
        WritePDOS(outFile+'.PDOSL.gz', options, DevGF, DOSL, basis)
        WritePDOS(outFile+'.PDOSR.gz', options, DevGF, DOSR, basis)

        WriteMPSH(outFile+'.MPSH.gz', options, DevGF, MPSHL+MPSHR, ev0)
        WriteMPSH(outFile+'.MPSHL.gz', options, DevGF, MPSHL, ev0)
        WriteMPSH(outFile+'.MPSHR.gz', options, DevGF, MPSHR, ev0)

    CF.PrintMainFooter('pyTBT')


def WritePDOS(fn, options, DevGF, DOS, basis):
    """
    PDOS from the surface Green's function from electrode calculations

    NOTE! The DOS is a sum over the atoms of the unitcell.
    NOTE! The outfile contains the DOS divided into s,p,d,f shells.
    This decomposition is not perfect since polarized basis orbitals
    will end up in L+1 part, i.e., 6s polarized orbital = 6p
    """

    import xml.dom.minidom as xml
    import gzip

    # First, last orbital in full space and pyTBT folded space.
    devOrbSt = DevGF.HS.lasto[options.DeviceAtoms[0]-1]
    pyTBTdevOrbSt = devOrbSt-DevGF.HS.lasto[options.DeviceAtoms[0]-1]
    devOrbEnd = DevGF.HS.lasto[options.DeviceAtoms[1]]-1
    pyTBTdevOrbEnd = devOrbEnd-DevGF.HS.lasto[options.DeviceAtoms[0]-1]

    doc = xml.Document()
    pdos = doc.createElement('pdos')
    doc.appendChild(pdos)
    xmladd(doc, pdos, 'nspin', '%i'%DevGF.HS.nspin)
    xmladd(doc, pdos, 'norbitals', '%i'%(DevGF.nuo))
    xmladd(doc, pdos, 'energy_values', myprint(options.Elist+DevGF.HS.ef))
    xmladd(doc, pdos, 'E_Fermi', '%.8f'%DevGF.HS.ef)
    for ii in range(DevGF.nuo):
        orb = doc.createElement('orbital')
        pdos.appendChild(orb)
        io = devOrbSt+ii
        orb.setAttribute('index', '%i'%(io+1))
        orb.setAttribute('atom_index', '%i'%basis.ii[io])
        orb.setAttribute('species', basis.label[io])
        orb.setAttribute('position', '%f %f %f'%(basis.xyz[io, 0], basis.xyz[io, 1], basis.xyz[io, 2]))
        orb.setAttribute('n', '%i'%basis.N[io])
        orb.setAttribute('l', '%i'%basis.L[io])
        orb.setAttribute('m', '%i'%basis.M[io])
        xmladd(doc, orb, 'data', myprint(DOS[:, :, ii]))
    doc.writexml(gzip.GzipFile(fn, 'w'))

    # Make plot
    atoms = list(set(basis.label))
    lVals  = list(set(basis.L))
    plots = [[atoms, lVals, 'Tot']]
    plots += [[atoms, [lVal], 'Tot L=%i'%lVal] for lVal in lVals]
    plots += [[[atom], lVals, atom+' Tot'] for atom in atoms]
    plots += [[[atom], [lVal], atom+' L=%i'%lVal] for lVal in lVals for atom in atoms]

    g = XMGR.Graph()
    for atom, lVal, name in plots:
        nspin, ee, PDOS = SIO.ExtractPDOS(fn, None, FermiRef=False, llist=lVal,
                                          species=atom, Normalize=True)
        for iS in range(nspin):
            g.AddDatasets(XMGR.XYset(ee-DevGF.HS.ef, (-1)**iS*PDOS[iS], legend=name, Lwidth=2))

    # Set axes and write XMGR plot to file
    g.SetXaxis(label='E-E\sF\N (eV)', autoscale=True)
    g.SetYaxis(label='DOS (1/eV/atom)', autoscale=True)
    g.SetTitle(fn, size=1.3)
    g.ShowLegend()
    p = XMGR.Plot(fn+'.xmgr', g)
    p.WriteFile()


def WriteMPSH(fn, options, DevGF, DOS, ev0):
    """
    Projected density of states onto MPSH eigenstates (at Gamma)
    """

    import xml.dom.minidom as xml
    import gzip

    doc = xml.Document()
    mpsh = doc.createElement('mpsh')
    doc.appendChild(mpsh)
    xmladd(doc, mpsh, 'nspin', '%i'%DevGF.HS.nspin)
    xmladd(doc, mpsh, 'norbitals', '%i'%(DevGF.nuo))
    xmladd(doc, mpsh, 'energy_values', myprint(options.Elist+DevGF.HS.ef))
    xmladd(doc, mpsh, 'E_Fermi', '%.8f'%DevGF.HS.ef)
    for ii in range(DevGF.nuo):
        orb = doc.createElement('orbital')
        mpsh.appendChild(orb)
        orb.setAttribute('index', '%i'%ii)
        xmladd(doc, orb, 'data', myprint(DOS[:, :, ii]))
    doc.writexml(gzip.GzipFile(fn, 'w'))

    # Make plot
    g = XMGR.Graph()
    for ii in range(DevGF.nuo):
        for iS in range(DevGF.HS.nspin):
            g.AddDatasets(XMGR.XYset(options.Elist, (-1)**iS*DOS[iS, :, ii], legend='', Lwidth=2))

    # Set axes and write XMGR plot to file
    g.SetXaxis(label='E-E\sF\N (eV)', autoscale=True)
    g.SetYaxis(label='DOS (1/eV)', autoscale=True)
    # Add MPSH eigenvalues to plot after axis scaling
    g.AddDatasets(XMGR.XYset(ev0, 0*ev0+1, Ltype=0, Stype=3))
    g.SetTitle(fn, size=1.3)
    g.ShowLegend()
    p = XMGR.Plot(fn+'.xmgr', g)
    p.WriteFile()


def myprint(x): # Do numpy vector or matrix to string
    s = ''
    dim = len(x.shape)
    if dim==1:
        rows, = x.shape
        for i in range(rows):
            s += '%s\n'%x[i]
    if dim==2:
        columns, rows = x.shape
        for i in range(rows):
            for j in range(columns):
                s += '%s '%x[j, i]
            s += '\n'
    return s


def xmladd(doc, parent, name, values):
    # Who came up with xml ... accountant moroons?
    elem = doc.createElement(name)
    parent.appendChild(elem)
    txt=doc.createTextNode(values)
    elem.appendChild(txt)
