print "SVN $Id$"

"""
################################################################

 python TBTrans 
 Magnus Paulsson magnus.paulsson@hik.se

 Requires: numpy (compile it linked with mkl, acml or atlas!)
           ScientificPython (vers. >= 2.8)
           For speed compile the fortran subroutines in F90 
           (cd F90;source compile.bat)

  UNITS! Always eV and Angstrom!
         k-values always given in range [0,1.0] (or [-0.5,0.5])
         They are not in reciprocal space. Instead they corresponds
         to the mathematical orthogonal space that is fourier 
         transformed.

################################################################
"""


import SiestaIO as SIO
import MiscMath as MM
import NEGF
import numpy as N
import Kmesh

def calc(options):
    GaussKronrod = False
    if options.Gk1>1:
        Nk1,t1 = options.Gk1,'GK'
        GaussKronrod = True
    else:
        Nk1,t1 = options.Nk1,'LIN'
    if options.Gk2>1:
        Nk2,t2 = options.Gk2,'GK'
        GaussKronrod = True
    else:
        Nk2,t2 = options.Nk2,'LIN'
    kmesh = Kmesh.kmesh(Nk1,Nk2,Nk3=1,meshtype=[t1,t2,'LIN'],invsymmetry=not options.skipsymmetry)
    NNk = len(kmesh.kpts)
    kPointList = kmesh.kpts[:,:2]
    kWeights = kmesh.wgts[:3]
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL,options.NA1L,options.NA2L,options.voltage/2.)
    elecL.scaling = options.scaleSigL
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR,options.NA1R,options.NA2R,-options.voltage/2.)
    elecR.scaling = options.scaleSigR
    myGF = NEGF.GF(options.TSHS,elecL,elecR,Bulk=options.UseBulk,DeviceAtoms=[options.devSt, options.devEnd])
    nspin = myGF.HS.nspin
    if options.devSt==0:
        options.devSt=GF.DeviceAtoms[0]
    if options.devEnd==0:
        options.devEnd=GF.DeviceAtoms[1]

    print """
##############################################################
pyTBT

Energy [eV]                     : %f:%f:%f
kpoints                         : %i, %i 
eta [eV]                        : %f
etaLead [eV]                    : %f
Device [Atoms Siesta numbering] : %i:%i 
Bulk                            : %s
SpinPolarization                : %i
Voltage                         : %f
##############################################################

"""%(options.Emin,options.dE,options.Emax,Nk1,Nk2,options.eta,options.etaLead,options.devSt,options.devEnd,options.UseBulk,nspin,options.voltage)
        
    if options.dos:
        DOSL=N.zeros((nspin,len(options.Elist),myGF.nuo),N.float)
        DOSR=N.zeros((nspin,len(options.Elist),myGF.nuo),N.float)
    # Loop over spin
    for iSpin in range(nspin):
        Tkpt=N.zeros((len(options.Elist),NNk,options.numchan+1),N.float)
        outFile = options.DestDir+'/%s.%ix%i'%(options.systemlabel,Nk1,Nk2)
        if nspin<2: thisspinlabel = outFile
        else: thisspinlabel = outFile+['.UP','.DOWN'][iSpin]
        fo=open(thisspinlabel+'.AVTRANS','write')
        fo.write('# Nk1(%s)=%i Nk2(%s)=%i eta=%.2e etaLead=%.2e\n'%(t1,Nk1,t2,Nk2,options.eta,options.etaLead))
        if GaussKronrod: fo.write('# E   Ttot(E)   Ti(E) (i=1-10) T_error(E)\n')
        else: fo.write('# E   Ttot(E)   Ti(E) (i=1-10)\n')
        # Loop over energy
        for ie, ee in enumerate(options.Elist):
            Tavg = N.zeros((options.numchan+1,len(kWeights)),N.float)
            AavL = N.zeros((myGF.nuo,myGF.nuo),N.complex)
            AavR = N.zeros((myGF.nuo,myGF.nuo),N.complex)
            # Loops over k-points
            for ik in range(NNk):
                kpt=N.array(kPointList[ik],N.float)
                myGF.calcGF(ee+options.eta*1.0j,kpt,ispin=iSpin,etaLead=options.etaLead,useSigNCfiles=options.signc)
                # Transmission:
                T = myGF.calcT(options.numchan)
                for iw in range(len(kWeights)):
                    Tavg[:,iw] += T*kWeights[iw][ik]
                Tkpt[ie,ik] = T
                # DOS calculation:
                if options.dos:
                    GamL, GamR, Gr = myGF.GamL, myGF.GamR, myGF.Gr
                    nuo, nuoL, nuoR = myGF.nuo, myGF.nuoL, myGF.nuoR
                    AL = MM.mm(Gr[:,0:nuoL],GamL,MM.dagger(Gr)[0:nuoL,:])
                    AR = MM.mm(Gr[:,nuo-nuoR:nuo],GamR,MM.dagger(Gr)[nuo-nuoR:nuo,:])
                    AavL += kWeights[0][ik]*MM.mm(AL,myGF.S)
                    AavR += kWeights[0][ik]*MM.mm(AR,myGF.S)
            # Print calculated quantities
            if GaussKronrod:
                err = (N.abs(Tavg[0,0]-Tavg[0,1])+N.abs(Tavg[0,0]-Tavg[0,2]))/2
                print ee, Tavg[:,0], err
            else:
                print ee, Tavg[:,0]
            
            transline = '\n%.10f '%ee
            for ichan in range(options.numchan+1):
                if ichan==0:
                    transline += '%.8e '%Tavg[ichan,0]
                else:
                    transline += '%.4e '%Tavg[ichan,0]
            if GaussKronrod: transline += '%.4e '%err
            fo.write(transline)
            # Partial density of states:
            if options.dos:
                DOSL[iSpin,ie,:] += N.diag(AavL).real/(2*N.pi)
                DOSR[iSpin,ie,:] += N.diag(AavR).real/(2*N.pi)
                print 'ispin= %i, e= %.4f, DOSL= %.4f, DOSR= %.4f'%(iSpin,ee,N.sum(DOSL[iSpin,ie,:]),N.sum(DOSR[iSpin,ie,:]))
        fo.close()
        
        # Write k-point-resolved transmission
        fo=open(thisspinlabel+'.TRANS','write')
        for ik in range(NNk):
            kpt = kPointList[ik]
            w = kWeights[0][ik]
            fo.write('\n\n# k = %f, %f    w = %f'%(kpt[0],kpt[1],w))
            for ie, ee in enumerate(options.Elist):
                transline = '\n%.10f '%ee
                for ichan in range(options.numchan+1):
                    if ichan==0:
                        transline += '%.8e '%Tkpt[ie,ik,ichan]
                    else:
                        transline += '%.4e '%Tkpt[ie,ik,ichan]
                fo.write(transline)
        fo.close()
    # End loop over spin
    NEGF.SavedSig.close() # Make sure saved Sigma is written to file

    if options.dos:
        # Read basis
        basis = SIO.BuildBasis('%s/%s.XV'%(options.head,options.systemlabel),1,myGF.HS.nua,myGF.HS.lasto)
        WritePDOS(outFile+'.PDOS.gz',options,myGF,DOSL+DOSR,basis)
        WritePDOS(outFile+'.PDOSL.gz',options,myGF,DOSL,basis)
        WritePDOS(outFile+'.PDOSR.gz',options,myGF,DOSR,basis)
    

def WritePDOS(fn,options,myGF,DOS,basis):
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
    devOrbSt = myGF.HS.lasto[options.devSt-1]
    pyTBTdevOrbSt = devOrbSt-myGF.HS.lasto[options.devSt-1]
    devOrbEnd = myGF.HS.lasto[options.devEnd]-1
    pyTBTdevOrbEnd = devOrbEnd-myGF.HS.lasto[options.devSt-1]

    doc = xml.Document()
    pdos = doc.createElement('pdos')
    doc.appendChild(pdos)
    xmladd(doc,pdos,'nspin','%i'%myGF.HS.nspin)
    xmladd(doc,pdos,'norbitals','%i'%(myGF.nuo))
    xmladd(doc,pdos,'energy_values',myprint(options.Elist+myGF.HS.ef))
    xmladd(doc,pdos,'E_Fermi','%.8f'%myGF.HS.ef)
    for ii in range(myGF.nuo):
        orb = doc.createElement('orbital')
        pdos.appendChild(orb)
        io = devOrbSt+ii
        orb.setAttribute('index','%i'%(io+1))
        orb.setAttribute('atom_index','%i'%basis.ii[io])
        orb.setAttribute('species',basis.label[io])
        orb.setAttribute('position','%f %f %f'%(basis.xyz[io,0],basis.xyz[io,1],basis.xyz[io,2]))
        orb.setAttribute('n','%i'%basis.N[io])
        orb.setAttribute('l','%i'%basis.L[io])
        orb.setAttribute('m','%i'%basis.M[io])
        xmladd(doc,orb,'data',myprint(DOS[:,:,ii]))
    doc.writexml(gzip.GzipFile(fn,'w'))

    # Make plot
    atoms = list(set(basis.label))
    lVals  = list(set(basis.L))
    plots = [[atoms,lVals,'Tot']]
    plots += [[atoms,[lVal],'Tot L=%i'%lVal] for lVal in lVals]
    plots += [[[atom],lVals,atom+' Tot'] for atom in atoms]
    plots += [[[atom],[lVal],atom+' L=%i'%lVal] for lVal in lVals for atom in atoms]

    import WriteXMGR as XMGR
    g = XMGR.Graph()
    for atom, lVal, name in plots:
        print atom
        print lVal
        print name
        print [ii for ii in lVal]
        nspin, ee, PDOS = SIO.ExtractPDOS(fn,None,\
                                      FermiRef=False,llist=lVal,\
                                      species=atom)        
        for iS in range(nspin):
            g.AddDatasets(
                XMGR.XYset(ee-myGF.HS.ef,(-1)**iS*PDOS[iS],legend=name,Lwidth=2))

    # Set axes and write XMGR plot to file
    g.SetXaxis(label='E-E\sF\N (eV)',autoscale=True)
    g.SetYaxis(label='DOS (1/eV)',autoscale=True)
    g.SetTitle(fn,size=1.3)
    g.ShowLegend()
    p = XMGR.Plot(fn+'.xmgr',g)
    p.WriteFile()

def myprint(x): # Do numpy vector or matrix to string
    s = ''
    dim = len(x.shape)
    if dim==1:
        rows, = x.shape
        for i in range(rows):
            s += '%s\n'%x[i]
    if dim==2:
        columns,rows = x.shape
        for i in range(rows):
            for j in range(columns):
                s += '%s '%x[j,i]
            s += '\n'
    return s

def xmladd(doc,parent,name,values):
    # Who came up with xml ... accountant moroons?
    elem = doc.createElement(name)
    parent.appendChild(elem)
    txt=doc.createTextNode(values)
    elem.appendChild(txt)

