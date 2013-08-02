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
import numpy.linalg as LA
import sys, string, os
from optparse import OptionParser, OptionGroup
try:
    import scipy.linalg as SLA 
    hasSciPy = True
except:
    hasSciPy = False 
        

def calc(general):
    kPointList, kWeights, NNk, Nk1, Nk2, GaussKronrod = getKpoints(general)
    elecL = NEGF.ElectrodeSelfEnergy(general.fnL,general.NA1L,general.NA2L,general.voltage/2.)
    elecR = NEGF.ElectrodeSelfEnergy(general.fnR,general.NA1R,general.NA2R,-general.voltage/2.)
    myGF = NEGF.GF(general.fnTSHS,elecL,elecR,Bulk=general.UseBulk,DeviceAtoms=[general.devSt, general.devEnd])
    nspin = myGF.HS.nspin
    if general.devSt==0:
        general.devSt=GF.DeviceAtoms[0]
    if general.devEnd==0:
        general.devEnd=GF.DeviceAtoms[1]

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

"""%(general.minE,general.dE,general.maxE,Nk1,Nk2,general.eta,general.etaLead,general.devSt,general.devEnd,general.UseBulk,nspin,general.voltage)
        
    if general.dos:
        DOSL=N.zeros((nspin,len(general.Elist),myGF.nuo),N.float)
        DOSR=N.zeros((nspin,len(general.Elist),myGF.nuo),N.float)
    # Loop over spin
    for iSpin in range(nspin):
        Tkpt=N.zeros((len(general.Elist),NNk,general.numchan+1),N.float)
        outFile = general.DestDir+'/%s.%ix%i'%(general.systemlabel,Nk1,Nk2)
        if nspin<2: thisspinlabel = outFile
        else: thisspinlabel = outFile+['.UP','.DOWN'][iSpin]
        fo=open(thisspinlabel+'.AVTRANS','write')
        fo.write('# Nk1=%i Nk2=%i eta=%.2e etaLead=%.2e\n'%(Nk1,Nk2,general.eta,general.etaLead))
        if GaussKronrod: fo.write('# E   Ttot(E)   Ti(E) (i=1-10) T_error(E)\n')
        else: fo.write('# E   Ttot(E)   Ti(E) (i=1-10)\n')
        # Loop over energy
        for ie, ee in enumerate(general.Elist):
            Tavg = N.zeros((general.numchan+1,len(kWeights)),N.float)
            AavL = N.zeros((myGF.nuo,myGF.nuo),N.complex)
            AavR = N.zeros((myGF.nuo,myGF.nuo),N.complex)
            # Loops over k-points
            for ik in range(NNk):
                kpt=N.array(kPointList[ik],N.float)
                myGF.calcGF(ee+general.eta*1.0j,kpt,ispin=iSpin,etaLead=general.etaLead,useSigNCfiles=general.signc)
                # Transmission:
                T = myGF.calcT(general.numchan)
                for iw in range(len(kWeights)):
                    Tavg[:,iw] += T*kWeights[iw][ik]
                Tkpt[ie,ik] = T
                # DOS calculation:
                if general.dos:
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
            for ichan in range(general.numchan+1):
                if ichan==0:
                    transline += '%.8e '%Tavg[ichan,0]
                else:
                    transline += '%.4e '%Tavg[ichan,0]
            if GaussKronrod: transline += '%.4e '%err
            fo.write(transline)
            # Partial density of states:
            if general.dos:
                DOSL[iSpin,ie,:] += N.diag(AavL).real/(2*N.pi)
                DOSR[iSpin,ie,:] += N.diag(AavR).real/(2*N.pi)
                print 'ispin= %i, e= %.4f, DOSL= %.4f, DOSR= %.4f'%(iSpin,ee,N.sum(DOSL[iSpin,ie,:]),N.sum(DOSR[iSpin,ie,:]))
        fo.close()
        
        # Write k-point-resolved transmission
        fo=open(thisspinlabel+'.TRANS','write')
        for ik in range(NNk):
            kpt=kPointList[ik]
            fo.write('\n\n# k = %f, %f '%(kpt[0],kpt[1]))
            for ie, ee in enumerate(general.Elist):
                transline = '\n%.10f '%ee
                for ichan in range(general.numchan+1):
                    if ichan==0:
                        transline += '%.8e '%Tkpt[ie,ik,ichan]
                    else:
                        transline += '%.4e '%Tkpt[ie,ik,ichan]
                fo.write(transline)
        fo.close()
    # End loop over spin
    NEGF.SavedSig.close() # Make sure saved Sigma is written to file

    if general.dos:
        WritePDOS(outFile+'.PDOS.gz',general,myGF,DOSL+DOSR)
        WritePDOS(outFile+'.PDOSL.gz',general,myGF,DOSL)
        WritePDOS(outFile+'.PDOSR.gz',general,myGF,DOSR)
    

def WritePDOS(fn,general,myGF,DOS):
    """
    PDOS from the surface Green's function from electrode calculations
    
    NOTE! The DOS is a sum over the atoms of the unitcell.
    NOTE! The outfile contains the DOS divided into s,p,d,f shells.
    This decomposition is not perfect since polarized basis orbitals
    will end up in L+1 part, i.e., 6s polarized orbital = 6p
    """

    import xml.dom.minidom as xml
    import gzip

    # Read basis
    basis=SIO.BuildBasis(general.systemlabel+'.XV',1,myGF.HS.nua,myGF.HS.lasto)
    
    # First, last orbital in full space and pyTBT folded space.
    devOrbSt = myGF.HS.lasto[general.devSt-1]
    pyTBTdevOrbSt = devOrbSt-myGF.HS.lasto[general.devSt-1]
    devOrbEnd = myGF.HS.lasto[general.devEnd]-1
    pyTBTdevOrbEnd = devOrbEnd-myGF.HS.lasto[general.devSt-1]

    doc = xml.Document()
    pdos = doc.createElement('pdos')
    doc.appendChild(pdos)
    xmladd(doc,pdos,'nspin','%i'%myGF.HS.nspin)
    xmladd(doc,pdos,'norbitals','%i'%(myGF.nuo))
    xmladd(doc,pdos,'energy_values',myprint(general.Elist+myGF.HS.ef))
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

def getKpoints(general):
    # Do 1-D k-points
    if general.Gk1>1 and general.Gk2>1: # Method fails with fewer points
        GaussKronrod=True
        kl1, kwl1, kwle1 = MM.GaussKronrod(general.Gk1)
        kl2, kwl2, kwle2 = MM.GaussKronrod(general.Gk2)        
    else:
        GaussKronrod=False
        kl1 = [(ii*1.0)/general.Nk1-0.5+0.5/general.Nk1 for ii in range(general.Nk1)]
        kwl1 = [1.0/general.Nk1 for ii in range(general.Nk1)]
        kl2 = [(ii*1.0)/general.Nk2-0.5+0.5/general.Nk2 for ii in range(general.Nk2)]
        kwl2 = [1.0/general.Nk2 for ii in range(general.Nk2)]
        kl1, kwl1, kl2, kwl2 = N.array(kl1), N.array(kwl1), N.array(kl2), N.array(kwl2)

    # Shift away from gamma for normal k-points
    if general.skipgamma and not GaussKronrod:
        if general.Nk1%2==1: kl1 += 0.5/general.Nk1
        if general.Nk2%2==1: kl2 += 0.5/general.Nk2

    # Repeat out for 2D
    kl, kwl = N.zeros((len(kl1)*len(kl2),2)), N.zeros((len(kl1)*len(kl2),))
    if GaussKronrod:
        TDkwle1, TDkwle2 = N.zeros((len(kl1)*len(kl2),)), N.zeros((len(kl1)*len(kl2),))
    jj=0
    for i1 in range(len(kl1)):
        for i2 in range(len(kl2)):
            kl[jj,:]=[kl1[i1],kl2[i2]]
            kwl[jj]=kwl1[i1]*kwl2[i2]
            if GaussKronrod:
                TDkwle1[jj]=kwle1[i1]*kwl2[i2]
                TDkwle2[jj]=kwl1[i1]*kwle2[i2]
            jj+=1

    # Remove duplicates for symmetry
    # INVERSION SYMMETRY:
    # If the Bloch function
    #    \psi(k) = exp(ikr)u(k),
    # with crystal momentum k, is an eigenstate of the Schroedinger equation then also
    #    \psi^\dagger(k) = exp(-ikr)u^\dagger(k)
    # with crystal momentum -k, is an eigenstate with same eigenvalue.
    # Hence E(k) = E(-k).
    # TIME REVERSAL SYMMETRY:
    # t,\psi(r,t) --> -t,\psi^\dagger(r,-t). T(k) = T(-k).
    # (Elastic) propagation from L to R is always identical to propagation from R to L.
    if not general.skipsymmetry:
        print 'pyTBT: Applying inversion (time-reversal) symmetry reduction to list of k-points'
        indx = []
        for i1 in range(len(kl)):
            for i2 in range(len(indx)):
                if N.allclose(-kl[i1],kl[indx[i2][0]],atol=1e-7,rtol=1e-7):
                    indx[i2][1]+=1
                    break
            else:
                indx+=[[i1,1]]
        indx, weight   = N.array([ii[0] for ii in indx]), N.array([ii[1] for ii in indx])
        kl, kwl = kl[indx], kwl[indx]*weight
        if GaussKronrod:
            TDkwle1, TDkwle2 = TDkwle1[indx]*weight, TDkwle2[indx]*weight
    
    s = 'pyTBT.getKpoints: i, k1[i], k2[i], kwl[i]'
    if GaussKronrod: s += ', TDkwle1[i], TDkwle2[i]'
    print s
    for i in range(len(kl)):
        s = '... %i   %.8f %.8f   %.8f'%(i,kl[i,0],kl[i,1],kwl[i])
        if GaussKronrod:
            s += '   %.8e %.8e'%(TDkwle1[i],TDkwle2[i])
        print s
    print 'Nk = %i, Sum(kw) = %.4f\n'%(len(kl),N.sum(kwl))
    
    if GaussKronrod:
        return kl, [kwl, TDkwle1, TDkwle2], len(kl), general.Gk1, general.Gk1, GaussKronrod
    else:
        return kl, [kwl], len(kl), general.Nk1, general.Nk1, GaussKronrod


