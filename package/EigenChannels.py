"""
Eigenchannels:
1: Eigenchannels, method from Paulsson and Brandbyge PRB 2007
2: Analyze IETS spectra, Paulsson et al PRL 2008
3: Calculate "bond" currents
"""
fudgeEnergy = 5e-5         # Fudge energy in eV used to see if Ef is close to the calculated energy

import pyTBT
import SiestaIO as SIO
import MakeGeom as MG
#import EigenChannels as EC
import WriteXMGR as XMGR
import numpy as N
import numpy.linalg as LA
import Scientific.IO.NetCDF as NC
import sys, string, struct, glob,  os
from optparse import OptionParser, OptionGroup
import PhysicalConstants as PC

# Dummy classes for global variables
class HS:       pass
class IETS:     pass
# general, basis and geom are a global variables

########################################################
##################### Main routine #####################
########################################################
def main():
    setupParameters()

    readxv()
    readHS()
    readbasis()
    calcT()

########################################################
def readxv():
    global geom
    # Read geometry from first .XV file found in dir
    fns=glob.glob('*.XV')

    if len(fns)>1:
        print "ERROR: Eigenchannels: More than one .XV file ... which geometry to choose???"
        sys.exit(1)
    elif len(fns)<1:
        print "ERROR: Eigenchannels: Error ... No .XV file found!"
        sys.exit(1)

    print('Reading geometry from "%s" file' % fns[0])
    geom = MG.Geom(fns[0])

########################################################
def readbasis():
    global basis
    fn=glob.glob('*.XV')
    basis = SIO.BuildBasis(fn[0], general.from_atom, general.to_atom, HS.GF.HS.lasto)

########################################################
def readHS():
    # Setup H, S and self-energies
    HS.L, HS.R, HS.GF, devSt, devEnd, general.Elist, general.eta, general.SiestaOutFile = \
        pyTBT.main(pyTBT=False, fn=general.fdfFile, 
                   deviceRegion=[general.from_atom, general.to_atom])
    if general.PhononNetCDF!=None:
        # Check device region so that it crashes fast!
        NNTS = len(HS.GF.H)
        NCfile = NC.NetCDFFile(general.PhononNetCDF,'r')
        NNPH = len(N.array(NCfile.variables['S0'][:,:]))
        NCfile.close()
        if NNPH!=NNTS:
            print "ERROR: Eigenchannels : ERROR the device regions of pyTBT and PhononNetCDF do not match!!!"
            sys.exit(1)

########################################################
def calcT(calledFromInelastica=False,InelasticaEf=0.0):
    # Calculate eigenchannels for specific k-point
    if calledFromInelastica:
        E=N.array([InelasticaEf],N.float)
    else:
        E = general.Elist
    HS.GF.calcGF(E[0]+general.eta*1j, general.kPoint[0:2], ispin=general.iSpin)

    # Orthogonalize everything
    S=HS.GF.S.copy()
    Usi=mysqrt(S)
    Us=LA.inv(Usi)

    S, H = mm(Us,S,Us), mm(Us,HS.GF.H,Us)


    # Matrix to save total and eigenchannel transmissions
    T=N.zeros((len(E),1+general.numchan),N.complex)
    IETS.T = T

    PDOSL, PDOSR = N.zeros((len(E),len(H)),N.float), N.zeros((len(E),len(H)),N.float)

    # Loop over energy points
    for iiE in range(len(E)):
        SIO.printDone(iiE,len(E),'Transmission') 

        # Setup the normal quantites
        EE=E[iiE]
        if iiE!=0:
            # Already done for first point, see above
            HS.GF.calcGF(EE+general.eta*1j, general.kPoint[0:2], 
                         ispin=general.iSpin)
        
        # Gammas in pyTBT GF can be smaller than device region
        nnL,  nnR = len(HS.GF.GamL),  len(HS.GF.GamR)
        G1,  G2 = N.zeros(HS.GF.H.shape, N.complex), N.zeros(HS.GF.H.shape, N.complex)
        G1[0:nnL, 0:nnL] ,  G2[-nnR:, -nnR:] = HS.GF.GamL, HS.GF.GamR
        G1, G2 = mm(Us,G1,Us), mm(Us,G2,Us)

        G=mm(Usi,HS.GF.Gr,Usi)
        Gd=dagger(G)

        A1, A2 =mm(G,G1,Gd), mm(G,G2,Gd)

        T[iiE,0]=HS.GF.calcT()[0]

        if not calledFromInelastica:
            # Calculate PDOS from left and right
            A1NO, A2NO = mm(Us,A1,Us), mm(Us,A2,Us)
            for ij in range(len(H)):
                PDOSL[iiE,ij]=abs(A1NO[ij,ij])
                PDOSR[iiE,ij]=abs(A2NO[ij,ij])

        # Calculate total "bond currents"
        if abs(EE-general.energy)<fudgeEnergy and not calledFromInelastica:
            Curr=-calcCurrent(A1NO)
            general.iChan, general.iSide = 0, 0
            writeCurrent(Curr)
            Curr=-calcCurrent(A2NO)
            general.iSide = 1
            writeCurrent(Curr)

        if not calledFromInelastica:
            # Calculate Eigenchannels from left
            EigC, EigT = calcEigChan(A1,G2,Us)
            for jj in range(general.numchan):
                T[iiE,1+jj]=EigT[len(EigT)-jj-1]

        # Save to file
        if abs(EE-general.energy)<fudgeEnergy and not calledFromInelastica:
            ECleft=EigC
            # Write wave functions
            for jj in range(general.numchan):
                general.iSide, general.iChan = 0, jj+1
                writeWavefunction(ECleft[jj])
                Curr=calcCurrent(ECleft[jj])
                writeCurrent(Curr)
                    
        # Calculate eigenchannels from right
        if abs(EE-general.energy)<fudgeEnergy and (general.bothsides or general.PhononNetCDF!=None):
            if not calledFromInelastica:
                ECright, EigT = calcEigChan(A2,G1,Us)
            
            if general.bothsides and not calledFromInelastica:
                # Write wave functions
                for jj in range(general.numchan):
                    general.iSide, general.iChan = 1, jj+1
                    writeWavefunction(ECright[jj])
                    Curr=calcCurrent(ECright[jj])
                    writeCurrent(Curr)
            
            # Calculate inelastic scattering rates
            if general.PhononNetCDF!=None:
                calcIETS(T,G,Gd,G1,G2,A1,A2,Us,Usi,iiE,calledFromInelastica)
                if not calledFromInelastica:
                    IETS.ECleft, IETS.ECright = ECleft, ECright
                    writeFGRrates()

    # Calculate eigenstates of device Hamiltonian
    if general.MolStates>0.0:
        ev, es = LA.eigh(H)
        for ii in range(len(ev)):
            if N.abs(ev[ii])<general.MolStates:
                fn=general.DestDir+'/'+general.SiestaOutFile+'.ES.%i.%2.2f'%(ii,ev[ii])
                writeWavefunction(mm(Us,es[:,ii]),fn=fn)

    if not calledFromInelastica:
        writeTrans(T,H,E)

########################################################
def calcEigChan(A1,G2,Us):
    # Calculate Eigenchannels using recipie from PRB
    # For right eigenchannels, A1=A2, G2=G1 !!!
    ev, U = LA.eigh(A1)
    U = N.transpose(U)

    Utilde=0.0*U

    for jj in range(len(ev)): # Problems with negative numbers
        if ev[jj]<0: ev[jj]=0
        Utilde[jj,:]=N.sqrt(ev[jj]/(2*N.pi))*U[jj,:]
    Utilde=N.transpose(Utilde)
        
    tt=mm(dagger(Utilde),2*N.pi*G2,Utilde)
    evF, UF = LA.eigh(tt)
    UF = N.transpose(UF)

    NumChans=max(N.sum(evF>1e-6),general.numchan)

    EC=[]
    for jj in range(NumChans):
        tmp=mm(Us,Utilde,UF[len(evF)-jj-1,:])
        EC.append(tmp.copy())

    return EC, evF

########################################################
def calcIETS(T,G,Gd,G1,G2,A1,A2,Us,Usi,iiE,calledFromInelastica=False):                
    # Calculate inelastic scattering rates, total and per eigenchannel
    IETS.totTrans = T[iiE,0]

    NCfile = NC.NetCDFFile(general.PhononNetCDF,'r')
    print 'Reading ',general.PhononNetCDF
    hw = N.array(NCfile.variables['hw'][:]) 
    bA1 = mm(Gd,G1,G) 
    
    # Structures to save Power, e-h damping, (non-)Hilbert term
    P1T,P2T = [],[]
    ehDamp1, ehDamp2 = [], []
    nHT,nHTin,nHTel,HT = [],[],[],[]
    
    # Save bond current changes due to inelastic scattering 
    dGnout, dGnin = [], []
    
    for ihw in range(len(hw)):
        SIO.printDone(ihw,len(hw),'IETS FGR')
        M = N.array(NCfile.variables['He_ph'][ihw,general.iSpin,:,:])
        H1=mm(Us,M,Us)
        MA1M, MA2M = mm(H1,A1,H1), mm(H1,A2,H1)
        MAM, MA2mA1M = MA1M+MA2M, MA2M-MA1M
        
        # Cryptic? Changes in G^lesser due to e-ph interaction at high and low energies, i.e.,
        # the changes in occupation due to out(in)-scattering
        tmp1, tmp2 = mm(G,MA2M,A1), mm(G,MA1M,A1)
        if not calledFromInelastica:
            dGnout.append(calcCurrent(mm(Us,-0.5j*(tmp1-dagger(tmp1)),Us)))
            dGnin.append(calcCurrent(mm(Us,mm(G,MA1M,Gd)-0.5j*(tmp2-dagger(tmp2)),Us)))

        # Power, damping and current rates
        P1T.append(checkImPart(N.trace(mm(MAM,A1+A2))))
        P2T.append(checkImPart(N.trace(mm(MA1M,A2))))        
        ehDamp1.append(checkImPart(N.trace(mm(MA1M,A1))))
        ehDamp2.append(checkImPart(N.trace(mm(MA2M,A2))))
        nHT.append(checkImPart(N.trace(mm(-2*MA2M+1.0j*(mm(MAM,G,G2)-mm(G2,Gd,MAM)),bA1))/2.0))
        nHTin.append(checkImPart(N.trace(mm(-2*MA2M,bA1))/2.0))
        nHTel.append(checkImPart(N.trace(mm(1.0j*(mm(MAM,G,G2)-mm(G2,Gd,MAM)),bA1))/2.0))
        HT.append(checkImPart(N.trace(mm((mm(G2,Gd,MA2mA1M)+mm(MA2mA1M,G,G2)),bA1))))
    NCfile.close()

    IETS.dGnin, IETS.dGnout = dGnin, dGnout
    IETS.P1T, IETS.P2T = P1T, P2T
    IETS.ehDamp1, IETS.ehDamp2 = ehDamp1, ehDamp2
    IETS.nHT, IETS.HT = nHT, HT
    IETS.nHTin, IETS.nHTel = nHTin, nHTel
    IETS.hw = hw

########################################################
def checkImPart(x):
    if abs(x.imag)>0.0000001:
        print "LOE : Imaginary part (%.3e) too big"%x.imag
        kuk
    return x.real   

########################################################
def writeFGRrates():
    unitConv=1.602177e-19/N.pi/1.054572e-34

    NCfile = NC.NetCDFFile(general.PhononNetCDF,'r')
    print 'Reading ',general.PhononNetCDF

    general.iChan, general.iSide = 0, 2
    fn = fileName()
    outFile = file(fn+'.FGR.out','w')
    outFile.write('Total transmission [in units of (1/s/eV)] : %e\n' % (unitConv*IETS.totTrans.real,))

    tmp=N.sort(abs(N.array(IETS.nHT[:])))
    SelectionMin=tmp[-general.NumPhCurr]        
    
    for ihw in range(len(IETS.hw)):
        SIO.printDone(ihw,len(IETS.hw),'Golden Rate') 
        M = N.array(NCfile.variables['He_ph'][ihw,general.iSpin,:,:])
        rate=N.zeros((len(IETS.ECleft),len(IETS.ECright)),N.float)
        totrate=0.0
        for iL in range(len(IETS.ECleft)):
            for iR in range(len(IETS.ECright)):
                tmp=N.dot(N.conjugate(IETS.ECleft[iL]),mm(M,IETS.ECright[iR]))
                rate[iL,iR]=(2*N.pi)**2*abs(tmp)**2
                totrate+=rate[iL,iR]

        if abs(IETS.nHT[ihw])>=SelectionMin:
            general.iChan = ihw
            currOut, currIn = IETS.dGnout[ihw], IETS.dGnin[ihw]
            general.iSide = 3
            writeCurrent(currIn)
            general.iSide = 4
            writeCurrent(currOut)
        
        outFile.write('\nPhonon mode %i : %f eV [Rates in units of (1/s/eV)]\n' % (ihw,IETS.hw[ihw]))
        outFile.write('eh-damp : %e (1/s) , heating %e (1/(sV)))\n' % (IETS.P1T[ihw]*unitConv*IETS.hw[ihw],IETS.P2T[ihw]*unitConv))
        outFile.write('eh-damp 1, 2 (MA1MA1, MA2MA2): %e (1/s) , %e (1/(s)))\n' % (IETS.ehDamp1[ihw]*unitConv*IETS.hw[ihw],IETS.ehDamp2[ihw]*unitConv*IETS.hw[ihw]))
        outFile.write('SymI : %e (1/(sV)) , AsymI %e (?))\n' % (IETS.nHT[ihw]*unitConv,IETS.HT[ihw]*unitConv))
        outFile.write('Elast : %e (1/(sV)) , Inelast %e (1/(sV)))\n' % (IETS.nHTel[ihw]*unitConv,IETS.nHTin[ihw]*unitConv))
        outFile.write('down=left EC, right=right EC\n')
        if IETS.P2T[ihw]>0.0:
            if abs(totrate/(IETS.P2T[ihw])-1)<0.05:
                outFile.write('Sum/Tr[MA1MA2] , Tr: %1.3f  %e\n'%(totrate/(IETS.P2T[ihw]),unitConv*IETS.P2T[ihw]))
            else:
                outFile.write('WARNING: !!!! Sum/Tr[MA1MA2] , Tr: %2.2e  %e\n'%(totrate/(IETS.P2T[ihw]),unitConv*IETS.P2T[ihw]))
        else:
            outFile.write(' Tr:  %e\n'%(unitConv*IETS.P2T[ihw]))
        for iL in range(len(IETS.ECleft)):
            for iR in range(len(IETS.ECright)):
                outFile.write('%e ' % (unitConv*rate[iL,iR],))
            outFile.write('\n')
    outFile.close()
    NCfile.close()
    
########################################################
def TODO():
    # Divide PDOS into S, P, D orbitals ... TODO! Crap with selecting interesting region
    PDOS2=N.zeros((len(E),3),N.float)    
    for ii in range(len(E)):
        for jj in range(len(H)):
            if basis.atomnum[jj]==27:
                if basis.M[jj]==0:
                    PDOS2[ii,basis.L[jj]]+=PDOS[ii,jj]

    for jj in range(3):
        tmp=['S','P','D']
        foS=file('Co.'+tmp[jj]+'.txt','w')
        for ii in range(len(E)):
            foS.write('%f %f\n'%(E[ii],PDOS2[ii,jj]))
        foS.close()

    PDOS2=N.zeros((len(E),3),N.float)
    
    for ii in range(len(E)):
        for jj in range(len(H)):
            if basis.atomnum[jj]==79:
                PDOS2[ii,basis.L[jj]]+=PDOS[ii,jj]

    for jj in range(3):
        tmp=['S','P','D']
        foS=file('Au.'+tmp[jj]+'.txt','w')
        for ii in range(len(E)):
            foS.write('%f %f\n'%(E[ii],PDOS2[ii,jj]))
        foS.close()

########################################################
def writeTrans(T,H,E):
    T = N.transpose(T)

    # Make plot
    g = XMGR.Graph()
    
    # Add total and eigenchannel transmissions
    g.AddDatasets(
        XMGR.XYset(E,T[0].real,legend='T\stot',Lwidth=2)
        )
    for i in range(1,len(T)):
        g.AddDatasets(
            XMGR.XYset(E,T[i].real,legend='T\s%i'%i)
            )

    # Add eigenvalues of the periodic system
    ev, U = LA.eigh(H)
    first,last = 0,0
    for jj in range(len(ev)):
        if ev[jj] <= min(E):
            first = jj
        if ev[jj] <= max(E):
            last = jj
    g.AddDatasets(
        XMGR.XYset(ev[first:last+1],0.0*ev[first:last+1],
                   legend='EV',Stype=4,Lstyle=0)
        )

    # Set axes and write XMGR plot to file
    KP = general.kPoint
    fn = general.DestDir+'/'+general.SiestaOutFile \
         + ('.KP_%.3f_%.3f_%.3f'%(KP[0],KP[1],KP[2]) ) \
         + '.TRANS'
    g.SetXaxis(label='E-E\sF\N (eV)',autoscale=True)
    g.SetYaxis(label='Transmission',autoscale=True)
    g.SetTitle(fn,size=1.3)
    g.ShowLegend()
    p = XMGR.Plot(fn+'.xmgr',g)
    p.WriteFile()

    # Write ASCII data file
    f = open(fn,'w')
    for i in range(len(E)):
        f.write('\n%.6f   '%E[i])
        for j in range(len(T)):
            f.write('%.6f '%T[j][i].real)
    f.close()


########################################################
def calcWF(Y):
    """
    Calculate wavefunction, returns:
    YY : complex wavefunction on regular grid
    dstep : stepsize
    origo : vector
    nx, ny, nz : number of grid points
    """

    xyz=N.array(geom.xyz[general.from_atom-1:general.to_atom])
    atomnum=geom.anr[general.from_atom-1:general.to_atom]

    # Size of cube
    xmin, xmax = min(xyz[:,0])-5.0, max(xyz[:,0])+5.0
    ymin, ymax = min(xyz[:,1])-5.0, max(xyz[:,1])+5.0
    zmin, zmax = min(xyz[:,2])-5.0, max(xyz[:,2])+5.0
    xl, yl, zl = xmax-xmin, ymax-ymin, zmax-zmin
    dx, dy, dz = general.res, general.res, general.res
    nx, ny, nz = int(xl/dx)+1, int(yl/dy)+1, int(zl/dz)+1

    origo = N.array([xmin,ymin,zmin],N.float)

    # Def cube
    YY=N.zeros((nx,ny,nz),N.complex)
    rx=N.array(range(nx),N.float)*dx+origo[0]
    ry=N.array(range(ny),N.float)*dy+origo[1]
    rz=N.array(range(nz),N.float)*dz+origo[2]

    for ii in range(len(Y)):
        if ii>0 and ii%(int(len(Y)/10))==0:
            SIO.printDone(ii,len(Y),'Wavefunction')  

        rax,ray,raz=basis.xyz[ii,0],basis.xyz[ii,1],basis.xyz[ii,2]
        # Only calulate in subset
        ixmin, ixmax = int((rax-origo[0]-basis.coff[ii])/dx), \
                       int((rax-origo[0]+basis.coff[ii])/dx)
        iymin, iymax = int((ray-origo[1]-basis.coff[ii])/dy), \
                       int((ray-origo[1]+basis.coff[ii])/dy)
        izmin, izmax = int((raz-origo[2]-basis.coff[ii])/dz), \
                       int((raz-origo[2]+basis.coff[ii])/dz)
        
        ddx,ddy,ddz=rx[ixmin:ixmax]-rax,ry[iymin:iymax]-ray,rz[izmin:izmax]-raz

        dr=N.sqrt(outerAdd(ddx*ddx,ddy*ddy,ddz*ddz))
        drho=N.sqrt(outerAdd(ddx*ddx,ddy*ddy,0*ddz))
        
        imax=(basis.coff[ii]-2*basis.delta[ii])/basis.delta[ii]
        ri=dr/basis.delta[ii]            
        ri=N.where(ri<imax,ri,imax)
        ri=ri.astype(N.int)
        costh, sinth = outerAdd(0*ddx,0*ddy,ddz)/dr, drho/dr
        cosfi, sinfi = outerAdd(ddx,0*ddy,0*ddz)/drho, outerAdd(0*ddx,ddy,0*ddz)/drho

        # Numpy has changed the choose function to crap!
        RR=N.take(basis.orb[ii],ri)

        SphHar=sphericalHarmonic(ii)

        YY[ixmin:ixmax,iymin:iymax,izmin:izmax]=YY[ixmin:ixmax,iymin:iymax,izmin:izmax]+\
                                                 RR*SphHar(sinth,costh,sinfi,cosfi)*Y[ii]

    return YY, general.res, origo, nx, ny, nz

###################### Current densities #########################

def calcCurrent(Y):
    """
    Calculate current density in atomic bonds
    Y : complex scattering state or
    Y : A_l or A_r! (for total current)
    """
    
    NN=len(HS.GF.H)
    NN2=general.to_atom-general.from_atom+1
    Curr=N.zeros((NN2,NN2),N.float)
    
    if len(Y.shape)==2:
        for ii in range(NN):
            a1=basis.ii[ii]-general.from_atom
            for jj in range(NN):
                a2=basis.ii[jj]-general.from_atom
                tmp=HS.GF.H[jj,ii]*Y[ii,jj]/2/N.pi
                Curr[a1,a2]=Curr[a1,a2]+4*N.pi*tmp.imag
    else:
        for ii in range(NN):
            a1=basis.ii[ii]-general.from_atom
            for jj in range(NN):
                a2=basis.ii[jj]-general.from_atom
                tmp=HS.GF.H[ii,jj]*N.conjugate(Y[ii])*Y[jj]
                Curr[a1,a2]=Curr[a1,a2]+4*N.pi*tmp.imag
    
    return Curr

########################################################
def writeCurrent(Curr):
    """
    Write bond currents to file
    fn : Filename
    Curr : Bond current matrix
    """
    fn=fileName()
    xyz=N.array(geom.xyz[general.from_atom-1:general.to_atom])
    atomnum=geom.anr[general.from_atom-1:general.to_atom]

    foC=file(fn+'.curr','w')
    foC.write('%i\n'%(general.to_atom-general.from_atom+1))

    for ii in range(len(xyz)):
        foC.write('%i %e %e %e\n'%(atomnum[ii],xyz[ii,0],xyz[ii,1],xyz[ii,2]))
        
    for ii in range(len(Curr)):
        for jj in range(len(Curr)):
            foC.write('%e '%(Curr[ii,jj]))
        foC.write('\n')
    foC.close()

########################################################
def writenetcdf(fn,YY,nx,ny,nz,origo,dstep):
    """
    THF: Write eigenchannels to netcdf format
    """
    import time
    file = NC.NetCDFFile(fn,'w','Created '+time.ctime(time.time()))
    file.createDimension('nx',nx)
    file.createDimension('ny',ny)
    file.createDimension('nz',nz)
    file.createDimension('natoms',len(geom.xyz))
    file.createDimension('naxes',3)
    file.createDimension('number',1)
    #file.createDimension('pair',2)
        
    # Grid
    #grid = file.createVariable('grid','d',('naxes','pair'))
    #tmp  = [origo,[dstep,dstep,dstep]]
    #grid[:] = N.transpose(N.array(tmp))

    # Fields
    varRe = file.createVariable('Re-Psi','d',('nx','ny','nz'))
    varRe[:] = YY.real
    
    varIm = file.createVariable('Im-Psi','d',('nx','ny','nz'))
    varIm[:] = YY.imag

    vardstep = file.createVariable('dstep','d',('number',))
    vardstep[:]  = dstep

    vargeom = file.createVariable('xyz','f',('natoms','naxes'))
    # OpenDX needs float for positions
    tmp = []
    for i in range(len(geom.xyz)):
        tmp2 = (geom.xyz[i]-origo)/dstep
        tmp.append([float(tmp2[0]),
                    float(tmp2[1]),
                    float(tmp2[2])])
    vargeom[:] = tmp

    varanr = file.createVariable('anr','i',('natoms',))
    varanr[:] = geom.anr
    
    # Set attributes
    setattr(varanr,'field','anr')
    setattr(varanr,'positions','xyz')

    #setattr(varRe,'field','Re-Psi')
    #setattr(varRe,'positions','grid, compact')

    #setattr(varIm,'field','Im-Psi')
    #setattr(varIm,'positions','grid, compact')

    file.close()
    
################# Write cube file ######################
def writecube(fn,YY,nx,ny,nz,origo,dstep):
    """
    Write wavefunction to cube file.
    """
    global geom

    xyz=N.array(geom.xyz)
    anr=geom.anr

    foR=file(fn,'w')
    foR.write('Eigenchannel wavefunction\n%s\n'%fn)
    foR.write('%i %f %f %f\n'% (len(xyz),origo[0]/PC.Bohr2Ang,origo[1]/PC.Bohr2Ang,origo[2]/PC.Bohr2Ang))
    foR.write('%i %f %f %f\n'% (nx,dstep/PC.Bohr2Ang,0.0,0.0))
    foR.write('%i %f %f %f\n'% (ny,0.0,dstep/PC.Bohr2Ang,0.0))
    foR.write('%i %f %f %f\n'% (nz,0.0,0.0,dstep/PC.Bohr2Ang))
    # Write atom coordinates
    for ii in range(len(xyz)):
        foR.write('%i %f '% (anr[ii],0.0))
        tmp=xyz[ii,:]/PC.Bohr2Ang
        foR.write('%f %f %f\n'% (tmp[0],tmp[1],tmp[2]))
    # Write wavefunction
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                foR.write('%1.3e\n' % YY[ix,iy,iz].real)
    foR.close()

################# Write wave function in macu format ######################
def writemacubin(fn,YY,nx,ny,nz,origo,dstep):
    """
    Write molekel binary format
    """
    
    fo=file(fn,'w')
    fo.write(struct.pack('i',36))
    xmin,xmax = origo[0], origo[0]+dstep*(nx-1)
    ymin,ymax = origo[1], origo[1]+dstep*(ny-1)
    zmin,zmax = origo[2], origo[2]+dstep*(nz-1)
    fo.write(struct.pack('6f4i',
                         xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,nx*ny*nz))
    
    for ii in range(nz):
        bin=struct.pack('i',nx*ny*4)
        for kk in range(ny):
            for jj in range(nx):
                bin+=struct.pack('f',YY[jj,kk,ii])
        bin+=struct.pack('i',nx*ny*4)
        fo.write(bin)

    fo.close()

################ Write wave function in XSF format ##################################
def writeXSF(fn,YY,nx,ny,nz,origo,dstep):
    """
    Write XSF datagrid for XCrysden
    """
    fo=file(fn,'w')
    vectors=geom.pbc
    speciesnumber=geom.snr
    atomnumber=geom.anr
    xyz=geom.xyz
    xmin,xmax = origo[0], origo[0]+dstep*(nx-1)
    ymin,ymax = origo[1], origo[1]+dstep*(ny-1)
    zmin,zmax = origo[2], origo[2]+dstep*(nz-1)
    
    fo.write(' ATOMS\n')
    numberOfAtoms = len(speciesnumber)

    # Write out the position for the atoms
    for i in range(numberOfAtoms):
        line = '   %i  ' %atomnumber[i]
        for j in range(3):
            line += string.rjust('%.9f'%xyz[i][j],16)
        fo.write(line+'\n')

    #Write the datagrid
    fo.write('BEGIN_BLOCK_DATAGRID_3D\n')
    fo.write(' DATA_from:Inelastica\n')
    # Real part
    fo.write(' BEGIN_DATAGRID_3D_REAL\n')
    fo.write('    %3.0i    %3.0i    %3.0i\n'%(nx,ny,nz))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (origo[0],origo[1],origo[2]))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (xmax-xmin,0.0000,0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000,ymax-ymin,0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000,0.0000,zmax-zmin))
    data=[]
    ll=0
    for ii in range(nz):
        for kk in range(ny):
            for jj in range(nx):
    		data.append(YY.real[jj,kk,ii])
    for iii in range((nx*ny*nz)) :
	if ((iii+1)%6==0) :	
            fo.write('  %1.5E\n'% (data[iii]))
	else :
            fo.write('  %1.5E'% (data[iii]))
    fo.write('\n END_DATAGRID_3D\n')
    # Imaginary part
    fo.write(' BEGIN_DATAGRID_3D_IMAG\n')
    fo.write('    %3.0i    %3.0i    %3.0i\n'%(nx,ny,nz))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (origo[0],origo[1],origo[2]))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (xmax-xmin,0.0000,0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000,ymax-ymin,0.0000))
    fo.write('  %1.7E  %1.7E  %1.7E\n'% (0.0000,0.0000,zmax-zmin))
    data=[]
    ll=0
    for ii in range(nz):
        for kk in range(ny):
            for jj in range(nx):
                data.append(YY.imag[jj,kk,ii])
    for iii in range((nx*ny*nz)) :
        if ((iii+1)%6==0) :	
            fo.write('  %1.5E\n'% (data[iii]))
        else :
            fo.write('  %1.5E'% (data[iii]))
    fo.write('\n END_DATAGRID_3D\n')
    fo.write('END_BLOCK_DATAGRID_3D')
    fo.close()


########################################################
def writeWavefunction(Y,fn=None):
    """
    Writes wave function to the specified output type
    Y: vector for wavefunction
                
    """
    if fn==None: fn = fileName()

    # Rotate in complex space
    max_amp=-1.0
    phase=1.0+0.0j

    for kk in range(len(Y)):
        if abs(Y[kk])>max_amp:
            max_amp=abs(Y[kk])
            phase=Y[kk]/max_amp
    Y=Y/phase

    foT=file(fn+'.abs.txt','w')
    foT.write('Atom nr M L abs(Y)\n')
    for ii in range(len(Y)):
        foT.write('%3.0i %3.0i %3.1i %3.1i %1.8f \n'%
        (basis.ii[ii],basis.atomnum[ii],basis.M[ii],basis.L[ii],abs(Y[ii])))
    foT.close()

    YY, dstep, origo, nx, ny, nz = calcWF(Y)

    # Write wave function in specified file format
    if general.format == 'macu':
        writemacubin(fn+'.Re.macu',YY.real,nx,ny,nz,origo,dstep)
        writemacubin(fn+'.Im.macu',YY.imag,nx,ny,nz,origo,dstep)

    if general.format == 'cube':
        writecube(fn+'.Re.cube',YY.real,nx,ny,nz,origo,dstep)
        writecube(fn+'.Im.cube',YY.imag,nx,ny,nz,origo,dstep)
        
    if general.format == 'XSF':
        writeXSF(fn+'.XSF',YY,nx,ny,nz,origo,dstep)

    if general.format == 'nc':
        writenetcdf(fn+'.nc',YY,nx,ny,nz,origo,dstep)


###########################################################
def setupParameters(calledFromInelastica=False):
    # Note: can be called from Inelastica as well!
    global general
    
    usage = "usage: %prog [options] DestinationDirectory"
    if calledFromInelastica:
        description = """
Inelastica script that calculates the IETS signal using the lowest order expansion, 
see Frederiksen et al. PRB 75, 205413 (2007) 

For help use --help!
"""
    else:
        description = """
Eigenchannels script that calculates:
1) Eigenchannels, Paulsson et al. PRB 76, 115117 (2007).
2) Fermi golden rule scattering rates, Paulsson et al. PRL 100, 226604 (2008).
3) Bond currents.

For help use --help!
"""
    parser = OptionParser(usage,description=description)
    EC = OptionGroup(parser, "Options for Eigenchannels")
    ECandIN = OptionGroup(parser, "Options for Eigenchannels and Inelastica")
    EC.add_option("-n", "--NumChan", dest="numchan", help="Number of eigenchannels [%default]", 
                  type='int', default=4)
    EC.add_option("-B", "--BothSides", dest='bothsides', default=False,action='store_true',
                  help="Calculate eigenchannels from both sides [%default]")
    EC.add_option("-M", "--MolecularStates", dest='MolStates', default=0.0, type='float',
                  help="Calculate eigenstates ov device Hamiltonian within [%default] eV from Ef")
    EC.add_option("-r", "--Res", dest='res', default=0.4,type='float',
                  help="Resolution [%default Ang]")
    EC.add_option("-w", "--format", dest='format', default='macu',type='string',
                  help="Wavefunction format (macu, cube, XSF, or nc) [%default]")
    EC.add_option("-N", "--NumPhCurr", dest='NumPhCurr', default=10,type='int',
                  help="Max number of changes in bond currents from inelastic scattering  [%default]")
    parser.add_option_group(EC)
    
    ECandIN.add_option("-e", "--Ef", dest='Ef', default=0.0,type='float',
                       help="Fermi energy [%default eV]")
    ECandIN.add_option("-f", "--fdf", dest='fdfFile', default='./RUN.fdf',type='string',
                       help="fdf file used for transiesta calculation [%default]")
    ECandIN.add_option("-p", "--PhononNetCDF", dest='PhononNetCDF', default=None,type='string',
                       help="Electron-phonon coupling NetCDF [%default]")
    ECandIN.add_option("-s", "--iSpin", dest='iSpin', default=1,type='int',
                       help="Spin channel [%default]")
    ECandIN.add_option("-k", "--kPoint", dest='kPoint', default='[0,0]',type='string',
                       help="2D k-point in range [0,1] given as string! default='[0.0,0.0]'")
    parser.add_option_group(ECandIN)

    if calledFromInelastica:
        IN = OptionGroup(parser, "Options for Inelastica")
        IN.add_option("-t", "--Temp", dest='Temp', default=4.2,type='float',
                      help="Temperature [%default K]")
        IN.add_option("-b", "--BiasPoints", dest='biasPoints', default=801,type='int',
                      help="Number of bias points [%default]")
        IN.add_option("-v", "--MinMaxVoltage", dest='MinMaxVoltage', default='-0.4:0.4',type='string',
                      help="Voltage range ['%default' V]")
        IN.add_option("-c", "--ModeCutoff", dest='modeCutoff', default='0.0025',type='float',
                      help="Ignore phonon modes with lower hw [%default eV]")
        IN.add_option("-V", "--Vrms", dest='Vrms', default='0.005',type='float',
                      help="Lock in amplifier broadening [%default V]")
        IN.add_option("-H", "--Heating", dest='PhHeating', default=False,action='store_true',
                      help="Include heating of vibrational modes [%default]")
        IN.add_option("-x", "--PhExtDamp", dest='PhExtDamp', default=1e-15,type='float',
                      help="External damping [%default (?) TODO check unit!]")
        parser.add_option_group(IN)
    
    (general, args) = parser.parse_args()
    print description

    general.energy, general.iSpin = general.Ef, general.iSpin-1
    
    if general.kPoint[0]!='[' or general.kPoint[-1]!=']':
        parser.error("ERROR: --kPoint='[0.0,0.0]' not --kPoint=%s"%general.kPoint)
    else:
        try:
            tmp=string.split(general.kPoint[1:-1],',')
            general.kPoint = N.array([float(tmp[0]),float(tmp[1]),0],N.float)
        except:
            parser.error("ERROR: --kPoint='[0.0,0.0]' not --kPoint=%s"%general.kPoint)

    if not os.path.exists(general.fdfFile):
        parser.error("No input fdf file found, specify with --fdf=file.fdf (default RUN.fdf)")

    try:
        general.from_atom = SIO.GetFDFlineWithDefault(
            general.fdfFile,'TS.TBT.PDOSFrom', int, None, 'Eigenchannels')
        general.to_atom = SIO.GetFDFlineWithDefault(
            general.fdfFile,'TS.TBT.PDOSTo', int, None, 'Eigenchannels')
    except:
        parser.error("Specify device region with TS.TBT.PDOS[To/From] keyword.")
    if calledFromInelastica:
        if general.PhononNetCDF==None:
            parser.error("ERROR: Inelastica needs a PhononNetCDF file!")
        try:
            tmp=string.split(general.MinMaxVoltage,':')
            general.minBias = float(tmp[0])
            general.maxBias = float(tmp[1])
        except:
            parser.error("ERROR: Inelastica failed to parse bias voltage!")

    print args
    if len(args)!=1:
        parser.error('ERROR: You need to specify destination directory')
    general.DestDir = args[0]
    if not os.path.isdir(general.DestDir):
        print '\nEigenchannels : Creating folder %s' %general.DestDir
        os.mkdir(general.DestDir)
    else:
        parser.error('ERROR: destination directory %s already exist!'%general.DestDir)

    def myprint(arg,file):
        # Save in parameter file
        print arg
        file.write(arg+'\n')

    class myopen:
        # Double stdout to RUN.out and stdout
        def write(self,x):
            self.stdout.write(x)
            self.file.write(x)

    fo = myopen()
    fo.stdout, fo.file = sys.stdout, open(general.DestDir+'/RUN.out','w',0)
    sys.stdout = fo

    file = open(general.DestDir+'/Parameters','w')    
    argv=""
    for ii in sys.argv: argv+=" "+ii
    myprint(argv,file)
    myprint('##################################################################################',file)
    myprint('## Eigenchannel options',file)
    myprint('fdfFile          : %s'%general.fdfFile,file)
    myprint('Ef [eV]          : %f'%general.energy,file)
    myprint('NumChan          : %i'%general.numchan,file)
    myprint('Res [Ang]        : %f'%general.res,file)
    myprint('PhononNetCDF     : %s'%general.PhononNetCDF,file)
    myprint('iSpin            : %i'%(general.iSpin+1),file)
    myprint('kPoint           : [%f,%f]'%(general.kPoint[0],general.kPoint[1]),file)
    myprint('BothSides        : %s'%general.bothsides,file)
    myprint('Device [from,to] : [%i,%i]'%(general.from_atom, general.to_atom),file)
    myprint('DestDir          : %s'%general.DestDir,file)
    myprint('##################################################################################',file)
    if calledFromInelastica:
        myprint('## Inelastica options',file)
        myprint('PhononNetCDF     : %s'%general.PhononNetCDF,file)
        myprint('Ef [eV]          : %f'%general.energy,file)
        myprint('Vrms [V]         : %f'%general.Vrms,file)
        myprint('PhHeating        : %s'%general.PhHeating,file)
        myprint('External damp [?]: %f'%general.PhExtDamp,file)
        myprint('mode cutoff [eV] : %f'%general.modeCutoff,file)
        myprint('Temperature [K]  : %f'%general.Temp,file)
        myprint('Bias points      : %i'%general.biasPoints,file)
        myprint('Min bias [V]     : %f'%general.minBias,file)
        myprint('Max bias [V]     : %f'%general.maxBias,file)
        myprint('##################################################################################',file)
    file.close()

################# File names ##################################
def fileName():
    # Generate filename '.EC.{1,Tot}{L,R,,In,Out}[UP][Ef=1.0].'
    if general.iChan==0:
        fn=general.SiestaOutFile+\
            '.EC.Tot%s%s'%(['L','R','','In','Out'][general.iSide],\
                               ['','UP','DOWN'][general.iSpin+HS.GF.HS.nspin-1])
    else:
        fn=general.SiestaOutFile+\
            '.EC.%i%s%s'%(general.iChan,['L','R','','In','Out'][general.iSide],\
                              ['','UP','DOWN'][general.iSpin+HS.GF.HS.nspin-1])
    if general.energy!=0.0:
        fn=fn+'Ef=%2.3f'%general.energy
    return general.DestDir+'/'+fn

################# Math helpers ################################
def mm(* args):
    # mm with arbitrary number of arguments
    tmp=args[0].copy()
    for mat in args[1:]:
        tmp=N.dot(tmp,mat)
    return tmp

def outerAdd(* args):
    # A_ijk=B_i+C_j+D_k
    tmp=args[0].copy()
    for ii in range(1,len(args)):
        tmp=N.add.outer(tmp,args[ii])
    return tmp

def dist(x):
    return N.sqrt(N.dot(x,x))

def mysqrt(x):
    # Square root of matrix    
    ev,U = LA.eig(x)
    U = N.transpose(U)

    tmp=N.zeros((len(ev),len(ev)),N.complex)
    for ii in range(len(ev)):
        tmp[ii,ii]=N.sqrt(ev[ii])
        
    return mm(LA.inv(U),tmp,U)

def dagger(x):
    return N.conjugate(N.transpose(x))

def sphericalHarmonic(ii):
    pi=3.141592654

    def Y00(sinth,costh,sinfi,cosfi):
        # l=0 m=0
        return 1/(2.*N.sqrt(pi))
    def Y1m1(sinth,costh,sinfi,cosfi):
        # l=1 m=-1
        return -(N.sqrt(3/pi)*sinfi*sinth)/2.
    def Y10(sinth,costh,sinfi,cosfi):
        # l=1 m=0
        return (costh*N.sqrt(3/pi))/2.
    def Y11(sinth,costh,sinfi,cosfi):
        # l=1 m=1
        return -(cosfi*N.sqrt(3/pi)*sinth)/2.
    def Y2m2(sinth,costh,sinfi,cosfi):
        # l=2 m=-2
        return (cosfi*N.sqrt(15/pi)*sinfi)/4. - \
               (cosfi*costh**2*N.sqrt(15/pi)*sinfi)/4. + \
               (cosfi*N.sqrt(15/pi)*sinfi*sinth**2)/4.
    def Y2m1(sinth,costh,sinfi,cosfi):
        # l=2 m=-1
        return -(costh*N.sqrt(15/pi)*sinfi*sinth)/2.
    def Y20(sinth,costh,sinfi,cosfi):
        # l=2 m=0
        return N.sqrt(5/pi)/8. + (3*costh**2*N.sqrt(5/pi))/8. - (3*N.sqrt(5/pi)*sinth**2)/8.
    def Y21(sinth,costh,sinfi,cosfi):
        # l=2 m=1
        return -(cosfi*costh*N.sqrt(15/pi)*sinth)/2.
    def Y22(sinth,costh,sinfi,cosfi):
        # l=2 m=2
        return (cosfi**2*N.sqrt(15/pi))/8. - (cosfi**2*costh**2*N.sqrt(15/pi))/    8. - (N.sqrt(15/pi)*sinfi**2)/8. + (costh**2*N.sqrt(15/pi)*sinfi**2)/    8. + (cosfi**2*N.sqrt(15/pi)*sinth**2)/8. - (N.sqrt(15/pi)*sinfi**2*sinth**2)/    8.
    def Y3m3(sinth,costh,sinfi,cosfi):
        # l=3 m=-3
        return (-9*cosfi**2*N.sqrt(35/(2.*pi))*sinfi*sinth)/    16. + (9*cosfi**2*costh**2*N.sqrt(35/(2.*pi))*sinfi*sinth)/    16. + (3*N.sqrt(35/(2.*pi))*sinfi**3*sinth)/    16. - (3*costh**2*N.sqrt(35/(2.*pi))*sinfi**3*sinth)/    16. - (3*cosfi**2*N.sqrt(35/(2.*pi))*sinfi*sinth**3)/    16. + (N.sqrt(35/(2.*pi))*sinfi**3*sinth**3)/16.
    def Y3m2(sinth,costh,sinfi,cosfi):
        # l=3 m=-2
        return (cosfi*costh*N.sqrt(105/pi)*sinfi)/8. - (cosfi*costh**3*N.sqrt(105/pi)*sinfi)/    8. + (3*cosfi*costh*N.sqrt(105/pi)*sinfi*sinth**2)/8.
    def Y3m1(sinth,costh,sinfi,cosfi):
        # l=3 m=-1
        return -(N.sqrt(21/(2.*pi))*sinfi*sinth)/    16. - (15*costh**2*N.sqrt(21/(2.*pi))*sinfi*sinth)/    16. + (5*N.sqrt(21/(2.*pi))*sinfi*sinth**3)/16.
    def Y30(sinth,costh,sinfi,cosfi):
        # l=3 m=0
        return (3*costh*N.sqrt(7/pi))/16. + (5*costh**3*N.sqrt(7/pi))/    16. - (15*costh*N.sqrt(7/pi)*sinth**2)/16.
    def Y31(sinth,costh,sinfi,cosfi):
        # l=3 m=1
        return -(cosfi*N.sqrt(21/(2.*pi))*sinth)/    16. - (15*cosfi*costh**2*N.sqrt(21/(2.*pi))*sinth)/    16. + (5*cosfi*N.sqrt(21/(2.*pi))*sinth**3)/16.
    def Y32(sinth,costh,sinfi,cosfi):
        # l=3 m=2
        return (cosfi**2*costh*N.sqrt(105/pi))/16. - (cosfi**2*costh**3*N.sqrt(105/pi))/    16. - (costh*N.sqrt(105/pi)*sinfi**2)/    16. + (costh**3*N.sqrt(105/pi)*sinfi**2)/    16. + (3*cosfi**2*costh*N.sqrt(105/pi)*sinth**2)/    16. - (3*costh*N.sqrt(105/pi)*sinfi**2*sinth**2)/16.
    def Y33(sinth,costh,sinfi,cosfi):
        # l=3 m=3
        return (-3*cosfi**3*N.sqrt(35/(2.*pi))*sinth)/    16. + (3*cosfi**3*costh**2*N.sqrt(35/(2.*pi))*sinth)/    16. + (9*cosfi*N.sqrt(35/(2.*pi))*sinfi**2*sinth)/    16. - (9*cosfi*costh**2*N.sqrt(35/(2.*pi))*sinfi**2*sinth)/    16. - (cosfi**3*N.sqrt(35/(2.*pi))*sinth**3)/    16. + (3*cosfi*N.sqrt(35/(2.*pi))*sinfi**2*sinth**3)/16.
      
    l=basis.L[ii]
    m=basis.M[ii]

    tmp=[[Y00],[Y1m1,Y10,Y11],[Y2m2,Y2m1,Y20,Y21,Y22],[Y3m3,Y3m2,Y3m1,Y30,Y31,Y32,Y33]]

    if l==3:
        print 'f-shell : l=%i, m=%i (NOT TESTED!!)'%(l,m)

    return tmp[l][m+l]

    
##################### Start main routine #####################

if __name__ == '__main__':
    main()

