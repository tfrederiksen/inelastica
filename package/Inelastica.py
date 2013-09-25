print "SVN $Id$"

"""
Inelastica:
1: Analyze IETS spectra, Paulsson et al PRL 2008
"""

import EigenChannels as EC
import NEGF
import SiestaIO as SIO
import MakeGeom as MG
import MiscMath as MM
import WriteNetCDF as NCDF
import numpy as N
import numpy.linalg as LA
import Scientific.IO.NetCDF as NC
import sys
import PhysicalConstants as PC
import time

########################################################
##################### Main routine #####################
########################################################
def main(options):
    options.XV = '%s/%s.XV'%(options.head,options.systemlabel)
    options.geom = MG.Geom(options.XV)
    # Set up electrodes and device Greens function
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL,options.NA1L,options.NA2L,options.voltage/2.)
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR,options.NA1R,options.NA2R,-options.voltage/2.)
    myGF = NEGF.GF(options.TSHS,elecL,elecR,Bulk=True,DeviceAtoms=[options.devSt, options.devEnd])
    # Read phonons
    NCfile = NC.NetCDFFile(options.PhononNetCDF,'r')
    print 'Inelastica: Reading ',options.PhononNetCDF
    hw = N.array(NCfile.variables['hw'][:])
    # Prepare lists for various trace factors
    #myGF.dGnout = []
    #myGF.dGnin = []
    myGF.P1T = []     # M.A.M.A (total e-h damping)
    myGF.P2T = []     # M.AL.M.AR (emission)
    myGF.ehDampL = [] # M.AL.M.AL (L e-h damping)
    myGF.ehDampR = [] # M.AR.M.AR (R e-h damping)
    myGF.nHT = []     # non-Hilbert/Isym factor
    myGF.nHTin = []
    myGF.nHTel = []
    myGF.HT = []      # Hilbert/Iasym factor
    # Calculate transmission at Fermi level
    myGF.calcGF(options.energy+options.eta*1.0j,options.kPoint[0:2],ispin=options.iSpin,etaLead=options.etaLead,useSigNCfiles=options.signc)
    basis = SIO.BuildBasis(options.XV,options.devSt,options.devEnd,myGF.HS.lasto)
    myGF.TeF = myGF.calcT(options.numchan)[0]
    # Check consistency of PHrun vs TSrun inputs
    IntegrityCheck(options,myGF,basis,NCfile)   
    # Calculate trace factors one mode at a time
    for ihw in range(len(hw)):
        calcTraces(options,myGF,basis,NCfile,ihw)
    # Multiply traces with voltage-dependent functions
    calcIETS(options,myGF,basis,hw)
    NCfile.close()


########################################################
def IntegrityCheck(options,myGF,basis,NCfile):
    # Perform consistency checks for device region in 
    # PH and TS calculations by comparing coordinates
    # and atom numbers
    PH_dev = N.array(NCfile.variables['DeviceAtoms'][:])
    PH_xyz = N.array(NCfile.variables['GeometryXYZ'][:])
    PH_anr = N.array(NCfile.variables['AtomNumbers'][:])
    TS_dev = range(options.devSt,options.devEnd+1)
    TS_anr = options.geom.anr
    TS_xyz = options.geom.xyz
    print '\nInelastica.IntegrityCheck:'
    print 'A = %s'%options.PhononNetCDF
    print 'B = %s'%options.XV
    print ' idxA    xA       yA       zA   anrA   ',
    print ' idxB    xB       yB       zB   anrB'
    for i in range(max(len(PH_dev),len(TS_dev))):
        # Geom A
        if PH_dev[0]+i in PH_dev:
            s = ('%i'%(PH_dev[0]+i)).rjust(5)
            for j in range(3):
                s += ('%.4f'%(PH_xyz[PH_dev[0]-1+i,j])).rjust(9)
            s += ('%i'%PH_anr[PH_dev[0]-1+i]).rjust(4)
        else:
            s = ('---').center(36)
        s += '  vs'
        # Geom B
        if options.devSt+i in TS_dev:
            s += ('%i'%(options.devSt+i)).rjust(5)
            for j in range(3):
                s += ('%.4f'%(TS_xyz[options.devSt-1+i,j])).rjust(9)
            s += ('%i'%TS_anr[options.devSt-1+i]).rjust(4)
        else:
            s += ('---').center(36)
        print s
    # - check 1: Matrix sizes
    PH_H0 = N.array(NCfile.variables['H0'][:])
    if N.shape(PH_H0[0])==N.shape(myGF.Gr):
        print '... Check 1 passed: Device orb. space matches'
        check1 = True
    else:
        print '... Check 1 failed: Device orb. space do not match!!!'
        check1 = False
    # - check 2&3: Geometry and atom number
    dist_xyz = 0.0
    dist_anr = 0.0
    for i in range(len(PH_dev)):
        # Geometric distance between atoms
        if i==0:
            # Allow for a global offset of coordinates R
            R = PH_xyz[PH_dev[0]-1+i]-TS_xyz[options.devSt-1+i]
            print 'Global offset R = [%.3f %.3f %.3f]'%(R[0],R[1],R[2])
        d = PH_xyz[PH_dev[0]-1+i]-TS_xyz[options.devSt-1+i] - R
        dist_xyz += N.dot(d,d)**.5
        # Difference between atom numbers
        a = PH_anr[PH_dev[0]-1+i]-TS_anr[options.devSt-1+i]
        dist_anr += abs(a)
    if dist_xyz<1e-3:
        print '... Check 2 passed: Atomic coordinates consistent'
        check2 = True
    elif dist_xyz<0.1:
        print '... Check 2 WARNING: Atomic coordinates deviate by %.3f Ang!!!'%dist_xyz
        check2 = True
    else:
        print '... Check 2 failed: Atomic coordinates deviate by %.3f Ang!!!'%dist_xyz
        check2 = False
    if dist_anr<1e-3:
        print '... Check 3 passed: Atomic numbers consistent'
        check3 = True
    else:
        print '... Check 3 failed: Atomic numbers inconsistent!!!'
        check3 = False
    if (not check1) or (not check2) or (not check3):
        sys.exit('Inelastica: Error - inconsistency detected for device region.\n')


def calcTraces(options,myGF,basis,NCfile,ihw):
    # Calculate various traces over the electronic structure
    G = myGF.Gr
    Gd = MM.dagger(G)
    nuo, nuoL, nuoR = myGF.nuo, myGF.nuoL, myGF.nuoR
    GL = myGF.GamL
    GR = myGF.GamR
    # Spectral functions
    AL = MM.mm(G[:,0:nuoL],GL,Gd[0:nuoL,:])
    AR = MM.mm(G[:,nuo-nuoR:nuo],GR,Gd[nuo-nuoR:nuo,:])
    A = AL + AR
    ALT = MM.mm(Gd[:,0:nuoL],GL,G[0:nuoL,:])
    # Electron-phonon couplings
    M = N.array(NCfile.variables['He_ph'][ihw,options.iSpin,:,:])
    # Calculation of intermediate quantity
    MARGLGM = MM.mm(M,AR[:,0:nuoL],GL,G[0:nuoL,:],M)
    # LOE expressions in compact form
    t1 = MM.mm(MARGLGM,AR)
    t2 = MM.mm(MARGLGM,AL)
    K23 = N.trace(t1+t2).imag
    K4 = N.trace(MM.mm(M,ALT,M,AR))
    aK23 = 2*N.trace(t1-t2).real # asymmetric part
    # Non-Hilbert term defined here with a minus sign
    myGF.nHTin.append(-1.0*checkImPart(K4))
    myGF.nHTel.append(-1.0*checkImPart(K23))
    myGF.nHT.append(-1.0*checkImPart(K23+K4))
    myGF.HT.append(checkImPart(aK23))
    # Power, damping and current rates
    myGF.P1T.append(checkImPart(N.trace(MM.mm(M,A,M,A))))
    myGF.P2T.append(checkImPart(N.trace(MM.mm(M,AL,M,AR))))
    myGF.ehDampL.append(checkImPart(N.trace(MM.mm(M,AL,M,AL))))
    myGF.ehDampR.append(checkImPart(N.trace(MM.mm(M,AR,M,AR))))
    # Remains from older version (see before rev. 219):
    #myGF.dGnout.append(EC.calcCurrent(options,basis,myGF.HNO,mm(Us,-0.5j*(tmp1-dagger(tmp1)),Us)))
    #myGF.dGnin.append(EC.calcCurrent(options,basis,myGF.HNO,mm(Us,mm(G,MA1M,Gd)-0.5j*(tmp2-dagger(tmp2)),Us)))
    # NB: TF Should one use myGF.HNO (nonorthogonal) or myGF.H (orthogonalized) above?

def calcIETS(options,myGF,basis,hw):
    # Calculate product of electronic traces and voltage functions
    print 'Inelastica.calcIETS: nHT =',N.array(myGF.nHT)*PC.unitConv # OK
    print 'Inelastica.calcIETS: HT  =',N.array(myGF.HT) # OK
    
    # Set up grid and Hilbert term
    kT = options.Temp/11604.0 # (eV)

    # Generate grid for numerical integration of Hilbert term    
    max_hw=max(hw)
    max_win=max(-options.minBias,max_hw)+20*kT+4*options.Vrms
    min_win=min(-options.maxBias,-max_hw)-20*kT-4*options.Vrms
    pts=int(N.floor((max_win-min_win)/kT*3))
    Egrid=N.array(range(pts),N.float)/pts*(max_win-min_win)+min_win
    print "Inelastica.calcIETS: Hilbert integration grid : %i pts [%f,%f]" % (pts,min(Egrid),max(Egrid))
    
    # Calculate the prefactors for the Hilbert and non-Hilbert terms
    # Also calculate the Hilbert transfer of the box on the grid for each mode
    hilb=[] # List of hilbert transforms
    ker=None
    for ii in range(len(hw)):
        hwi=hw[ii]
        tmp=MM.box(hwi,-hwi,Egrid,kT)
        tmp2, ker = MM.Hilbert(tmp,ker)
        hilb.append(tmp2)
        SIO.printDone(ii,len(hw),'Inelastica.calcIETS: Hilbert transform')
    
    NN = options.biasPoints
    print 'Inelastica.calcIETS: Biaspoints =',NN

    # Add some points for the Lock in broadening
    approxdV=(options.maxBias-options.minBias)/NN
    NN+=int(((8*options.Vrms)/approxdV)+.5)
    
    Vl=options.minBias-4*options.Vrms+ \
        (options.maxBias-options.minBias+8*options.Vrms)/NN*N.array(range(NN),N.float)
    
    InH=N.zeros((NN,),N.complex) # Non-Hilb-term current
    IH=N.zeros((NN,),N.complex) # Hilb-term current
    Pow=N.zeros((NN,),N.float) # Power
    nPhtot=N.zeros((NN,),N.float) # Number of phonons (sum over modes)
    nPhvsBias=N.zeros((NN,len(hw)),N.float) # Number of phonons
    nPh=N.zeros((len(hw),),N.float)         # Number of phonons                                  

    for iV in range(len(Vl)):
        SIO.printDone(iV,len(Vl),'Inelastica.calcIETS: Calculating voltage-dependent functions')
        
        InH[iV], IH[iV], Pow[iV] = 0.0, 0.0, 0.0
        V=Vl[iV]
        
        kasse=MM.box(0,-V,Egrid,kT)

        for iOmega in range(len(hw)):
            ihw = hw[iOmega]
            if ihw>options.modeCutoff:
                PV=abs(V)
                
                # Power
                exphw=N.exp(N.clip(ihw/kT,-300,300))
                expV=N.exp(N.clip(PV/kT,-300,300))
                exphwmV=N.exp(N.clip((ihw-PV)/kT,-300,300))
                if ihw/kT>300:
                    mexphw=0.0
                else:
                    mexphw=N.exp(N.clip(-ihw/kT,-300,300))
                if PV/kT>300:
                    mexpV=0.0
                else:
                    mexpV=N.exp(N.clip(-PV/kT,-300,300))

                # Re-written third time lucky? Numerical problems galore!
                if (ihw-PV)/kT>150:
                    tmpheat=0
                else:
                    if abs((ihw-PV)/kT)<1e-3:
                        tmpheat= (2*mexphw*ihw+kT*mexphw**2-kT)/(mexphw**2-1)-\
                                 (PV-ihw)*(4*ihw*mexphw**2+kT*mexphw**4-kT)/(2*kT*(1-mexphw**2)**2)
                    else:
                        tmpheat=(1-mexpV)/(1-mexphw)* \
                                 ((1+mexphw)*(1-mexpV)*ihw-(1-mexphw)*(1+mexpV)*PV)/ \
                                 (exphwmV+mexphw*mexpV-mexpV**2-1)
                        
                heat=tmpheat*myGF.P2T[iOmega]/N.pi # Heating term / (hbar hw)
                
                if options.PhHeating: # Heating?
                    nPh[iOmega]=heat/(ihw*myGF.P1T[iOmega]/N.pi+options.PhExtDamp)+1.0/(exphw-1.0)
                nPhvsBias[iV,iOmega]=nPh[iOmega]
                    
                # Damping term /(hbar hw)
                damp = ihw*(1/(exphw-1)-nPh[iOmega])*myGF.P1T[iOmega]/N.pi
                
                # Power in units of 1/(hbar)
                Pow[iV]=ihw*(heat+damp)
                nPhtot[iV]=nPhtot[iV]+nPh[iOmega]
                tmp=0.0
                if abs(V-ihw)/kT<1e-7:
                    tmp=-kT
                else:
                    tmp=(V-ihw)/(N.exp(N.clip((ihw-V)/kT,-70.0,70.0))-1)
                if abs(V+ihw)/kT<1e-7:
                    tmp+=kT
                else:
                    tmp+=(V+ihw)/(N.exp(N.clip((ihw+V)/kT,-70.0,70.0))-1)
                    InH[iV]+=(tmp-V*nPh[iOmega])*myGF.nHT[iOmega]

                # Finite temp Hilbert
                IH[iV]-=myGF.HT[iOmega]*MM.trapez(Egrid,kasse*hilb[iOmega],\
                                                     equidistant=True)/2
                
        InH[iV]+=V*myGF.TeF # Last to reduce cancelation errors

    # Get the right units for gamma_eh, gamma_heat
    gamma_eh=N.zeros((len(hw),),N.float)
    gamma_heat=N.zeros((len(hw),),N.float)
    for iOmega in range(len(hw)):
        # Units [Phonons per Second per dN where dN is number extra phonons]
        gamma_eh[iOmega]=myGF.P1T[iOmega]*hw[iOmega]*PC.unitConv
        # Units [Phonons per second per eV [eV-ihw]
        gamma_heat[iOmega]=myGF.P2T[iOmega]*PC.unitConv

    print 'Inelastica.calcIETS: gamma_eh =',gamma_eh # OK
    print 'Inelastica.calcIETS: gamma_heat =',gamma_heat # OK

    #hw, T, nHT, HT, lLOE, nPhtot, nPh, = hw, myGF.TeF, myGF.nHT, myGF.HT, [Vl, InH, IH], nPhtot, nPhvsBias
    V, I, dI, ddI, BdI, BddI, NnPhtot,NnPh = Broaden(options,Vl,InH+IH,nPhtot,nPhvsBias)

    print 'Inelastica.calcIETS: BdI[:10] =',BdI[:10] # OK
    print 'Inelastica.calcIETS: BddI[:10] =',BddI[:10] # OK

    datafile = '%s/%s.IN'%(options.DestDir,options.systemlabel)
    initncfile(datafile,hw)
    writeLOEData2Datafile(datafile,hw,myGF.TeF,myGF.nHT,myGF.HT)
    writeLOE2ncfile(datafile,hw,myGF.nHT,myGF.HT,V,I,NnPhtot,NnPh,\
                    dI,ddI,BdI,BddI,gamma_eh,gamma_heat)

    
########################################################
def checkImPart(x):
    if abs(x.imag)>0.0000001:
        print "Inelastica.checkImPart: Imaginary part (%.3e) too big"%x.imag
        kuk
    return x.real   

########################################################
def writeFGRrates():
    # This does not work at the moment....
    NCfile = NC.NetCDFFile(options.PhononNetCDF,'r')
    print 'Reading ',options.PhononNetCDF

    options.iChan, options.iSide = 0, 2
    outFile = file('%s/%s.IN.FGR'%(options.DestDir,options.systemlabel,'w'))
    outFile.write('Total transmission [in units of (1/s/eV)] : %e\n' % (PC.unitConv*myGF.totTrans.real,))

    tmp=N.sort(abs(N.array(myGF.nHT[:])))
    SelectionMin=tmp[-options.NumPhCurr]        
    
    for ihw in range(len(myGF.hw)):
        SIO.printDone(ihw,len(myGF.hw),'Golden Rate') 
        M = N.array(NCfile.variables['He_ph'][ihw,options.iSpin,:,:])
        rate=N.zeros((len(myGF.ECleft),len(myGF.ECright)),N.float)
        totrate=0.0
        inter,intra = 0.0, 0.0 # splitting total rate in two
        for iL in range(len(myGF.ECleft)):
            for iR in range(len(myGF.ECright)):
                tmp=N.dot(N.conjugate(myGF.ECleft[iL]),mm(M,myGF.ECright[iR]))
                rate[iL,iR]=(2*N.pi)**2*abs(tmp)**2
                totrate+=rate[iL,iR]
                if iL==iR: intra += rate[iL,iR]
                else: inter += rate[iL,iR]

        if abs(myGF.nHT[ihw])>=SelectionMin:
            options.iChan = ihw
            currOut, currIn = IETS.dGnout[ihw], IETS.dGnin[ihw]
            options.iSide = 3
            writeCurrent(currIn)
            options.iSide = 4
            writeCurrent(currOut)
            
        outFile.write('\nPhonon mode %i : %f eV [Rates in units of (1/s/eV)]\n' % (ihw,IETS.hw[ihw]))
        outFile.write('eh-damp : %e (1/s) , heating %e (1/(sV)))\n' % (IETS.P1T[ihw]*PC.unitConv*IETS.hw[ihw],IETS.P2T[ihw]*PC.unitConv))
        outFile.write('eh-damp 1, 2 (MALMAL, MARMAR): %e (1/s) , %e (1/(s)))\n' % (IETS.ehDampL[ihw]*PC.unitConv*IETS.hw[ihw],IETS.ehDampR[ihw]*PC.unitConv*IETS.hw[ihw]))
        outFile.write('SymI : %e (1/(sV)) , AsymI %e (?))\n' % (IETS.nHT[ihw]*PC.unitConv,IETS.HT[ihw]*PC.unitConv))
        outFile.write('Elast : %e (1/(sV)) , Inelast %e (1/(sV)))\n' % (IETS.nHTel[ihw]*PC.unitConv,IETS.nHTin[ihw]*PC.unitConv))
        outFile.write('down=left EC, right=right EC\n')
        if IETS.P2T[ihw]>0.0:
            if abs(totrate/(IETS.P2T[ihw])-1)<0.05:
                outFile.write('Sum/Tr[MALMAR] , Tr: %1.3f  %e\n'%(totrate/(IETS.P2T[ihw]),PC.unitConv*IETS.P2T[ihw]))
            else:
                outFile.write('WARNING: !!!! Sum/Tr[MALMAR] , Tr: %2.2e  %e\n'%(totrate/(IETS.P2T[ihw]),PC.unitConv*IETS.P2T[ihw]))
        else:
            outFile.write(' Tr:  %e\n'%(PC.unitConv*IETS.P2T[ihw]))
        inter = inter/IETS.P2T[ihw]
        intra = intra/IETS.P2T[ihw]
        outFile.write('Interchannel ratio: Sum(inter)/Tr[MALMAR]      = %.4f \n'%inter)
        outFile.write('Intrachannel ratio: Sum(intra)/Tr[MALMAR]      = %.4f \n'%intra)
        outFile.write('Inter+intra ratio: Sum(inter+intra)/Tr[MALMAR] = %.4f \n'%(inter+intra))
        for iL in range(len(IETS.ECleft)):
            for iR in range(len(IETS.ECright)):
                outFile.write('%e ' % (PC.unitConv*rate[iL,iR],))
            outFile.write('\n')
        outFile.close()
        NCfile.close()
                        


################################################################
# Broadening due to Vrms
################################################################

def Broaden(options,VV,II,nPhtot,nPh):
    """
    Broadening corresponding to Lock in measurements for the
    conductance and IETS spectra. Also resample II, nPh and nPhtot
    to match a common voltage list
    """

    II=II.copy()
    II=II.real
    
    # First derivative dI and bias list dV
    dI=(II[1:len(II)]-II[:-1])/(VV[1]-VV[0])
    dV=(VV[1:len(VV)]+VV[:-1])/2
    # Second derivative and bias ddV
    ddI=(dI[1:len(dI)]-dI[:-1])/(VV[1]-VV[0])
    ddV=(dV[1:len(dV)]+dV[:-1])/2

    # Modulation amplitude
    VA=N.sqrt(2.0)*options.Vrms 

    # New bias ranges for broadening
    tmp=int(N.floor(VA/(dV[1]-dV[0]))+1)
    BdV=dV[tmp:-tmp]
    BddV=ddV[tmp:-tmp]

    # Initiate derivatives
    BdI=0*BdV
    BddI=0*BddV
    
    # Calculate first derivative with Vrms broadening
    for iV, V in enumerate(BdV):
        SIO.printDone(iV,len(BdV),'Inelastica.Broaden: First-derivative Vrms broadening')
        wt=(N.array(range(200))/200.0-0.5)*N.pi
        VL=V+VA*N.sin(wt)
        dIL=MM.interpolate(VL,dV,dI)
        BdI[iV]=2/N.pi*N.sum(dIL*(N.cos(wt)**2))*(wt[1]-wt[0])

    # Calculate second derivative with Vrms broadening    
    for iV, V in enumerate(BddV):
        SIO.printDone(iV,len(BddV),'Inelastica.Broaden: Second-derivative Vrms broadening')
        wt=(N.array(range(200))/200.0-0.5)*N.pi
        VL=V+VA*N.sin(wt)
        ddIL=MM.interpolate(VL,ddV,ddI)
        BddI[iV]=8.0/3.0/N.pi*N.sum(ddIL*(N.cos(wt)**4))*(wt[1]-wt[0])

    # Reduce to one voltage grid
    NN=options.biasPoints 
    V=options.minBias+(options.maxBias-options.minBias)/NN*N.array(range(NN)) 

    NI=MM.interpolate(V,VV,II)
    NdI=MM.interpolate(V,dV,dI)
    NddI=MM.interpolate(V,ddV,ddI)
    NBdI=MM.interpolate(V,BdV,BdI)
    NBddI=MM.interpolate(V,BddV,BddI)
    NnPhtot=MM.interpolate(V,VV,nPhtot)
    NnPh=N.zeros((len(V),len(nPh[0,:])),N.float)
    for ii in range(len(nPh[0,:])):
        NnPh[:,ii]=MM.interpolate(V,VV,nPh[:,ii])
    
    return V, NI ,NdI, NddI, NBdI, NBddI, NnPhtot, NnPh

################################################################
# Output to NetCDF file
################################################################

def initncfile(filename,hw):
    'Initiate netCDF file'
    print 'Inelastica: Initializing nc-file'
    ncfile = NC.NetCDFFile(filename+'.nc','w','Created '+time.ctime(time.time()))
    ncfile.title = 'Inelastica Output'
    ncfile.version = 1
    ncfile.createDimension('Nph',len(hw))
    ncfile.createDimension('Bias',None)
    nchw = ncfile.createVariable('hw','d',('Nph',))
    nchw[:] = N.array(hw)
    nchw.units = 'eV'
    ncfile.close()
    
def writeLOE2ncfile(filename,hw,nHT,HT,V,I,nPhtot,nPh,\
                     dI,ddI,BdI,BddI,gamma_eh,gamma_heat):
    'Write LOE data to netCDF file'
    print 'Inelastica: Write LOE data to nc-file'
    ncfile = NC.NetCDFFile(filename+'.nc','a')
    ncfile.createDimension('LOE_V',len(V))
    tmp=ncfile.createVariable('LOE_V','d',('LOE_V',))
    tmp[:]=N.array(V)
    tmp.units='V'
    tmp=ncfile.createVariable('LOE_IETS','d',('LOE_V',))
    tmpn=BddI/BdI
    tmp[:]=N.array(tmpn.real)
    tmp.units='Broadened ddI/dI [1/V]'
    tmp=ncfile.createVariable('LOE_dI','d',('LOE_V',))
    tmp[:]=N.array(dI.real)
    tmp.units='dI LOE, G0'
    tmp=ncfile.createVariable('LOE_ddI','d',('LOE_V',))
    tmp[:]=N.array(ddI.real)
    tmp.units='ddI LOE, G0/V'
    tmp=ncfile.createVariable('LOE_BdI','d',('LOE_V',))
    tmp[:]=N.array(BdI.real)
    tmp.units='Broadened dI LOE, G0'
    tmp=ncfile.createVariable('LOE_BddI','d',('LOE_V',))
    tmp[:]=N.array(BddI.real)
    tmp.units='Broadened ddI LOE, G0/V'
    tmp=ncfile.createVariable('LOE_I','d',('LOE_V',))
    tmp[:]=N.array(I.real)
    tmp.units='G0 V'
    tmp=ncfile.createVariable('LOE_tot_nPh','d',('LOE_V',))
    tmp[:]=N.array(nPhtot.real)
    tmp.units='Total number of phonons '
    tmp=ncfile.createVariable('LOE_nPh','d',('LOE_V','Nph'))
    tmp[:]=N.array(nPh.real)
    tmp.units='Number of phonons'
    tmp=ncfile.createVariable('LOE_gamma_eh','d',('Nph',))
    tmp[:]=N.array(gamma_eh)
    tmp.units='e-h damping [*deltaN=1/Second]'
    tmp=ncfile.createVariable('LOE_gamma_heat','d',('Nph',))
    tmp[:]=N.array(gamma_heat)
    tmp.units='Phonon heating [*(bias-hw) (eV) = 1/Second]'
    ncTotR = ncfile.createVariable('LOE_ISymTr','d',('Nph',))
    ncTotR[:] = N.array(nHT)*PC.unitConv
    ncTotR.units = 'Trace giving Symmetric current contribution (same units as eigenchannels [1/s/eV])'
    ncTotL = ncfile.createVariable('LOE_IAsymTr','d',('Nph',))
    ncTotL[:] = N.array(HT)
    ncTotL.units = 'Trace giving Asymmetric current contribution'
    ncfile.close()
    
def writeLOEData2Datafile(file,hw,T,nHT,HT):
    f = open(file,'a')
    f.write("## Almost First Order Born Approximation (kbT=0)\n")
    f.write("## hw(eV)      Trans        non-H-term (e/pi/hbar)   H-term\n")
    for ii in range(len(hw)):
        f.write("## %e %e %e %e\n" % (hw[ii],T,nHT[ii],HT[ii]))
    f.write("\n")
    f.close()

