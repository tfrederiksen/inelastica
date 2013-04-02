print "SVN $Id$"

"""
Eigenchannels:
1: Eigenchannels, method from Paulsson and Brandbyge PRB 2007
2: Analyze IETS spectra, Paulsson et al PRL 2008
3: Calculate "bond" currents
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
import sys, string, struct, glob,  os
from optparse import OptionParser, OptionGroup
import PhysicalConstants as PC
import time

########################################################
##################### Main routine #####################
########################################################
def main():
    general = setupParameters()
    geom = readxv()
    myGF = readHS(general)
    basis = readbasis(general,myGF.HS)
    calcTeig(general,myGF)
    calcIETS(general,myGF,basis)

########################################################
def readxv():
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
    return geom

########################################################
def readbasis(param,GF):
    fn=glob.glob('*.XV')
    basis = SIO.BuildBasis(fn[0], param.from_atom, param.to_atom, GF.lasto)
    return basis

########################################################
def readHS(general):
    # Setup H, S and self-energies
    fnL  = SIO.GetFDFlineWithDefault(general.fdfFile,'TS.HSFileLeft', str, None, 'Error Eigenchannels')
    NA1L = SIO.GetFDFlineWithDefault(general.fdfFile,'TS.ReplicateA1Left', int, 1, 'Error Eigenchannels')
    NA2L = SIO.GetFDFlineWithDefault(general.fdfFile,'TS.ReplicateA2Left', int, 1, 'Error Eigenchannels')
    fnR  = SIO.GetFDFlineWithDefault(general.fdfFile,'TS.HSFileRight', str, None, 'Error Eigenchannels')
    NA1R = SIO.GetFDFlineWithDefault(general.fdfFile,'TS.ReplicateA1Right', int, 1, 'Error Eigenchannels')
    NA2R = SIO.GetFDFlineWithDefault(general.fdfFile,'TS.ReplicateA2Right', int, 1, 'Error Eigenchannels')
    voltage  =SIO.GetFDFlineWithDefault(general.fdfFile,'TS.Voltage', float, 0.0, 'Error Eigenchannels')

    elecL = NEGF.ElectrodeSelfEnergy(fnL,NA1L,NA2L,voltage/2.)
    elecR = NEGF.ElectrodeSelfEnergy(fnR,NA1R,NA2R,-voltage/2.)

    general.systemlabel = SIO.GetFDFlineWithDefault(general.fdfFile,'SystemLabel', str, None, 'pyTBT')

    myGF = NEGF.GF(general.systemlabel+'.TSHS',elecL,elecR,Bulk=True,DeviceAtoms=[general.from_atom, general.to_atom])

    myGF.calcGF(general.energy+general.eta*1j, general.kPoint[0:2], ispin=general.iSpin)

    return myGF

########################################################
def calcTeig(general,myGF):
    # Matrix to save total and eigenchannel transmissions
    # BEFORE ORTHOGO
    T = myGF.calcT(general.numchan)
    print 'Transmission T(E=%.4f) [Ttot, T1, T2, ... Tn]:'%general.energy
    for t in T:
        print '%.9f '%t,
    print

    # Now orthogonalize in device region
    #HNO = myGF.H.copy() # nonorthogonal device Hamiltonian (needed later)
    myGF.orthogonalize()

    # Check that we get the same transmission:
    T = myGF.calcT(general.numchan)
    print 'Transmission T(E=%.4f) [Ttot, T1, T2, ... Tn]:'%general.energy
    for t in T:
        print '%.9f '%t,
    print

    # Calculate Eigenchannels
    myGF.calcEigChan(general.numchan)
    print 'Left eigenchannel transmissions:',myGF.EigTleft
    print 'Right eigenchannel transmissions:',myGF.EigTright
    myGF.trans0 = T[0]

########################################################
def calcIETS(general,myGF,basis):
    # Calculate inelastic scattering rates, total and per eigenchannel

    NCfile = NC.NetCDFFile(general.PhononNetCDF,'r')
    print 'Reading ',general.PhononNetCDF
    hw = N.array(NCfile.variables['hw'][:]) 
    # Write matrices to file?
    WriteOrtho = False
    if WriteOrtho:
        ncdf = NCDF.NCfile(general.DestDir+'/ortho.nc')
        ncdf.write(hw,'hw')
    G = myGF.Gr
    G1 = myGF.GamL
    G2 = myGF.GamR
    Gd = dagger(G)
    A1 = myGF.A1
    A2 = myGF.A2
    Us = myGF.Us
    
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
        
        # Write matrices to file?
        if WriteOrtho:
            ncdf.write(H1.real,'ReM%.3i'%ihw)
            ncdf.write(H1.imag,'ImM%.3i'%ihw)
        
        # Cryptic? Changes in G^lesser due to e-ph interaction at high and low energies, i.e.,
        # the changes in occupation due to out(in)-scattering
        tmp1, tmp2 = mm(G,MA2M,A1), mm(G,MA1M,A1)
        dGnout.append(EC.calcCurrent(general,basis,myGF.HNO,mm(Us,-0.5j*(tmp1-dagger(tmp1)),Us)))
        dGnin.append(EC.calcCurrent(general,basis,myGF.HNO,mm(Us,mm(G,MA1M,Gd)-0.5j*(tmp2-dagger(tmp2)),Us)))
        # NB: TF Should one use myGF.HNO or myGF.H above?

        # Power, damping and current rates
        P1T.append(checkImPart(N.trace(mm(MAM,A1+A2))))
        P2T.append(checkImPart(N.trace(mm(MA1M,A2))))        
        ehDamp1.append(checkImPart(N.trace(mm(MA1M,A1))))
        ehDamp2.append(checkImPart(N.trace(mm(MA2M,A2))))
        nHT.append(checkImPart(N.trace(mm(-2*MA2M+1.0j*(mm(MAM,G,G2)-mm(G2,Gd,MAM)),bA1))/2.0))
        nHTin.append(checkImPart(N.trace(mm(-2*MA2M,bA1))/2.0))
        nHTel.append(checkImPart(N.trace(mm(1.0j*(mm(MAM,G,G2)-mm(G2,Gd,MAM)),bA1))/2.0))
        HT.append(checkImPart(N.trace(mm((mm(G2,Gd,MA2mA1M)+mm(MA2mA1M,G,G2)),bA1))))
    # Write matrices to file?
    if WriteOrtho:
        ncdf.write(-1.0*N.array(nHT),'SymTrace')
        ncdf.write(HT,'AsymTrace')
        ncdf.close()
    
    NCfile.close()

    myGF.dGnin, myGF.dGnout = dGnin, dGnout
    myGF.P1T, myGF.P2T = P1T, P2T
    myGF.ehDamp1, myGF.ehDamp2 = ehDamp1, ehDamp2
    myGF.nHT, myGF.HT = nHT, HT
    myGF.nHTin, myGF.nHTel = nHTin, nHTel
    myGF.hw = hw

    unitConv=1.602177e-19/N.pi/1.054572e-34
    print 'checking: nHT',N.array(myGF.nHT)*unitConv # OK
    print 'checking: HT',N.array(myGF.HT) # OK
    
    # Set up grid and Hilbert term
    kT = general.Temp/11604.0 # (eV)

    # Generate grid for numerical integration of Hilbert term    
    max_hw=max(hw)
    max_win=max(-general.minBias,max_hw)+20*kT+4*general.Vrms
    min_win=min(-general.maxBias,-max_hw)-20*kT-4*general.Vrms
    pts=int(N.floor((max_win-min_win)/kT*3))
    Egrid=N.array(range(pts),N.float)/pts*(max_win-min_win)+min_win
    print "LOE: Hilbert integration grid : %i pts [%f,%f]" % (pts,min(Egrid),max(Egrid))
    
    # Calculate the prefactors for the Hilbert and non-Hilbert terms
    # Also calculate the Hilbert transfer of the box on the grid for each mode
    hilb=[] # List of hilbert transforms
    ker=None
    for ii in range(len(hw)):
        hwi=hw[ii]
        tmp=MM.box(hwi,-hwi,Egrid,kT)
        tmp2, ker = MM.Hilbert(tmp,ker)
        hilb.append(tmp2)
        SIO.printDone(ii,len(hw),'LOE : Hilbert transform')
    
    NN = general.biasPoints
    print 'biaspoints',NN

    # Add some points for the Lock in broadening
    approxdV=(general.maxBias-general.minBias)/NN
    NN+=int(((8*general.Vrms)/approxdV)+.5)
    
    Vl=general.minBias-4*general.Vrms+ \
        (general.maxBias-general.minBias+8*general.Vrms)/NN*N.array(range(NN),N.float)
    
    InH=N.zeros((NN,),N.complex) # Non-Hilb-term current
    IH=N.zeros((NN,),N.complex) # Hilb-term current
    Pow=N.zeros((NN,),N.float) # Power
    nPhtot=N.zeros((NN,),N.float) # Number of phonons (sum over modes)
    nPhvsBias=N.zeros((NN,len(hw)),N.float) # Number of phonons
    nPh=N.zeros((len(hw),),N.float)         # Number of phonons                                  

    for iV in range(len(Vl)):
        SIO.printDone(iV,len(Vl),'LOE : Calc current')
        
        InH[iV], IH[iV], Pow[iV] = 0.0, 0.0, 0.0
        V=Vl[iV]
        
        kasse=MM.box(0,-V,Egrid,kT)

        for iOmega in range(len(hw)):
            ihw = hw[iOmega]
            if ihw>general.modeCutoff:
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
                
                if general.PhHeating: # Heating?
                    nPh[iOmega]=heat/(ihw*myGF.P1T[iOmega]/N.pi+general.PhExtDamp)+1.0/(exphw-1.0)
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
                
        InH[iV]+=V*myGF.trans0 # Last to reduce cancelation errors

    # Get the right units for gamma_eh, gamma_heat
    hbar=1.0545727e-34/1.6021773e-19 # hbar in eV S
    gamma_eh=N.zeros((len(hw),),N.float)
    gamma_heat=N.zeros((len(hw),),N.float)
    for iOmega in range(len(hw)):
        # Units [Phonons per Second per dN where dN is number extra phonons]
        gamma_eh[iOmega]=myGF.P1T[iOmega]*hw[iOmega]/N.pi/hbar
        # Units [Phonons per second per eV [eV-ihw]
        gamma_heat[iOmega]=myGF.P2T[iOmega]/N.pi/hbar

    print 'checking: gamma_eh',gamma_eh # OK
    print 'checking: gamma_heat',gamma_heat # OK

    #hw, T, nHT, HT, lLOE, nPhtot, nPh, = hw, myGF.trans0, myGF.nHT, myGF.HT, [Vl, InH, IH], nPhtot, nPhvsBias
    V, I, dI, ddI, BdI, BddI, NnPhtot,NnPh = Broaden(general,Vl,InH+IH,nPhtot,nPhvsBias)

    print 'checking: BddI',BddI[:10] # OK


    datafile=general.systemlabel+'.IN'
    initncfile(datafile,hw)
    writeLOEData2Datafile(datafile,hw,myGF.trans0,nHT,HT)
    writeLOE2ncfile(datafile,hw,nHT,HT,V,I,NnPhtot,NnPh,\
                    dI,ddI,BdI,BddI,gamma_eh,gamma_heat)
                            

#    return hw, myGF.trans0, myGF.nHT, myGF.HT, [Vl, InH, IH], nPhtot, \
#           nPhvsBias, gamma_eh, gamma_heat
    


    
########################################################
def checkImPart(x):
    if abs(x.imag)>0.0000001:
        print "LOE : Imaginary part (%.3e) too big"%x.imag
        kuk
    return x.real   

########################################################
def writeFGRrates():
    # This does not work at the moment....
    unitConv=1.602177e-19/N.pi/1.054572e-34

    NCfile = NC.NetCDFFile(general.PhononNetCDF,'r')
    print 'Reading ',general.PhononNetCDF

    general.iChan, general.iSide = 0, 2
    outFile = file(general.systemlabel+'.IN.FGR','w')
    outFile.write('Total transmission [in units of (1/s/eV)] : %e\n' % (unitConv*myGF.totTrans.real,))

    tmp=N.sort(abs(N.array(myGF.nHT[:])))
    SelectionMin=tmp[-general.NumPhCurr]        
    
    for ihw in range(len(myGF.hw)):
        SIO.printDone(ihw,len(myGF.hw),'Golden Rate') 
        M = N.array(NCfile.variables['He_ph'][ihw,general.iSpin,:,:])
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
        inter = inter/IETS.P2T[ihw]
        intra = intra/IETS.P2T[ihw]
        outFile.write('Interchannel ratio: Sum(inter)/Tr[MA1MA2]      = %.4f \n'%inter)
        outFile.write('Intrachannel ratio: Sum(intra)/Tr[MA1MA2]      = %.4f \n'%intra)
        outFile.write('Inter+intra ratio: Sum(inter+intra)/Tr[MA1MA2] = %.4f \n'%(inter+intra))
        for iL in range(len(IETS.ECleft)):
            for iR in range(len(IETS.ECright)):
                outFile.write('%e ' % (unitConv*rate[iL,iR],))
            outFile.write('\n')
        outFile.close()
        NCfile.close()
                        
    
    


###########################################################
def setupParameters():
    # Note: can be called from Inelastica as well!
    usage = "usage: %prog [options]"
    description = """
Inelastica script calculates and writes LOE quantities in ascii (Systemlabel.IN) and NetCDF (Systemlabel.IN.nc) 

For help use --help!
"""
    parser = OptionParser(usage,description=description)
    parser.add_option("-n", "--NumChan", dest="numchan", help="Number of eigenchannels [%default]",
                      type='int', default=4)
    parser.add_option("-e", "--Energy", dest='energy', default=0.0,type='float',
                      help="Energy where eigenchannel scattering states are evaluated [%default eV]")
    parser.add_option("--eta", dest='eta', default=0.000001,type='float',
                      help="Tiny imag. part in Green's functions etc. [%default eV]")
    parser.add_option("-f", "--fdf", dest='fdfFile', default='./RUN.fdf',type='string',
                      help="fdf file used for transiesta calculation [%default]")
    parser.add_option("-s", "--iSpin", dest='iSpin', default=0,type='int',
                      help="Spin channel [%default]")
    parser.add_option("-k", "--kPoint", dest='kPoint', default='[0.0,0.0]',type='string',
                      help="2D k-point in range [0,1[ given as string, default=%default")
    parser.add_option("-p", "--PhononNetCDF", dest='PhononNetCDF', default=None,type='string',
                      help="Electron-phonon coupling NetCDF [%default]")
    parser.add_option("-t", "--Temp", dest='Temp', default=4.2,type='float',
                  help="Temperature [%default K]")
    parser.add_option("-b", "--BiasPoints", dest='biasPoints', default=801,type='int',
                  help="Number of bias points [%default]")
    parser.add_option("-v", "--MinMaxVoltage", dest='MinMaxVoltage', default='-0.4:0.4',type='string',
                  help="Voltage range ['%default' V]")
    parser.add_option("-c", "--ModeCutoff", dest='modeCutoff', default='0.0025',type='float',
                  help="Ignore phonon modes with lower hw [%default eV]")
    parser.add_option("-V", "--Vrms", dest='Vrms', default='0.005',type='float',
                  help="Lock in amplifier broadening [%default V]")
    parser.add_option("-H", "--Heating", dest='PhHeating', default=False,action='store_true',
                  help="Include heating of vibrational modes [%default]")
    parser.add_option("-x", "--PhExtDamp", dest='PhExtDamp', default=1e-15,type='float',
                  help="External damping [%default (?) TODO check unit!]")
    (general, args) = parser.parse_args()
    print description
    
    if general.kPoint[0]!='[' or general.kPoint[-1]!=']':
        parser.error("ERROR: please specify --kPoint='[x.x,y.y]' not --kPoint=%s"%general.kPoint)
    else:
        try:
            tmp=string.split(general.kPoint[1:-1],',')
            general.kPoint = N.array([float(tmp[0]),float(tmp[1]),0],N.float)
        except:
            parser.error("ERROR: please specify --kPoint='[x.x,y.y]' not --kPoint=%s"%general.kPoint)

    if not os.path.exists(general.fdfFile):
        parser.error("No input fdf file found, specify with --fdf=file.fdf (default RUN.fdf)")

    try:
        general.from_atom = SIO.GetFDFlineWithDefault(
            general.fdfFile,'TS.TBT.PDOSFrom', int, None, 'Eigenchannels')
        general.to_atom = SIO.GetFDFlineWithDefault(
            general.fdfFile,'TS.TBT.PDOSTo', int, None, 'Eigenchannels')
    except:
        parser.error("Specify device region with TS.TBT.PDOS[To/From] keyword.")
    try:
        tmp=string.split(general.MinMaxVoltage,':')
        general.minBias = float(tmp[0])
        general.maxBias = float(tmp[1])
    except:
        parser.error("ERROR: Inelastica failed to parse bias voltage!")

    return general


################################################################
# Broadening due to Vrms
################################################################

def Broaden(general,VV,II,nPhtot,nPh):
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
    VA=N.sqrt(2.0)*general.Vrms 

    # New bias ranges for broadening
    tmp=int(N.floor(VA/(dV[1]-dV[0]))+1)
    BdV=dV[tmp:-tmp]
    BddV=ddV[tmp:-tmp]

    # Initiate derivatives
    BdI=0*BdV
    BddI=0*BddV
    
    # Calculate first derivative with Vrms broadening
    for iV, V in enumerate(BdV):
        wt=(N.array(range(200))/200.0-0.5)*N.pi
        VL=V+VA*N.sin(wt)
        dIL=MM.interpolate(VL,dV,dI)
        BdI[iV]=2/N.pi*N.sum(dIL*(N.cos(wt)**2))*(wt[1]-wt[0])

    # Calculate second derivative with Vrms broadening    
    for iV, V in enumerate(BddV):
        wt=(N.array(range(200))/200.0-0.5)*N.pi
        VL=V+VA*N.sin(wt)
        ddIL=MM.interpolate(VL,ddV,ddI)
        BddI[iV]=8.0/3.0/N.pi*N.sum(ddIL*(N.cos(wt)**4))*(wt[1]-wt[0])

    # Reduce to one voltage grid
    NN=general.biasPoints 
    V=general.minBias+(general.maxBias-general.minBias)/NN*N.array(range(NN)) 

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
    tmp = ncfile.createVariable('IelasticL','d',('Bias',))
    tmp.units = 'Conductance quantum'
    tmp1 = ncfile.createVariable('IelasticR','d',('Bias',))
    tmp1.units = 'Conductance quantum'
    tmp2 = ncfile.createVariable('IinelasticL','d',('Bias',))
    tmp2.units = 'Conductance quantum'
    tmp3 = ncfile.createVariable('IinelasticR','d',('Bias',))
    tmp3.units = 'Conductance quantum'
    tmp4 = ncfile.createVariable('Bias','d',('Bias',))
    tmp4.units = 'V'
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
    unitConv=1.602177e-19/N.pi/1.054572e-34
    ncTotR = ncfile.createVariable('LOE_ISymTr','d',('Nph',))
    ncTotR[:] = N.array(nHT)*unitConv
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

def writeData2ncfile(filename,dict,ibias):
    'Append data to netCDF file'
    print 'Inelastica: Writing SCBA data to nc-file'
    bias = dict['bias']
    IelasticL = dict['NCL'] # Elastic current [Units, e/Pi/hbar, i.e.,
    #                         energy in eV, bias in V -> Current * 77 muA]
    IelasticR = (-dict['NCR'])
    IinelasticL = (dict['ICL']) # Left Current
    IinelasticR = (-dict['ICR']) # Right Current
    Nph=dict['nPh'] # Number of phonons
    # Add new data (some nicer way should exist)
    
    ncfile = NC.NetCDFFile(filename+'.nc','a')
    ncIelasticL = ncfile.variables['IelasticL']
    ncIelasticR = ncfile.variables['IelasticR']
    ncIinelasticL = ncfile.variables['IinelasticL']
    ncIinelasticR = ncfile.variables['IinelasticR']
    ncBias = ncfile.variables['Bias']
    ncNph = ncfile.variables['SCBA_nPh']

    tmp = ncBias[:]
    ncBias[len(tmp)] = bias

    tmp = ncIelasticL[:]
    ncIelasticL[len(tmp)]=IelasticL
    ncIelasticR[len(tmp)]=IelasticR
    ncIinelasticL[len(tmp)]=IinelasticL
    ncIinelasticR[len(tmp)]=IinelasticR
    ncNph[len(tmp),:]=Nph
    ncfile.close()

def writeData2Datafile(file,dict,hw):
    c1 = dict['bias']
    c2 = dict['NCL'] # Currents in units of e/pi/hbar*energy, i.e.,
    #                  G=77 muA/V * dI/dV if energy [eV], bias [V]
    c3 = (dict['ICL']-dict['ICR'])/2 # Average of left and right
    c4 = dict['NPL']+dict['NPR']
    c5 = dict['IPL']+dict['IPR']
    data = [c1,c2,c3,c4,c5]
    form = '%e   %e   %e   %e   %e   '
    for i in range(len(hw)):
        data.append(dict['nPh'][i])
        form += '%f   '
    form += '\n'
    data = tuple(data)
    f = open(file,'a')
    f.write(form %data)
    f.close()




################# Math helpers ################################
mm = MM.mm
outerAdd = MM.outerAdd
dagger = MM.dagger
    
##################### Start main routine #####################

if __name__ == '__main__':
    main()

