print "SVN $Id$"

"""
Eigenchannels:
1: Eigenchannels, method from Paulsson and Brandbyge PRB 2007
2: Analyze IETS spectra, Paulsson et al PRL 2008
3: Calculate "bond" currents
"""

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


########################################################
##################### Main routine #####################
########################################################
def main():
    general = setupParameters()
    geom = readxv()
    myGF = readHS(general)
    basis = readbasis(general,myGF.HS)
    calcTeig(general,myGF)
    calcIETS(general,myGF)

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

    systemlabel = SIO.GetFDFlineWithDefault(general.fdfFile,'SystemLabel', str, None, 'pyTBT')

    myGF = NEGF.GF(systemlabel+'.TSHS',elecL,elecR,Bulk=True,DeviceAtoms=[general.from_atom, general.to_atom])

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

    # Calculate Eigenchannels from left
    myGF.A1 = mm(myGF.Gr,myGF.GamL,myGF.Ga)
    myGF.ECleft, EigT = myGF.calcEigChan(myGF.A1,myGF.GamR,general.numchan)
    print 'Left eigenchannel transmissions:',EigT

    # Calculate Eigenchannels from right
    myGF.A2 = mm(myGF.Gr,myGF.GamR,myGF.Ga)
    myGF.ECright, EigT = myGF.calcEigChan(myGF.A2,myGF.GamL,general.numchan)
    print 'Right eigenchannel transmissions:',EigT

########################################################
def calcIETS(general)
    # Calculate inelastic scattering rates, total and per eigenchannel
    IETS.totTrans = T[iiE,0]

    NCfile = NC.NetCDFFile(general.PhononNetCDF,'r')
    print 'Reading ',general.PhononNetCDF
    hw = N.array(NCfile.variables['hw'][:]) 
    # Write matrices to file?
    WriteOrtho = False
    if WriteOrtho:
        ncdf = NCDF.NCfile(general.DestDir+'/ortho.nc')
        ncdf.write(hw,'hw')
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
    # Write matrices to file?
    if WriteOrtho:
        ncdf.write(-1.0*N.array(nHT),'SymTrace')
        ncdf.write(HT,'AsymTrace')
        ncdf.close()
    
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
        inter,intra = 0.0, 0.0 # splitting total rate in two
        for iL in range(len(IETS.ECleft)):
            for iR in range(len(IETS.ECright)):
                tmp=N.dot(N.conjugate(IETS.ECleft[iL]),mm(M,IETS.ECright[iR]))
                rate[iL,iR]=(2*N.pi)**2*abs(tmp)**2
                totrate+=rate[iL,iR]
                if iL==iR: intra += rate[iL,iR]
                else: inter += rate[iL,iR]

        if abs(IETS.nHT[ihw])>=SelectionMin:
            general.iChan = ihw
          
    
    


###########################################################
def setupParameters():
    # Note: can be called from Inelastica as well!
    usage = "usage: %prog [options] DestinationDirectory"
    description = """
Inelastica script calculates:
1) ...

For help use --help!
"""
    parser = OptionParser(usage,description=description)
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

    if len(args)!=1:
        parser.error('ERROR: You need to specify destination directory')
    general.DestDir = args[0]
    if not os.path.isdir(general.DestDir):
        print '\nInelastica : Creating folder %s' %general.DestDir
        os.mkdir(general.DestDir)
#    else:
#        parser.error('ERROR: destination directory %s already exist!'%general.DestDir)
    return general


################# Math helpers ################################
mm = MM.mm
outerAdd = MM.outerAdd
dagger = MM.dagger
    
##################### Start main routine #####################

if __name__ == '__main__':
    main()

