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

################### Help functions ############################
try:
    import F90helpers as F90
    F90imported = True
except:
    F90imported = False
    print "########################################################"
    print "Perhaps time to compile F90/setkpointhelper"
    print "Try:" 
    print "        cd F90;source compile.bat"
    print "########################################################"



################### Main program ############################
def main():
    usage = "usage: %prog [options] DestinationDirectory"
    description = """pyTBT is the Python version of TBtrans originally developed by Mads Brandbyge.
For help use --help!
 """
    parser = OptionParser(usage,description=description)

    # Determine keywords provided
    parser.add_option("-n", "--NumChan", dest="numchan", help="Number of eigenchannels [%default]",
                      type='int', default=10)
    parser.add_option("-e","--eta", dest="eta", help="Imaginary part in self-energies [%default eV]",
                      type='float', default=0.000001)
    parser.add_option("-x","--Nk1", dest='Nk1', default=1,type='int',
                      help="k-points Nk1 [%default]")
    parser.add_option("-y","--Nk2", dest='Nk2', default=1,type='int',
                      help="k-points Nk2 [%default]")
    parser.add_option("-s", "--sym", dest='symmetry',default=False,action='store_true',
                      help="Use time reversal symmetry to reduce number of k-points (Nk2 points sample only the range [0,0.5]) [%default]")
    parser.add_option("-g", "--avoid-gamma", dest='skipgamma',default=False,action='store_true',
                      help="Avoid gamma point in k-point sampling [%default]")
    parser.add_option("-f", "--fdf", dest='fn',default='./RUN.fdf',type='string',
                      help="Input fdf-file for TranSIESTA calculations [%default]")
    parser.add_option("-d", "--skip-dos", dest='dos',default=True,action='store_true',
                      help="Calculate DOS calculation? [%default]")

    (general, args) = parser.parse_args()
    print description

    if len(args)!=1:
        parser.error('ERROR: You need to specify destination directory')
    general.DestDir = args[0]

    if not os.path.isdir(general.DestDir):
        print '\npyTBT: Creating folder %s' %general.DestDir
        os.mkdir(general.DestDir)

    # Make sure to avoid Gamma point if time-reversal symmetry is used
    if general.symmetry:
        print 'pyTBT: Applying time-reversal symmetry (Nk2=%i sampling in the range [0,0.5])'%general.Nk2
        if not general.skipgamma:
            print '... needs to avoid Gamma point. Using flag -g (skipgamma)'
            general.skipgamma = True
    
    # Read options from fdf files
    ##############################################################################
    
    fn = general.fn
    head,tail = os.path.split(fn)
    print "pyTBT: Reading keywords from %s \n"%fn
    
    # Electrodes
    fnL  =head+'/'+SIO.GetFDFlineWithDefault(fn,'TS.HSFileLeft', str, None, 'pyTBT')
    NA1L =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA1Left', int, 1, 'pyTBT')
    NA2L =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA2Left', int, 1, 'pyTBT')
    fnR  =head+'/'+SIO.GetFDFlineWithDefault(fn,'TS.HSFileRight', str, None, 'pyTBT')
    NA1R =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA1Right', int, 1, 'pyTBT')
    NA2R =SIO.GetFDFlineWithDefault(fn,'TS.ReplicateA2Right', int, 1, 'pyTBT')

    # Device region
    devSt =SIO.GetFDFlineWithDefault(fn,'TS.TBT.PDOSFrom', int, 0, 'pyTBT')
    devEnd=SIO.GetFDFlineWithDefault(fn,'TS.TBT.PDOSTo', int, 0, 'pyTBT')
    
    # Voltage
    voltage  =SIO.GetFDFlineWithDefault(fn,'TS.Voltage', float, 0.0, 'pyTBT')

    # Energy range
    nE  =SIO.GetFDFlineWithDefault(fn,'TS.TBT.NPoints', int, 21, 'pyTBT')
    minE=SIO.GetFDFlineWithDefault(fn,'TS.TBT.Emin', float, -1.0, 'pyTBT')
    maxE=SIO.GetFDFlineWithDefault(fn,'TS.TBT.Emax', float, 1.0, 'pyTBT')
    if nE>1:
        dE = (maxE-minE)/float(nE-1)
        Elist = N.array(range(int((maxE-minE+1e-9)/dE)+1),N.float)*dE+minE
    else:
        dE=0.0
        Elist=N.array((minE,),N.float)

    UseBulk=SIO.GetFDFlineWithDefault(fn,'TS.UseBulkInElectrodes', bool, True, 'pyTBT')

    #eta=SIO.GetFDFlineWithDefault(fn,'pyTBT.eta', float, 0.000001, 'pyTBT')
    #Nk1=SIO.GetFDFlineWithDefault(fn,'pyTBT.K_A1', int, 1, 'pyTBT')
    #Nk2=SIO.GetFDFlineWithDefault(fn,'pyTBT.K_A2', int, 1, 'pyTBT')
    eta = general.eta
    Nk1 = general.Nk1
    Nk2 = general.Nk2

    general.systemlabel = SIO.GetFDFlineWithDefault(fn,'SystemLabel', str, 'Systemlabel', 'pyTBT')

    outFile=general.DestDir+'/'+general.systemlabel + '.%ix%i'%(Nk1,Nk2)

    ##############################################################################
    # Define electrodes and device

    elecL = NEGF.ElectrodeSelfEnergy(fnL,NA1L,NA2L,voltage/2.)
    elecR = NEGF.ElectrodeSelfEnergy(fnR,NA1R,NA2R,-voltage/2.)
    myGF = NEGF.GF(head+'/%s.TSHS'%general.systemlabel,elecL,elecR,Bulk=UseBulk,DeviceAtoms=[devSt, devEnd])
    nspin = myGF.HS.nspin
    if devSt==0:
        devSt=GF.DeviceAtoms[0]
    if devEnd==0:
        devEnd=GF.DeviceAtoms[1]
        
    print """
##############################################################
pyTBT

Energy [eV]                     : %f:%f:%f
kpoints                         : %i, %i 
eta [eV]                        : %f 
Device [Atoms Siesta numbering] : %i:%i 
Bulk                            : %s
SpinPolarization                : %i
Voltage                         : %f
##############################################################

"""%(minE,dE,maxE,Nk1,Nk2,eta,devSt,devEnd,UseBulk,nspin,voltage)

    
    channels = general.numchan
    Tkpt=N.zeros((len(Elist),Nk1,Nk2,channels+1),N.float)
    DOSL=N.zeros((len(Elist),myGF.nuo),N.float)
    DOSR=N.zeros((len(Elist),myGF.nuo),N.float)
    # Loop over spin
    for iSpin in range(nspin):
        if nspin<2:
            fo=open(outFile+'.AVTRANS','write')
        else:
            fo=open(outFile+['.UP','.DOWN'][iSpin]+'.AVTRANS','write')
        fo.write('# Nk1=%i Nk2=%i eta=%.2e\n'%(Nk1,Nk2,eta))
        fo.write('# E   Ttot(E)   Ti(E) (i=1-10)\n')
        # Loop over energy
        for ie, ee in enumerate(Elist):
            Tavg = N.zeros(channels+1,N.float)
            AavL = N.zeros((myGF.nuo,myGF.nuo),N.complex)
            AavR = N.zeros((myGF.nuo,myGF.nuo),N.complex)
            # Loops over k-points
            for ik1 in range(Nk1):
                for ik2 in range(Nk2):
                    kpt=N.array([ik1/float(Nk1),ik2/float(Nk2)],N.float)
                    if general.skipgamma:
                        # Add \Delta k/2 to shift away from Gamma
                        kpt += N.array([1./float(2*Nk1),1./float(2*Nk2)],N.float)
                    if general.symmetry:
                        # Let Nk2 points sample only the range [0,0.5]
                        kpt[1] = kpt[1]/2
                    myGF.calcGF(ee+eta*1.0j,kpt,ispin=iSpin)
                    # Transmission:
                    T = myGF.calcT(channels)
                    Tavg += T/(Nk1*Nk2)
                    Tkpt[ie,ik1,ik2] = T
                    # DOS calculation:
                    if general.dos:
                        GamL, GamR, Gr = myGF.GamL, myGF.GamR, myGF.Gr
                        nuo, nuoL, nuoR = myGF.nuo, myGF.nuoL, myGF.nuoR
                        AL = MM.mm(Gr[:,0:nuoL],GamL,MM.dagger(Gr)[0:nuoL,:])
                        AR = MM.mm(Gr[:,nuo-nuoR:nuo],GamR,MM.dagger(Gr)[nuo-nuoR:nuo,:])
                        AavL += MM.mm(AL,myGF.S)
                        AavR += MM.mm(AR,myGF.S)
            # Print calculated quantities
            print ee, Tavg
            transline = '\n%.10f '%ee
            for ichan in range(channels+1):
                if ichan==0:
                    transline += '%.8e '%Tavg[ichan]
                else:
                    transline += '%.4e '%Tavg[ichan]
            fo.write(transline)
            # Partial density of states:
            if general.dos:
                DOSL[ie,:] += N.diag(AavL).real/(Nk1*Nk2*2*N.pi)
                DOSR[ie,:] += N.diag(AavR).real/(Nk1*Nk2*2*N.pi)
                print ee," ",N.sum(DOSL[ie,:]),N.sum(DOSR[ie,:]), '# DOS'
        fo.close()
        
        # Write k-point-resolved transmission
        if nspin<2:
            fo=open(outFile+'.TRANS','write')
        else:
            fo=open(outFile+['.UP','.DOWN'][iSpin]+'.TRANS','write')
        for ik1 in range(Nk1):
            for ik2 in range(Nk2):
                kpt=N.array([ik1/float(Nk1),ik2/float(Nk2)],N.float)
                if general.skipgamma:
                    kpt += N.array([1./float(2*Nk1),1./float(2*Nk2)],N.float)
                if general.symmetry:
                    kpt[1] = kpt[1]/2
                fo.write('\n\n# k = %f, %f '%(kpt[0],kpt[1]))
                for ie, ee in enumerate(Elist):
                    transline = '\n%.10f '%ee
                    for ichan in range(channels+1):
                        if ichan==0:
                            transline += '%.8e '%Tkpt[ie,ik1,ik2,ichan]
                        else:
                            transline += '%.4e '%Tkpt[ie,ik1,ik2,ichan]
                    fo.write(transline)
        fo.close()

    # Things only needed for SDOS script
    return elecL, elecR, myGF, devSt, devEnd, Elist,eta, general.systemlabel

if __name__ == '__main__':
    main()
