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
    general.nspin = myGF.HS.nspin
    basis = readbasis(general,myGF.HS)
    calcT(general,geom,myGF,basis)

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

    myGF = NEGF.GF(general.systemlabel+'.TSHS',elecL,elecR,Bulk=True,DeviceAtoms=[general.from_atom, general.to_atom])

    myGF.calcGF(general.energy+general.eta*1j, general.kPoint[0:2], ispin=general.iSpin)

    return myGF

########################################################
def calcT(general,geom,myGF,basis):
    # Matrix to save total and eigenchannel transmissions
    # BEFORE ORTHOGO
    T = myGF.calcT(general.numchan)
    print 'Transmission T(E=%.4f) [Ttot, T1, T2, ... Tn]:'%general.energy
    for t in T:
        print '%.9f '%t,
    print

    # Now orthogonalize in device region
    myGF.orthogonalize()

    # Check that we get the same transmission:
    T = myGF.calcT(general.numchan)
    print 'Transmission T(E=%.4f) [Ttot, T1, T2, ... Tn]:'%general.energy
    for t in T:
        print '%.9f '%t,
    print

    # Calculate Eigenchannels
    myGF.calcEigChan(general.numchan)

    # Calculate Eigenchannels from left
    ECleft, EigT = myGF.ECleft, myGF.EigTleft
    print 'Left eigenchannel transmissions:',EigT
    for jj in range(general.numchan):
        T[jj+1]=EigT[len(EigT)-jj-1]

    # Write wave functions
    for jj in range(general.numchan):
        general.iSide, general.iChan = 0, jj+1
        writeWavefunction(general,geom,basis,ECleft[jj])
        Curr=calcCurrent(general,basis,myGF.HNO,ECleft[jj])
        writeCurrent(general,geom,Curr)
                    
    # Calculate eigenchannels from right
    if general.bothsides:
        ECright, EigT = myGF.ECright, myGF.EigTright
        print 'Right eigenchannel transmissions:',EigT
        for jj in range(general.numchan):
            general.iSide, general.iChan = 1, jj+1
            writeWavefunction(general,geom,basis,ECright[jj])
            Curr=calcCurrent(general,basis,myGF.HNO,ECright[jj])
            writeCurrent(general,geom,Curr)
            
    # Calculate total "bond currents"
    A1NO = mm(myGF.Us,myGF.A1,myGF.Us)
    A2NO = mm(myGF.Us,myGF.A2,myGF.Us)
    Curr=-calcCurrent(general,basis,myGF.HNO,A1NO)
    general.iChan, general.iSide = 0, 0
    writeCurrent(general,geom,Curr)
    Curr=-calcCurrent(general,basis,myGF.HNO,A2NO)
    general.iSide = 1
    writeCurrent(general,geom,Curr)

    # Calculate eigenstates of device Hamiltonian
    if general.MolStates>0.0:
        print 'calculating molecular eigenstates of folded, orthogonalized device hamiltonian'
        ev, es = LA.eigh(myGF.H)
        for ii in range(len(ev)):
            if N.abs(ev[ii])<general.MolStates:
                fn=general.DestDir+'/'+general.systemlabel+'.S%.3i.E%.3f'%(ii,ev[ii])
                writeWavefunction(general,geom,basis,mm(myGF.Us,es[:,ii]),fn=fn)



########################################################
def calcWF(general,geom,basis,Y):
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

        # Calculate spherical harmonics
        l = basis.L[ii]
        m = basis.M[ii]
        if l==3:
            print 'f-shell : l=%i, m=%i (NOT TESTED!!)'%(l,m)
        SphHar = MM.sphericalHarmonics(sinth,costh,sinfi,cosfi)
        thisSphHar = SphHar[l][m+l]

        YY[ixmin:ixmax,iymin:iymax,izmin:izmax]=YY[ixmin:ixmax,iymin:iymax,izmin:izmax]+\
                                                 RR*thisSphHar*Y[ii]

    return YY, general.res, origo, nx, ny, nz

###################### Current densities #########################

def calcCurrent(general,basis,H,Y):
    """
    Calculate current density in atomic bonds
    Y : complex scattering state or
    Y : A_l or A_r! (for total current)
    """
    
    NN=len(H)
    NN2=general.to_atom-general.from_atom+1
    Curr=N.zeros((NN2,NN2),N.float)
    
    if len(Y.shape)==2:
        for ii in range(NN):
            a1=basis.ii[ii]-general.from_atom
            for jj in range(NN):
                a2=basis.ii[jj]-general.from_atom
                tmp=H[jj,ii]*Y[ii,jj]/2/N.pi
                Curr[a1,a2]=Curr[a1,a2]+4*N.pi*tmp.imag
    else:
        for ii in range(NN):
            a1=basis.ii[ii]-general.from_atom
            for jj in range(NN):
                a2=basis.ii[jj]-general.from_atom
                tmp=H[ii,jj]*N.conjugate(Y[ii])*Y[jj]
                Curr[a1,a2]=Curr[a1,a2]+4*N.pi*tmp.imag
    
    return Curr

########################################################
def writeCurrent(general,geom,Curr):
    """
    Write bond currents to file
    fn : Filename
    Curr : Bond current matrix
    """
    fn=fileName(general)
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
def writenetcdf(geom,fn,YY,nx,ny,nz,origo,dstep):
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
def writecube(geom,fn,YY,nx,ny,nz,origo,dstep):
    """
    Write wavefunction to cube file.
    """
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
def writeXSF(geom,fn,YY,nx,ny,nz,origo,dstep):
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
def writeWavefunction(general,geom,basis,Y,fn=None):
    """
    Writes wave function to the specified output type
    Y: vector for wavefunction
                
    """
    if fn==None: fn = fileName(general)
    print 'Eigenchannels.writeWavefunction: Writing',fn
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

    YY, dstep, origo, nx, ny, nz = calcWF(general,geom,basis,Y)

    # Write wave function in specified file format
    if general.format.lower() == 'macu':
        writemacubin(fn+'.Re.macu',YY.real,nx,ny,nz,origo,dstep)
        writemacubin(fn+'.Im.macu',YY.imag,nx,ny,nz,origo,dstep)

    if general.format.lower() == 'cube':
        writecube(geom,fn+'.Re.cube',YY.real,nx,ny,nz,origo,dstep)
        writecube(geom,fn+'.Im.cube',YY.imag,nx,ny,nz,origo,dstep)
        
    if general.format.lower() == 'xsf':
        writeXSF(geom,fn+'.XSF',YY,nx,ny,nz,origo,dstep)

    if general.format.lower() == 'nc':
        writenetcdf(geom,fn+'.nc',YY,nx,ny,nz,origo,dstep)


###########################################################
def setupParameters():
    # Note: can be called from Inelastica as well!
    usage = "usage: %prog [options] DestinationDirectory"
    description = """
Eigenchannels script that calculates:
1) Eigenchannels, Paulsson et al. PRB 76, 115117 (2007).
2) Fermi golden rule scattering rates, Paulsson et al. PRL 100, 226604 (2008).
3) Bond currents.

For help use --help!
"""
    parser = OptionParser(usage,description=description)
    parser.add_option("-n", "--NumChan", dest="numchan", help="Number of eigenchannels [%default]", 
                  type='int', default=4)
    parser.add_option("-B", "--BothSides", dest='bothsides', default=False,action='store_true',
                  help="Calculate eigenchannels from both sides [%default]")
    parser.add_option("-M", "--MolecularStates", dest='MolStates', default=0.0, type='float',
                  help="Calculate eigenstates of device Hamiltonian within [%default] eV from Ef")
    parser.add_option("-r", "--Res", dest='res', default=0.4,type='float',
                  help="Resolution [%default Ang]")
    parser.add_option("-w", "--format", dest='format', default='macu',type='string',
                  help="Wavefunction format (macu, cube, XSF, or nc) [%default]")
    parser.add_option("-N", "--NumPhCurr", dest='NumPhCurr', default=10,type='int',
                  help="Max number of changes in bond currents from inelastic scattering  [%default]")
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

    (general, args) = parser.parse_args()
    print description
    
    # Sanity checks

    # Check k-point input
    if general.kPoint[0]!='[' or general.kPoint[-1]!=']':
        parser.error("ERROR: please specify --kPoint='[x.x,y.y]' not --kPoint=%s"%general.kPoint)
    else:
        try:
            tmp=string.split(general.kPoint[1:-1],',')
            general.kPoint = N.array([float(tmp[0]),float(tmp[1]),0],N.float)
        except:
            parser.error("ERROR: please specify --kPoint='[x.x,y.y]' not --kPoint=%s"%general.kPoint)

    # Check for Siesta input file
    if not os.path.exists(general.fdfFile):
        parser.error("No input fdf file found, specify with --fdf=file.fdf (default RUN.fdf)")

    # Find Systemlabel keyword
    general.systemlabel = SIO.GetFDFlineWithDefault(general.fdfFile,'SystemLabel', str, 'Systemlabel', 'Error Eigenchannels')
        
    # Find TS.TBT.PDOSTo/From keywords
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
        print '\nEigenchannels : Creating folder %s' %general.DestDir
        os.mkdir(general.DestDir)
#    else:
#        parser.error('ERROR: destination directory %s already exist!'%general.DestDir)
    return general

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
    myprint('iSpin            : %i'%(general.iSpin),file)
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
def fileName(general):
    systemlabel = general.systemlabel
    # Generate filename '.EC.{1,Tot}{L,R,,In,Out}[UP][Ef=1.0].'
    if general.iChan==0:
        fn=systemlabel+'.EC.Tot%s'%(['L','R','','In','Out'][general.iSide])
    else:
        fn=systemlabel+'.EC.%i%s'%(general.iChan,['L','R','','In','Out'][general.iSide])
    if general.nspin==2:
        fn += ['UP','DOWN'][general.iSpin]
    if general.energy!=0.0:
        fn += '_E%.3f'%general.energy
    if general.kPoint[0]!=0.0 or general.kPoint[1]!=0.0:
        fn += '_kx%.3f_ky%.3f'%(general.kPoint[0],general.kPoint[1])
    return general.DestDir+'/'+fn

################# Math helpers ################################
mm = MM.mm
outerAdd = MM.outerAdd
dagger = MM.dagger
    
##################### Start main routine #####################

if __name__ == '__main__':
    main()

