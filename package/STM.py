print "SVN $Id$"

"""
Eigenchannels:
1: Eigenchannels, method from Paulsson and Brandbyge PRB 2007
2: Analyze IETS spectra, Paulsson et al PRL 2008
3: Calculate "bond" currents
"""
fudgeEnergy = 5e-5         # Fudge energy in eV used to see if Ef is close to the calculated energy

import pyTBT as pyTBT
import SiestaIO as SIO
import NEGF
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
    calcSTM()

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
    class options:
        pass
    options.fn='RUN.fdf'
    options.fnL  = './'+SIO.GetFDFlineWithDefault(options.fn,'TS.HSFileLeft', str, None, 'pyTBT')
    options.NA1L = SIO.GetFDFlineWithDefault(options.fn,'TS.ReplicateA1Left', int, 1, 'pyTBT')
    options.NA2L = SIO.GetFDFlineWithDefault(options.fn,'TS.ReplicateA2Left', int, 1, 'pyTBT')
    options.fnR  = './'+SIO.GetFDFlineWithDefault(options.fn,'TS.HSFileRight', str, None, 'pyTBT')
    options.NA1R = SIO.GetFDFlineWithDefault(options.fn,'TS.ReplicateA1Right', int, 1, 'pyTBT')
    options.NA2R = SIO.GetFDFlineWithDefault(options.fn,'TS.ReplicateA2Right', int, 1, 'pyTBT')

    options.systemlabel = SIO.GetFDFlineWithDefault("RUN.fdf",'SystemLabel', str, 'Systemlabel', 'Eigenchannels')       
    options.TSHS = './%s.TSHS'%(options.systemlabel)

    # Setup H, S and self-energies
    # Setup self-energies and device GF
    elecL = NEGF.ElectrodeSelfEnergy(options.fnL,options.NA1L,options.NA2L,0)
    elecR = NEGF.ElectrodeSelfEnergy(options.fnR,options.NA1R,options.NA2R,0)
    myGF = NEGF.GF(options.TSHS,elecL,elecR,Bulk=False,DeviceAtoms=[general.from_atom, general.to_atom])

    HS.L, HS.R, HS.GF, devSt, devEnd, general.Elist, general.eta, general.SiestaOutFile = \
        elecL, elecR, myGF, myGF.DeviceOrbs[0], myGF.DeviceOrbs[1], [0.0], 1e-6, options.systemlabel

########################################################
def calcSTM():
    # Calculate STM picture
    Ef = general.Ef

    # Find division between L-R
    zcoord = N.array(geom.xyz)[:,2]
    zdiff = zcoord[1:]-zcoord[:-1]
    lastLeftAtom = N.where(zdiff==N.max(zdiff))[0][0]
    tipxyz = geom.xyz[lastLeftAtom+1] # Guess at tip possition
    lastLeftOrb = HS.GF.HS.lasto[lastLeftAtom+1]-HS.GF.DeviceOrbs[0]
    firstRightOrb = lastLeftOrb+1

    HS.GF.calcGF(Ef+general.eta*1j, general.kPoint[0:2], ispin=general.iSpin)

    # Calculate transmission with coupling L-R
    transmission = HS.GF.calcT()

    # Remove coupling L-R
    H, S = HS.GF.H, HS.GF.S
    H[:lastLeftOrb+1,firstRightOrb:]=0
    H[firstRightOrb:,:lastLeftOrb+1]=0
    S[:lastLeftOrb+1,firstRightOrb:]=0
    S[firstRightOrb:,:lastLeftOrb+1]=0

    nnL,  nnR = len(HS.GF.GamL),  len(HS.GF.GamR)
    Gam1,  Gam2 = N.zeros(HS.GF.H.shape, N.complex), N.zeros(HS.GF.H.shape, N.complex)
    Sig1,  Sig2 = N.zeros(HS.GF.H.shape, N.complex), N.zeros(HS.GF.H.shape, N.complex)
    Gam1[0:nnL, 0:nnL] ,  Gam2[-nnR:, -nnR:] = HS.GF.GamL, HS.GF.GamR
    Sig1[0:nnL, 0:nnL] ,  Sig2[-nnR:, -nnR:] = HS.GF.SigL, HS.GF.SigR

    Gr = LA.inv(Ef*N.eye(len(H))-H-Sig1-Sig2)

    rho1= mm(Gr,Gam1,dagger(Gr),S)/(2*N.pi)
    rho2= mm(Gr,Gam2,dagger(Gr),S)/(2*N.pi)
    eval1, evec1 = LA.eig(rho1)
    eval2, evec2 = LA.eig(rho2)

    for ii in range(len(eval1)):
        evec1[:,ii] = evec1[:,ii]/N.sqrt(N.dot(N.conjugate(evec1[:,ii]),evec1[:,ii]))
        evec2[:,ii] = evec2[:,ii]/N.sqrt(N.dot(N.conjugate(evec2[:,ii]),evec2[:,ii]))

    Y1, dstep, origo, nx, ny, nz = calcWF(evec1[:,0])
    rho1 =N.real(eval1[0]*Y1*N.conjugate(Y1))
    Y2, dstep, origo, nx, ny, nz = calcWF(evec2[:,0])
    rho2 =N.real(eval2[0]*Y2*N.conjugate(Y2))

    for ii in range(1,len(eval1)):
        SIO.printDone(ii,len(eval1),'Real space density generation')
        if eval1[ii]>1e-8:
            Y, dstep, origo, nx, ny, nz = calcWF(evec1[:,ii])
            rho1 += N.real(eval1[ii]*Y*N.conjugate(Y))
        if eval2[ii]>1e-8:
            Y, dstep, origo, nx, ny, nz = calcWF(evec2[:,ii])
            rho2 += N.real(eval2[ii]*Y*N.conjugate(Y))

    fn = general.DestDir+'/'+general.SiestaOutFile
    if general.format == 'macu':
        writemacubin(fn+'.rho1.macu',rho1.real,nx,ny,nz,origo,dstep)
        writemacubin(fn+'.rho2.macu',rho2.real,nx,ny,nz,origo,dstep)
    if general.format == 'cube':
        writecube(fn+'.rho1.cube',rho1.real,nx,ny,nz,origo,dstep)
        writecube(fn+'.rho2.cube',rho2.real,nx,ny,nz,origo,dstep)
    if general.format == 'XSF':
        writeXSF(fn+'.rho1.XSF',rho1,nx,ny,nz,origo,dstep)
        writeXSF(fn+'.rho2.XSF',rho2,nx,ny,nz,origo,dstep)
    if general.format == 'nc':
        writenetcdf(fn+'.rho1.nc',rho1,nx,ny,nz,origo,dstep)
        writenetcdf(fn+'.rho2.nc',rho2,nx,ny,nz,origo,dstep)

    upperLim1=N.max(N.where(rho1>0)[2])
    lowerLim2=N.min(N.where(rho2>0)[2])

    rho1=rho1[:,:,lowerLim2:upperLim1+1]
    rho2=rho2[:,:,lowerLim2:upperLim1+1]

    origo=origo+dstep*(lowerLim2+(upperLim1-lowerLim2+1)/2)*\
        N.array([0,0,1],N.float)

    mvx, mvy = tipxyz[0]-origo[0], tipxyz[1]-origo[1]
    rho1 = N.roll(rho1, -int(mvx/dstep), axis=0)
    rho1 = N.roll(rho1, -int(mvy/dstep), axis=1)
    origo[0], origo[1] = tipxyz[0]-int(mvx/dstep)*dstep, tipxyz[1]-int(mvy/dstep)*dstep
    Nz = general.NumZ

    image = N.zeros((nx,ny,Nz))
    for iz in range(Nz):
        for iy in range(ny):
            for ix in range(nx):
                SIO.printDone(iz*nx*ny+iy*nx+ix,Nz*nx*ny,'Image generation')
                image[ix,iy,iz]=N.sum(rho1*rho2)
                rho2 = N.roll(rho2, 1, axis=0)
            rho2 = N.roll(rho2, 1, axis=1)
        rho2 = N.roll(rho2, 1, axis=2) # Move tip up
        rho2[:,:,0]=0.0

    image=image/N.max(image)

    if general.format == 'macu':
        writemacubin(fn+'.STM.macu',image,nx,ny,Nz,origo,dstep)
    if general.format == 'cube':
        writecube(fn+'.STM.cube',image,nx,ny,Nz,origo,dstep)
    if general.format == 'XSF':
        writeXSF(fn+'.STM.XSF',image,nx,ny,Nz,origo,dstep)
    if general.format == 'nc':
        writenetcdf(fn+'.STM.nc',image,nx,ny,Nz,origo,dstep)
    


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



###########################################################
def setupParameters():
    # Note: can be called from Inelastica as well!
    global general
    
    usage = "usage: %prog [options] DestinationDirectory"
    description = """
STM image generation program
For help use --help!
"""
    parser = OptionParser(usage,description=description)
    EC = OptionGroup(parser, "Options for STM")
    EC.add_option("-r", "--Res", dest='res', default=0.4,type='float',
                  help="Resolution [%default Ang]")
    EC.add_option("-w", "--format", dest='format', default='macu',type='string',
                  help="Wavefunction format (macu, cube, XSF, or nc) [%default]")
    EC.add_option("-N", "--NumZ", dest='NumZ', default=10,type='int',
                  help="Also move the tip in the Z-direction away from sample with NumZ steps [%default]")
    parser.add_option_group(EC)
    
    EC.add_option("-e", "--Ef", dest='Ef', default=0.0,type='float',
                       help="Fermi energy [%default eV]")
    EC.add_option("-f", "--fdf", dest='fdfFile', default='./RUN.fdf',type='string',
                       help="fdf file used for transiesta calculation [%default]")
    EC.add_option("-s", "--iSpin", dest='iSpin', default=1,type='int',
                       help="Spin channel [%default]")
    EC.add_option("-k", "--kPoint", dest='kPoint', default='[0,0]',type='string',
                       help="2D k-point in range [0,1] given as string! default='[0.0,0.0]'")

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
    myprint('Res [Ang]        : %f'%general.res,file)
    myprint('iSpin            : %i'%(general.iSpin+1),file)
    myprint('kPoint           : [%f,%f]'%(general.kPoint[0],general.kPoint[1]),file)
    myprint('Device [from,to] : [%i,%i]'%(general.from_atom, general.to_atom),file)
    myprint('DestDir          : %s'%general.DestDir,file)
    myprint('##################################################################################',file)
    file.close()

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

