"""
Routines for IO in different formats:
1: Geometries in xyz, XV, ANI, fdf, mkl etc formats
2: Read Hamiltonian and Overlap from .TSHS file.
3: fdf manipulations
4: Real space Gaussian cube files
5: Obtain information of basis orbitals for calculation (.ion.nc)
"""
import numpy as N
import numpy.linalg as LA
import string, struct, os.path, sys
import gzip
import netCDF4 as NC4
import Inelastica.physics.constants as PC
import Inelastica.ValueCheck as VC

# For speed some routines can be linked as F90 code
try:
    import Inelastica.fortran.F90helpers as F90
    F90imported = True
except:
    F90imported = False
    print "########################################################"
    print "Problems encountered with F90helpers.so"
    print "Falling back on a pure python (slower) implementation"
    print "Try compiling manually following these steps:"
    print " $ cd Inelastica/fortran"
    print " $ source compile.bat (or compile_alternative.bat)"
    print " $ cp F90helpers.so <python>/site-packages/Inelastica/fortran"
    print "########################################################"

# Check length of int and long and use the one that has 8 bytes
if struct.calcsize('l')==8:
    fortranPrefix='='
    fortranuLong='I'
    fortranLong='i'    
else:
    fortranPrefix=''
    fortranuLong='I'
    fortranLong='i'    


def SIO_open(filename,mode='r'):
    "A io.siesta redefinition of the function open() to handle gzip format"
    try:
        if filename[-3:]=='.gz':
            # filename is explicitly a gzip file
            file = gzip.open(filename,mode)
        else:
            # filename is given as a non-zip file
            file = open(filename,mode)
    except:
        # if filename is not existing upon read, then try append the '.gz' ending
        file = gzip.open(filename+'.gz',mode)
    return file

#--------------------------------------------------------------------------------
# XV-format

def ReadXVFile(filename,InUnits='Bohr',OutUnits='Ang',ReadVelocity=False):
    "Returns tuple (vectors,speciesnumber,atomnumber,xyz,[v,]) from an XV-file"
    print 'io.siesta.ReadXVFile: Reading',filename
    if (InUnits=='Bohr') and (OutUnits=='Ang'): convFactor = PC.Bohr2Ang
    elif (InUnits=='Ang') and (OutUnits=='Bohr'): convFactor = PC.Ang2Bohr
    elif (((InUnits=='Ang') and (OutUnits=='Ang')) \
       or ((InUnits=='Bohr') and (OutUnits=='Bohr'))): convFactor = 1
    else: print 'io.siesta.ReadXVFile: Unit conversion error!'
    file = SIO_open(filename,'r')
    # Read cell vectors (lines 1-3)
    vectors = []
    for i in range(3):
        data = string.split(file.readline())
        vectors.append([string.atof(data[j])*convFactor for j in range(3)])
    # Read number of atoms (line 4)
    numberOfAtoms = string.atoi(string.split(file.readline())[0])
    # Read remaining lines
    speciesnumber, atomnumber, xyz, V = [], [], [], []
    for line in file.readlines():
        if len(line)>5: # Ignore blank lines
            data = string.split(line)
            speciesnumber.append(string.atoi(data[0]))
            atomnumber.append(string.atoi(data[1]))
            xyz.append([string.atof(data[2+j])*convFactor for j in range(3)])
            V.append([string.atof(data[5+j])*convFactor for j in range(3)])
    file.close()
    if len(speciesnumber)!=numberOfAtoms:
        print 'io.siesta.ReadXVFile: Inconstency in %s detected!' %filename
    if ReadVelocity:
        return N.array(vectors),N.array(speciesnumber),N.array(atomnumber),N.array(xyz), N.array(V)
    else:
        return N.array(vectors),N.array(speciesnumber),N.array(atomnumber),N.array(xyz)

def WriteXVFile(filename,vectors,speciesnumber,atomnumber,xyz,\
                InUnits='Ang',OutUnits='Bohr',Velocity=[]):
    "Writes (vectors,speciesnumber,atomnumber,xyz,[V],) to the XV-file format"
    print 'io.siesta.WriteXVFile: Writing',filename
    if (InUnits=='Bohr') and (OutUnits=='Ang'): convFactor = PC.Bohr2Ang
    elif (InUnits=='Ang') and (OutUnits=='Bohr'): convFactor = PC.Ang2Bohr
    elif (((InUnits=='Ang') and (OutUnits=='Ang')) \
       or ((InUnits=='Bohr') and (OutUnits=='Bohr'))): convFactor = 1
    else: print 'io.siesta.WriteXVFile: Unit conversion error!'
    file = SIO_open(filename,'w')
    # Write basis vectors (lines 1-3)
    for i in range(3):
        line = '  '
        for j in range(3):
            line += string.rjust('%.9f'%(vectors[i][j]*convFactor),16)
        line += '     0.00  0.00  0.00 \n'
        file.write(line)
    # Write number of atoms (line 4)
    numberOfAtoms = len(speciesnumber)
    file.write('     %i\n' %numberOfAtoms)
    # Go through the remaining lines
    for i in range(numberOfAtoms):
        line = '  %i' %speciesnumber[i]
        line += '  %i  ' %atomnumber[i]
        for j in range(3):
            line += string.rjust('%.9f'%(xyz[i][j]*convFactor),16)
        if len(Velocity)==0:
            line += '     0.00  0.00  0.00 '
        else:
            for j in range(3):
                line += string.rjust('%.9f'%(Velocity[i][j]*convFactor),16)
        file.write(line+'\n')
    file.close()

def ReadAXVFile(filename,MDstep,tmpXVfile="tmp.XV",ReadVelocity=False):
    "Read concatenated XV files from an MD simulation"    
    # Determine the number of atoms in the supercell
    f = SIO_open(filename,'r')
    [f.readline() for i in range(3)]
    Natoms = int(string.split(f.readline())[0])
    f.close()
    XVlines = Natoms + 4
    # Extract the XV file from MDstep and write tmpXVfile
    f = SIO_open(filename,'r') 
    g = open(tmpXVfile,'w')
    i = 0
    for line in f:
        if i >= MDstep*XVlines and i < (MDstep+1)*XVlines: g.write(line)
        if i == (MDstep+1)*XVlines: break
        i += 1
    f.close()
    g.close()
    # Read tmpXVfile corresponding to MDstep
    return ReadXVFile(tmpXVfile,InUnits='Bohr',OutUnits='Ang',ReadVelocity=ReadVelocity)
    


#--------------------------------------------------------------------------------
# AXSF-format
def WriteAXSFFiles(filename,geoms,forces=None):
    'Writes [geom1, geom2 ...] to AXSF format'
    f = open(filename,'w')
    f.write('ANIMSTEPS %i\nCRYSTAL\n'%len(geoms))    
    for i in range(len(geoms)):
        f.write('PRIMVEC %i\n'%(i+1))
        f.write('%.6f %.6f %.6f\n'%(geoms[i].pbc[0][0],geoms[i].pbc[0][1],geoms[i].pbc[0][2]))
        f.write('%.6f %.6f %.6f\n'%(geoms[i].pbc[1][0],geoms[i].pbc[1][1],geoms[i].pbc[1][2]))
        f.write('%.6f %.6f %.6f\n'%(geoms[i].pbc[2][0],geoms[i].pbc[2][1],geoms[i].pbc[2][2]))
        f.write('PRIMCOORD %i\n'%(i+1))
        f.write('%i 1\n'%(len(geoms[i].xyz)))
        
        for j in range(len(geoms[i].xyz)):
            ln = ' %i'%geoms[i].anr[j]
            for k in range(3):
                ln += ' %.6f'%geoms[i].xyz[j][k]
            if type(forces) != bool:
                for k in range(3):
                    ln += ' %.6f'%forces[i][j][k]
            ln += '\n'
            f.write(ln)
    f.close()


# ANI-format

def WriteANIFile(filename,Geom,Energy,InUnits='Ang',OutUnits='Ang'):
    " Write .ANI file from list of geometries with Energy=[E1,E2..]"
    print 'io.siesta.WriteANIFile: Writing',filename
    if (InUnits=='Bohr') and (OutUnits=='Ang'): convFactor = PC.Bohr2Ang
    elif (InUnits=='Ang') and (OutUnits=='Bohr'): convFactor = PC.Ang2Bohr
    elif (((InUnits=='Ang') and (OutUnits=='Ang')) \
       or ((InUnits=='Bohr') and (OutUnits=='Bohr'))): convFactor = 1
    else: print 'io.siesta.WriteANIFile: Unit conversion error!'
    file = open(filename,'w')
    for ii, iGeom in enumerate(Geom):
        file.write('%i \n'%iGeom.natoms)
        file.write('%f \n'%Energy[ii])
        for iixyz in range(iGeom.natoms):
            file.write('%s %2.6f %2.6f %2.6f\n'%\
                       (PC.PeriodicTable[abs(iGeom.anr[iixyz])],\
                        convFactor*iGeom.xyz[iixyz][0],\
                        convFactor*iGeom.xyz[iixyz][1],\
                        convFactor*iGeom.xyz[iixyz][2]))
    file.close()

def ReadANIFile(filename,InUnits='Ang',OutUnits='Ang'):
    "Returns tuple (Geometry,Energy[Ry?],) from an ANI-file"
    import Inelastica.MakeGeom as MG
    print 'io.siesta.ReadANIFile: Reading',filename
    if (InUnits=='Bohr') and (OutUnits=='Ang'): convFactor = PC.Bohr2Ang
    elif (InUnits=='Ang') and (OutUnits=='Bohr'): convFactor = PC.Ang2Bohr
    elif (((InUnits=='Ang') and (OutUnits=='Ang')) \
       or ((InUnits=='Bohr') and (OutUnits=='Bohr'))): convFactor = 1
    else: print 'io.siesta.ReadANIFile: Unit conversion error!'

    Energy, Geom = [], []
    file = SIO_open(filename,'r')
    newNN=file.readline()
    while newNN!='':
        NN=string.atoi(newNN)
        newG = MG.Geom()
        try:
            Energy.append(string.atof(file.readline()))
        except:
            Energy.append(0.0)
        for ii in range(NN):
            line=string.split(file.readline())
            xyz=[convFactor*string.atof(line[1]),
                 convFactor*string.atof(line[2]),\
                 convFactor*string.atof(line[3])]
            newG.addAtom(xyz,1,PC.PeriodicTable[line[0]])
        Geom.append(newG)
        newNN=file.readline()
    file.close()
    return Geom,Energy 

#--------------------------------------------------------------------------------
# Reading SIESTA FC ascii files

def ReadFCFile(filename):
    "Returns FC from an FC-file"
    print 'io.siesta.ReadFCFile: Reading',filename
    file = SIO_open(filename,'rb')
    # Read comment line (line 1)
    line = file.readline()
    if string.strip(line)!='Force constants matrix':
        print 'io.siesta.ReadFCFile: Inconstency in %s detected!' %filename
    # Read remaining lines
    FC = []
    for line in file.readlines():
        data = string.split(line)
        FC.append([string.atof(data[j]) for j in range(3)])
    file.close()
    return FC

#--------------------------------------------------------------------------------
# Reading SIESTA Fortan binary files

def ReadFortranBin(file,type,num,printLength=False,unpack=True):
    "Reads Fortran binary data structures"
    import struct
    fmt = ''
    for i in range(num): fmt += type
    bin = file.read(struct.calcsize(fortranPrefix+fortranuLong+fmt+fortranuLong))
    if unpack:
        data = struct.unpack(fortranPrefix+fortranLong+fmt+fortranLong,bin)
        if printLength:
            print 'io.siesta.ReadFortranBin: %i bytes read' %data[0]
        if data[0]!=data[-1] or data[0]!=struct.calcsize(fortranPrefix+fmt):
            print 'io.siesta.ReadFortranBin: Error reading Fortran formatted binary file'
            sys.exit(1)
        return data[1:-1]
    else:
        return bin
    
def ReadWFSFile(filename):
    """
    TF/071008
    Returns the WF coefficients etc. from SIESTA systemlabe.WFS files
    (see siesta/utils/readwf.f for details)
    """
    file = SIO_open(filename,'rb')
    nk, = ReadFortranBin(file,'I',1)
    nspin, = ReadFortranBin(file,'I',1)
    nuotot, = ReadFortranBin(file,'I',1)
    #print nk, nspin, nuotot

    PSIvectors = []
    for iik in range(nk):
        for iispin in range(nspin):
            ik,k1,k2,k3 =  ReadFortranBin(file,'Iddd',1)
            #print ik,k1,k2,k3
            ispin, = ReadFortranBin(file,'I',1)
            nwflist, = ReadFortranBin(file,'I',1)
            #print ispin,nwflist
            for iw in range(nwflist):
                REpsi,IMpsi = [],[]
                indwf, = ReadFortranBin(file,'I',1)
                energy, = ReadFortranBin(file,'d',1)
                #print indwf, energy
                for jj in range(nuotot):
                    label = ''
                    for a in range(20): label += 'c'
                    out =  ReadFortranBin(file,'I'+label+'III'+label+'dd',1)
                    iaorb,lab1,j,iphorb,cnfigfio,lab2,repsi,impsi = out[0],out[1:21],out[21],out[22],out[23],out[24:44],out[44],out[45]
                    labelfis,symfio = '',''
                    for a in range(20):
                        labelfis += lab1[a]
                        symfio += lab2[a]
                    #print labelfis,symfio
                    REpsi.append(repsi),IMpsi.append(impsi)
                PSIvectors.append(N.array(REpsi)+1j*N.array(IMpsi))
    file.close()
    return nk,nspin,nuotot,nwflist,PSIvectors


def printDone(i,n,mess):
    # Print progress report
    if n>10:
        if i%(int(n/10)+1)==0:
            print mess,": %3.0f %% done" %(10.0*int(10.0*float(i+1)/n))
        if i+1==n:
            print mess,": 100 % done"
        sys.stdout.flush()

#--------------------------------------------------------------------------------
# MKL-format IO

def WriteMKLFile(filename,atomnumber,xyz,freq,vec,FCfirst,FClast):
    "Writes a MKL-file"
    print 'io.siesta.WriteMKLFile: Writing',filename
    file = open(filename,'w')
    file.write('$MKL\n$COORD\n')
    for i in range(len(atomnumber)):
        line = str(atomnumber[i])
        for j in range(3):
            line += string.rjust('%.9f'%xyz[i][j],16)
        line +='\n'
        file.write(line)
    file.write('$END\n')
    if len(freq)>0:
        file.write('$FREQ\n')
        for i in range(len(freq)/3):
            file.write('C1 C1 C1\n')
            # Write frequencies
            line = ''
            for j in range(3):
                f = 1000*freq[3*i+j] # Write in meV
                try:
                    line += '%f '%f.real
                except:
                    line += '%f '%f
            line += '\n'
            file.write(line)
            # Write modes
            for j in range(FCfirst-1):
                file.write('0 0 0 0 0 0 0 0 0\n')
            for j in range(FClast-FCfirst+1):
                line = ''
                for k in range(3):
                    line += '%f %f %f '%(vec[3*i+k][3*j],vec[3*i+k][3*j+1],vec[3*i+k][3*j+2])
                line += '\n'
                file.write(line)
            for j in range(FClast,len(xyz)):
                file.write('0 0 0 0 0 0 0 0 0\n')
        file.write('$END\n\n')
    file.close()

#--------------------------------------------------------------------------------
# XYZ-format IO

def ReadXYZFile(filename):
    file = SIO_open(filename,'r')
    # Read number of atoms (line 4)
    numberOfAtoms = string.atoi(string.split(file.readline())[0])
    # Read remaining lines
    label, atomnumber, xyz = [], [], []
    for line in file.readlines():
        if len(line)>5: # Ignore blank lines
            data = string.split(line)
            label.append(data[0])
            atomnumber.append(PC.PeriodicTable[data[0]])
            xyz.append([string.atof(data[1+j]) for j in range(3)])
    file.close()
    if len(xyz)!=numberOfAtoms:
        print 'io.siesta.ReadXYZFile: Inconstency in %s detected!' %filename
    return label,N.array(atomnumber),N.array(xyz)


def WriteXYZFile(filename,atomnumber,xyz,write_ghosts=False):
    "Writes atomic geometry in xyz-file format"
    print 'io.siesta.WriteXYZFile: Writing',filename
    # Number of ghost atoms
    nga = len(N.where(N.array(atomnumber)<0)[0])
    # Write file
    file = open(filename,'w')
    if write_ghosts:
        file.write(str(len(xyz)))
    else:
        if nga>0:
            print '... skipped %i ghost atoms'%nga
        file.write(str(len(xyz)-nga))
    file.write('\n\n')
    for i in range(len(xyz)):
        try:
            element = PC.PeriodicTable[abs(atomnumber[i])]
        except:
            element = 'X' 
        line = string.ljust(element,5)
        for j in range(3):
            line += string.rjust('%.9f'%xyz[i][j],16)
        line +='\n'
        if atomnumber[i]>0 or write_ghosts:
            file.write(line)
    file.close()

#--------------------------------------------------------------------------------
# FDF format IO
    
def ReadFDFFile(infile):
    """ Reads an FDF file and gives the output values: pbc, xyz, snr, anr, natoms
        infile = FDF inputfile"""
    pbc = Getpbc(infile)
    xyz = Getxyz(infile,pbc)
    snr = Getsnr(infile)
    anr = Getanr(infile)
    natoms = Getnatoms(infile)
    if natoms != len(xyz):
        print 'Error! natoms != len(xyz)'
    return pbc, xyz, snr, anr, natoms

def WriteFDFFile(filename,vectors,speciesnumber,atomnumber,xyz):
    "Write STRUCT.fdf file"
    print 'io.siesta.WriteFDFFile: Writing',filename
    file = open(filename,'w')
    file.write('NumberOfAtoms '+str(len(xyz))+'\n')
    file.write('NumberOfSpecies '+str(max(speciesnumber))+'\n')
    file.write('LatticeConstant 1.0 Ang\n%block LatticeVectors\n')
    for ii in range(3):
        for jj in range(3):
            file.write(string.rjust('%.9f'%vectors[ii][jj],16)+' ')
        file.write('\n')
    file.write('%endblock LatticeVectors\nAtomicCoordinatesFormat  Ang'+
               '\n%block AtomicCoordinatesAndAtomicSpecies\n')
    for ii in range(len(xyz)):
        line=string.rjust('%.9f'%xyz[ii][0],16)+' '
        line+=string.rjust('%.9f'%xyz[ii][1],16)+' '
        line+=string.rjust('%.9f'%xyz[ii][2],16)+' '
        line+=str(int(speciesnumber[ii]))+' # %i\n'%(ii+1)
        file.write(line)
    file.write('%endblock AtomicCoordinatesAndAtomicSpecies\n')

def WriteFDFFileZmat(filename,vectors,speciesnumber,atomnumber,xyz,first=0,last=0,zmat=None):
    """Write STRUCT.fdf file using the z-matrix format
       xyz        : all cartesian coordinates
       first/last : first,last atoms in molecule block
       zmat[:,:3] : integer part of z-matrix (indices)
       zmat[:,3:] : fractional part of z-matrix (angles)"""
    # Sanity check
    if first > last or first > len(xyz) or first < 0 or last > len(xyz) or last < 0:
        print 'io.siesta.WriteFDFFileZmat: Meaningless first (%i) / last (%i) inputs. '%(first,last)
        first, last = 0,0
    # Writing zmatrix
    print 'io.siesta.WriteFDFFileZmat: Writing',filename
    file = open(filename,'w')
    file.write('NumberOfAtoms '+str(len(xyz))+'\n')
    file.write('NumberOfSpecies '+str(max(speciesnumber))+'\n')
    file.write('LatticeConstant 1.0 Ang\n%block LatticeVectors\n')
    for ii in range(3):
        for jj in range(3):
            file.write(string.rjust('%.9f'%vectors[ii][jj],16)+' ')
        file.write('\n')
    file.write('%endblock LatticeVectors\nAtomicCoordinatesFormat Ang'+
               '\n\nZM.UnitsLength Ang\nZM.UnitsAngle deg\n'+
               '\n%block Zmatrix\n')
    if first != 1:
        file.write('cartesian\n')
    for ii in range(len(xyz)):
        if ii+1 == first:
            file.write('molecule\n')
        if ii+1 >= first and ii+1 <= last:
            # We are within the molecular block
            line =string.rjust('%i'%speciesnumber[ii],2)
            a,b,c,d,e,f = zmat[ii+1-first]
            line+=' %i %i %i '%(a,b,c)
            line+=string.rjust('%.9f'%d,16)
            line+=string.rjust('%.9f'%e,16)
            line+=string.rjust('%.9f'%f,16)
            line+='   0 0 0\n'
            file.write(line)
        else:
            line =string.rjust('%i'%speciesnumber[ii],2)
            line+=string.rjust('%.9f'%xyz[ii][0],16)
            line+=string.rjust('%.9f'%xyz[ii][1],16)
            line+=string.rjust('%.9f'%xyz[ii][2],16)
            line+='   0 0 0\n'
            file.write(line)
        if ii+1 == last:
            file.write('cartesian\n')
    file.write('constants\n')
    file.write('variables\n')
    file.write('constraints\n')
    file.write('%endblock Zmatrix\n')

#--------------------------------------------------------------------------------                                                                                                                                                           
# Read systemlabel.STRUCT_OUT files

def ReadSTRUCT_OUTFile(filename):
    file = SIO_open(filename,'r')
    # Read cell vectors (lines 1-3)
    vectors = []
    for i in range(3):
        data = string.split(file.readline())
        vectors.append([string.atof(data[j]) for j in range(3)])
    # Read number of atoms (line 4)
    numberOfAtoms = string.atoi(string.split(file.readline())[0])
    # Read remaining lines
    speciesnumber, atomnumber, xyz = [], [], []
    for line in file.readlines():
        if len(line)>4: # Ignore blank lines
            data = string.split(line)
            speciesnumber.append(string.atoi(data[0]))
            atomnumber.append(string.atoi(data[1]))
            xyz.append([string.atof(data[2+j]) for j in range(3)])
    file.close()
    if len(speciesnumber)!=numberOfAtoms:
        print 'io.siesta.ReadSTRUCT_OUTFile: Inconstency in %s detected!' %filename
    xyz = N.array(xyz,N.float)
    vectors = N.array(vectors,N.float)
    for i in range(len(xyz)):
        xyz[i] = N.dot(vectors,xyz[i])
    return vectors,speciesnumber,atomnumber,xyz
          
#--------------------------------------------------------------------------------
# Writing SIESTA Fortran binary files

def WriteFortranBin(file,type,data):
    "Writes Fortran binary data structures"
    import struct
    try:
        L = len(data)
        if L == 1: data = data[0]
    except: L = 1
    fmt = ''
    for i in range(L): fmt += type
    bin = struct.pack(fortranPrefix+fortranLong,struct.calcsize(fmt))
    if L>1:
        for i in range(L):
            bin += struct.pack(type,data[i])
    else:
        bin += struct.pack(type,data)
    bin += struct.pack(fortranuLong,struct.calcsize(fmt))
    file.write(bin)


#--------------------------------------------------------------------------------
# "Low-level" FDF format functionality

def ReadFDFLines(infile,head='', printAlot=True):
    """ Returns an FDF file and all the %include files as split strings
        infile = input file"""
    infile = os.path.abspath(infile)
    if head == '':
        head,tail =  os.path.split(infile)
    if printAlot:
        print 'io.siesta.ReadFDFLines: Reading', infile
    file = SIO_open(infile,'r')
    lines = []
    tmp = file.readline()
    while tmp != '':
        if len(tmp)>3:
            tmp = tmp.replace(':',' ') # Remove ':' from fdf
            tmp = tmp.replace('=',' ') # Remove '=' from fdf
            tmp = string.split(tmp)
            for ii,s in enumerate(tmp):  # Remove comments
                if s[0]=="#":
                    break
            if s[0]=='#':
                tmp = tmp[0:ii]
            if len(tmp)>0:
                if tmp[0] == '%include':
                    subfile = head+'/'+tmp[1]
                    tmp2 = ReadFDFLines(subfile,head=head, printAlot=printAlot)
                    lines += tmp2
                else:
                    lines.append(tmp)
        tmp = file.readline()
    file.close()
    return lines


def Getnatoms(infile):
      """ Gives the number of atoms included in an FDF file
          infile = FDF input file"""
      natoms = GetFDFline(infile,'NumberOfAtoms')
      if natoms==None:
          natoms = GetFDFline(infile,'NumberOfAtoms:')
      return int(natoms[0])

   
def Getxyz(infile,pbc=[]):
    """ Gives a list of the xyz posistions in a FDF file
        infile = FDF input file"""
    latt_const = GetFDFline(infile,'LatticeConstant')
    AC_format=GetFDFline(infile,'AtomicCoordinatesFormat')
    if latt_const[1]=='Bohr' :
	latt_const=float(latt_const[0])*PC.Bohr2Ang
    else :
     	latt_const=float(latt_const[0])
    data = GetFDFblock(infile, 'AtomicCoordinatesAndAtomicSpecies')
    xyz = []
    if AC_format[0]=='Ang' or AC_format[0]=='NotScaledCartesianAng' :
    	for i in range(len(data)):
    	    xyz.append([string.atof(data[i][j]) for j in range(3)])
    elif AC_format[0]=='Bohr' or AC_format[0]=='NotScaledCartesianBohr' :
    	for i in range(len(data)):
    	    xyz.append([string.atof(data[i][j])*PC.Bohr2Ang for j in range(3)])
    elif AC_format[0]=='ScaledCartesian':    
    	for i in range(len(data)):
    	    xyz.append([string.atof(data[i][j])*latt_const for j in range(3)])
    elif AC_format[0]=='Fractional' or AC_format[0]=='ScaledByLatticeVectors' :
	for i in range(len(data)):
    	    xyz.append([string.atof(data[i][0])*pbc[0][j]+string.atof(data[i][1])*pbc[1][j]+string.atof(data[i][2])*pbc[2][j] for j in range(3)])
    else :
	print("Give correct AtomicCoordinates Format")
	sys.exit(1)
    return xyz


def Getpbc(infile):
    """ Gives a list of the lattice vectores in a FDF file
        infile = FDF input file"""
    data = GetFDFblock(infile,'LatticeVectors')
    pbc = []
    latt_const = GetFDFline(infile,'LatticeConstant')
    if latt_const[1]=='Bohr' :
	latt_const=float(latt_const[0])*PC.Bohr2Ang
    else :
     	latt_const=float(latt_const[0])
    for i in range(len(data)):
        pbc.append([string.atof(data[i][j])*latt_const for j in range(3)])
    return pbc


def Getsnr(infile):
     """ Gives a list of the species numbers in a FDF file
         infile = FDF input file"""
     data = GetFDFblock(infile,'AtomicCoordinatesAndAtomicSpecies')
     snr = []
     for i in range(len(data)):
         snr.append(string.atoi(data[i][3]))
     return snr


def Getanr(infile):
    """ Gives a list of the atomic numbers in a FDF file
        infile = FDF input file"""   
    data = GetFDFblock(infile,'ChemicalSpeciesLabel')
    tmp = []
    table ={} 
    for i in range(len(data)):
        tmp.append([string.atoi(data[i][j]) for j in range(2)])
        table[tmp[i][0]] = tmp[i][1]
    snr = Getsnr(infile)
    anr = []
    for i in range(len(snr)):
        anr.append(table[snr[i]])
    return anr
        
   
def GetFDFline(infile, KeyWord = '', printAlot=True):
    """ Finds a line and gives the value as a string
        infile = FDF input file
        KeyWord = line to find"""
    lines = ReadFDFLines(infile, printAlot=printAlot)
    kwl = KeyWord.lower()
    for s in ['-','_','.']: # these characters should be ignored in fdf keys
        kwl = kwl.replace(s,'')
    for line in lines:
        key = line[0].lower()
        for s in ['-','_','.']: # these characters should be ignored in fdf keys
            key = key.replace(s,'')
        if key == kwl:
            return line[1:]

def GetFDFlineWithDefault(infile, key, type, default, error):
    """ Finds a line and gives the value of type type.
        If not found, default returned.
        If default=None, print error and exit.
    """
    data = GetFDFline(infile, key, printAlot=False)
    if data is None:
        if default is None:
            raise LookupError('GetFDFlineWithDefault failed to find "{}" in file {}.'.format(key,infile))
        else:
            return default
    else:
        data = data[0]
        # Boolean is tricky!
        if type != bool:
            return type(data)
        else:
            data = data.lower()
            if data in ['true','t','.true.','yes','y']:
                return True
            elif data in ['false','f','.false.','no','n']:
                return False
            else:
                raise TypeError('GetFDFlineWithDefault failed to convert '+
                                '"{}" to boolean from key "{}" in file {}.'.format(data,key,infile))
            
def GetFDFblock(infile, KeyWord = ''):
    """Finds the values in a block as strings
       infile = FDF input file
       KeyWord = block to find"""
    lines = ReadFDFLines(infile)
    kwl = KeyWord.lower()
    data = []
    start = 0
    for i,line in enumerate(lines):
        if line[0].lower() == '%block':
            for s in ['-','_','.']: # these characters should be ignored in fdf keys
                line[1] = line[1].replace(s,'')
                kwl = kwl.replace(s,'')
            if line[1].lower() == kwl:
                start = i+1
                break
    if start > 0: # Only append data if block was found
        for i in range(len(lines)):
            tmp = lines[i+start]
            if tmp[0].lower() != '%endblock':
                data.append(tmp)
            else: break
    return data

#--------------------------------------------------------------------------------

def GetTotalEnergy(infile):
    # Find total energy from SIESTA stdout file
    f = SIO_open(infile,'rb')
    lines = f.readlines()
    f.close()
    E = 0.0
    for line in lines:
        words = string.split(line)
        if 'Total' in words and '=' in words:
            E = string.atof(words[-1])
            break
    return E

def GetFermiEnergy(infile):
    # Read Fermi energy from SIESTA stdout file
    print 'io.siesta.GetFermiEnergy: Reading',infile
    f = SIO_open(infile,'rb')
    lines = f.readlines()
    f.close()
    E = 0.0
    for line in lines:
        words = string.split(line)
        if 'Fermi' in words and 'energy' in words:
            E = string.atof(words[-2])
            print '... eF = %.4f eV'%E
            break
    if E==0.0:
        print '... did not find the Fermi energy'
    return E

def ReadEIGfile(infile,printing=False,FermiRef=True):
    # Read *EIG file and print eigenvalues with respect to eF.
    f = SIO_open(infile,'rb')
    eF = float(string.split(f.readline())[0])
    f.readline() # Skip second line
    EIG = []
    eig = string.split(f.readline())[1:] # Skip first entry in third line
    EIG += eig
    for line in f.readlines():
        eig = string.split(line)
        EIG += eig
    if not FermiRef:
        # Do not take Fermi level as energy reference
        eF = 0.0
    if printing:
        print '# State, eigenvalue wrt. eF'
        for i in range(len(EIG)):
            print i+1, float(EIG[i])-eF
    f.close()
    return eF

def ReadMullikenPop(infile,outfile,writeallblocks=False):
    # Read Mulliken populations from the *.out file
    mline = False
    mpop = []
    popsum = 0.0
    iter = 0
    block = 0
    spin = False
    f = SIO_open(infile,'r')
    print 'io.siesta.ReadMullikenPop: Reading',infile
    for line in f.readlines():
        if 'mulliken: Atomic and Orbital Populations:' in line:
            # Start of populations block
            mline = True
            iter += 1
        if 'mulliken: Spin' in line:
            spin = True
        if 'mulliken: Qtot =' in line:
            # Determine whether or not we are still within a populations block
            if spin and block%2==0: mline = True
            else: mline = False
	    mpop.sort()
            # Write data to file
            if writeallblocks:
                thisfile = outfile+'%.2i'%iter
            else:
                thisfile = outfile
            if spin and block%2==0:
                thisfile += '.UP'
            if spin and block%2==1:
                thisfile += '.DOWN'
            print 'io.siesta.ReadMullikenPop: Writing',thisfile
            f2 = open(thisfile,'w')
            f2.write('# Sum of Mulliken charges: %.6f\n'%popsum)
            f2.write('# Atomnr  Pop.  dPop   Cum.sum.\n')
            for i in range(len(mpop)):
                f2.write('  %.i  %.6f  %.6f  %.6f\n'%mpop[i])
            f2.close()
            block += 1
            mpop = []
            popsum = 0.0
        if mline:
            # We are inside a range of lines with mpop
            s = string.split(line)
            try:
                nr,pop = int(s[0]),float(s[1])
                dpop = pop-round(pop,0)
                popsum += pop
                mpop.append((nr,pop,dpop,1.0*popsum))
            except:
                pass
    f.close()

def ReadForces(infile):
    # Read forces from the *.out file
    data = []
    fline = False
    f = SIO_open(infile,'r')
    print 'io.siesta.ReadForces: Reading',infile
    for line in f.readlines():
        if fline:
            # We are inside a range of lines with forces
            try:
                s = string.split(line)
	        data.append([int(s[0]),float(s[1]),float(s[2]),float(s[3])])
            except:
                fline = False
        if 'Atomic forces' in line:
            # Start of forces block
            fline = True
    f.close()
    return data

def ReadFAFile(filename):
    "Returns forces from a FA-file"
    print 'io.siesta.ReadFAFile: Reading',filename
    file = SIO_open(filename,'rb')
    # Read comment line (line 1)
    line = file.readline()
    natoms = int(string.strip(line))
    # Read remaining lines
    FA = []
    for line in file.readlines():
        data = string.split(line)
        FA.append([string.atof(data[j]) for j in range(1,4)])
    file.close()
    return N.array(FA)


def ReadTRANSAVfile(infile):
    # Read (averaged *.TRANS.AV) transmission function
    data = []
    e = -1e10
    f = SIO_open(infile,'r')
    for line in f.readlines():
        s = string.split(line)
        if float(s[0])>e:
            data.append([float(s[0]),float(s[1])])
            e = float(s[0])
    data = N.array(data)
    [e,t] = N.transpose(data)
    f.close()
    return e,t
    
#--------------------------------------------------------------------------------
# Related to band structure calculations and density of states

def CrossProd(A,B):
    "Returns the cross product of two geometric vectors"
    [ax,ay,az] = A
    [bx,by,bz] = B
    return N.array([ay*bz-az*by,az*bx-ax*bz,ax*by-ay*bx])

def GetReciprocalLatticeVectors(infile):
    # Calculate reciprocal lattice vectors
    print 'io.siesta.GetReciprocalLatticeVectors: Reading',infile
    if infile.endswith('.fdf'):
        pbc = Getpbc(infile)
    elif infile.endswith('.XV'):
        pbc,speciesnumber,atomnumber,xyz = ReadXVFile(infile)
    a0 = N.array(pbc[0])
    a1 = N.array(pbc[1])
    a2 = N.array(pbc[2])
    b0 = 2*N.pi*CrossProd(a1,a2)/(N.dot(a0,CrossProd(a1,a2)))
    b1 = 2*N.pi*CrossProd(a2,a0)/(N.dot(a1,CrossProd(a2,a0)))
    b2 = 2*N.pi*CrossProd(a0,a1)/(N.dot(a2,CrossProd(a0,a1)))
    return b0,b1,b2

def WriteBandLinesFDF(infile,outfile,res=0.05):
    # Sets up an fdf-file which can be included in RUN.fdf for
    # calculating band structure and total density of states
    print 'io.siesta.WriteBandLinesFDF: Writing fdf-input file',outfile
    b0,b1,b2 = GetReciprocalLatticeVectors(infile)
    f = open(outfile,'w')
    f.write('BandLineScale  pi/a\n%block BandLines\n')
    f.write('1   0.00  0.00  0.00  \Gamma\n')
    n = int(N.dot(b0,b0)**.5/res)
    f.write('%i'%n+'  %.8f  %.8f  %.8f  b0\n'%tuple(b0/(2*N.pi)))
    n = int(N.dot(b1,b1)**.5/res)
    f.write('%i'%n+'  %.8f  %.8f  %.8f  b1\n'%tuple((b0+b1)/(2*N.pi)))
    n = int(N.dot(b0+b1,b0+b1)**.5/res)
    f.write('%i'%n+'  0.00  0.00  0.00  \Gamma\n')
    #n = int(N.dot(b2,b2)**.5/res)
    #f.write('%i'%n+'  %.8f  %.8f  %.8f  b2\n'%tuple(b2/(2*N.pi)))
    f.write('%endblock BandLines\nWriteBands true\n\n')
    f.write('%endblock BandLines\nWriteBands true\nPDOSBandbyBand true\n\n')
    f.write('%block ProjectedDensityOfStates\n')
    f.write('-10.0 5.0 0.10 500 eV\n')
    f.write('%endblock ProjectedDensityOfStates\n')
    f.close()
    
def ReadBlock(file,lines,type='float'):
    # Designed for reading, e.g., *.bands files from SIESTA
    data = ''
    for i in range(lines):
        data += file.readline()
    data = data.split()
    for i in range(len(data)):
        if type=='int':
            data[i] = int(data[i])
        elif type=='float':
            data[i] = float(data[i])
        elif type=='str':
            pass
    return data

def ReadBandsFile(filename,origformat=True):
    print 'io.siesta.ReadBandsFile: Reading',filename
    # Reads SIESTA *.bands files
    f = SIO_open(filename,'r')
    # line 1:
    eF = ReadBlock(f,1)[0]
    # lines 2-3:
    xmin,xmax = ReadBlock(f,1)
    ymin,ymax = ReadBlock(f,1)
    # lines 4:
    eval,spins,kpts = ReadBlock(f,1,type='int')
    
    if origformat:
        linesPerKpt = (spins*eval)/10 # SIESTA
        if (spins*eval)%10 != 0: linesPerKpt += 1
    else:
        linesPerKpt = 1 # Siesta-2.5?
    # Read k-points
    EvsK = []
    for i in range(kpts):
        EvsK.append(ReadBlock(f,linesPerKpt))
    # Read labels
    numlabels = ReadBlock(f,1,type='int')[0]
    labels = []
    for i in range(numlabels):
        lab,val = ReadBlock(f,1,type='str')
        labels.append([float(lab),val[1:-1]])
    f.close()
    return eF,N.array(EvsK),labels,spins

def WriteBandsFile(filename,eF,EvsK,labels,spins):
    "Writes SIESTA *.bands files"
    print 'io.siesta.WriteBandsFile: Writing',filename
    xmin,xmax,ymin,ymax = min(EvsK[:,0]),max(EvsK[:,0]),min(EvsK[0,1:]),max(EvsK[0,1:])
    for i in range(1,len(EvsK)):
        if min(EvsK[i,1:])<ymin: ymin = min(EvsK[i,1:])
        if max(EvsK[i,1:])>ymax: ymax = max(EvsK[i,1:])
    f = open(filename,'w')
    f.write('   %.9f       \n'%eF)
    f.write('   %.9f  %.9f \n'%(xmin,xmax)) 
    f.write('   %.9f  %.9f \n'%(ymin,ymax))
    f.write('   %i    %i    %i\n'%(len(EvsK[0])-1,spins,len(EvsK)))
    for i in range(len(EvsK)):
        f.write('%.6f '%EvsK[i,0])
        for j in range(1,len(EvsK[0])):
            f.write('%.6f '%EvsK[i,j])
            if (j)%10==0: f.write('\n            ')
        f.write('\n')
    f.write('  %i \n'%len(labels))
    for l in labels:
        f.write(' %.6f  %s \n'%(l[0],l[1]))

def ConvertBandsFile(filename):
    eF,EvsK,labels,spins = ReadBandsFile(filename)
    print 'io.siesta.ConvertBandsFile: Writing',filename+'.dat'
    f = open(filename+'.dat','w')
    for j in range(1,len(EvsK[0])):
        for i in range(len(EvsK)):
            f.write('%.6f   %.6f\n'%(EvsK[i][0],EvsK[i][j]-eF))
        f.write('\n')
    f.close()

def ReadDOSFile(filename):
    # Reads SIESTA *.DOS files
    f = SIO_open(filename,'r')
    print 'io.siesta.ReadDOSFile: Reading',filename
    DOS = []
    for line in f.readlines():
        data = line.split()
        for i in range(len(data)):
            data[i] = float(data[i])
        DOS.append(data)
    f.close()
    return N.array(DOS)


#--------------------------------------------------------------------------------
# XML-format

import xml.dom.minidom as xml

def GetXMLFermiEnergy(dom):
    "Looks for the Fermi energy in xml file"
    try:
        node = dom.getElementsByTagName('E_Fermi')[0] # First (and only) entry
        eF = float(node.childNodes[0].data)
        print 'io.siesta.GetXMLFermiEnergy: Found E_Fermi = %.8f eV'%eF
    except:
        eF = 0.0
    return eF

def GetPDOSnspin(dom):
    "Returns an integer for the number of spins (variable nspin)"
    node = dom.getElementsByTagName('nspin')[0] # First (and only) entry
    return int(node.childNodes[0].data)

def GetPDOSnorbitals(dom):
    # Read norbitals
    node = dom.getElementsByTagName('norbitals')[0] # First (and only) entry
    return int(node.childNodes[0].data)

def GetPDOSenergyValues(dom):
    # Read energy values
    node = dom.getElementsByTagName('energy_values')[0] # First (and only) entry
    data = node.childNodes[0].data.split()
    for i in range(len(data)):
        data[i] = float(data[i])
    return N.array(data)

def GetPDOSfromOrbitals(dom,index=[],atom_index=[],species=[],nlist=[],llist=[],mlist=[]):
    dim = len(GetPDOSenergyValues(dom))*GetPDOSnspin(dom)
    pdos = N.zeros(dim,N.float)
    nodes = dom.getElementsByTagName('orbital')
    usedOrbitals = []
    for node in nodes:
        ok = True
        i = int(node.attributes['index'].value)
        ai = int(node.attributes['atom_index'].value)
        s = node.attributes['species'].value
        n = int(node.attributes['n'].value)
        l = int(node.attributes['l'].value)
        m = int(node.attributes['m'].value)
        if i not in index and index!=[]:
            ok = False
        if ai not in atom_index and atom_index!=[]: 
            ok = False
        if s not in species and species!=[]:
            ok = False
        if n not in nlist and nlist!=[]:
            ok = False
        if l not in llist and llist!=[]:
            ok = False
        if m not in mlist and mlist!=[]:
            ok = False
        if ok:
            usedOrbitals.append([i,ai,s,n,l])
            data = node.getElementsByTagName('data')[0] # First (and only) entry
            data = data.childNodes[0].data.split()
            for i in range(len(data)):
                data[i] = float(data[i])
            pdos += N.array(data)
    # Generate some output-related information
    if atom_index!=[]: print '... Atom indices =',atom_index
    if species!=[]: print '... Species =',species
    if nlist!=[]: print '... Allowed n quantum numbers =',nlist
    if llist!=[]: print '... Allowed l quantum numbers =',llist
    if mlist!=[]: print '... Allowed m quantum numbers =',mlist
    print '... Orbitals included =',len(usedOrbitals)
    usedAtoms = []
    for orb in usedOrbitals:
        if orb[1] not in usedAtoms: usedAtoms.append(orb[1])
    print '... Atoms included =',len(usedAtoms)
    return pdos,usedOrbitals,usedAtoms

def ReadPDOSFile(filename,index=[],atom_index=[],species=[],nlist=[],llist=[],mlist=[]):
    # Reads SIESTA *.PDOS files summing up contributions from orbitals
    # belonging to a subset specified by the keywords    
    file = SIO_open(filename,mode='r')
    dom = xml.parse(file)
    nspin = GetPDOSnspin(dom)
    norb = GetPDOSnorbitals(dom)
    ev = GetPDOSenergyValues(dom)
    eF = GetXMLFermiEnergy(dom)
    pdos,usedOrbitals,usedAtoms = GetPDOSfromOrbitals(dom,index,atom_index,species,nlist,llist,mlist)
    file.close()
    return nspin,norb,ev,pdos,usedOrbitals,usedAtoms,eF

def ExtractPDOS(filename,outfile,index=[],atom_index=[],species=[],nlist=[],llist=[],mlist=[],FermiRef=True,Normalize=False):
    print 'io.siesta.ExtractPDOS: Reading',filename 
    head,tail =  os.path.split(filename)
    nspin,norb,ev,pdos,usedOrbitals,usedAtoms,eF = ReadPDOSFile(filename,index,atom_index,species,nlist,llist,mlist)
    if FermiRef:
        # Set energy reference to the Fermi energy
        if eF == 0.0:
            eF = GetFermiEnergy(head+'/RUN.out')
        if eF == 0.0:
            print 'io.siesta.ExtractPDOS: Reading',filename[:-5]+'.EIG'
            eF = ReadEIGfile(filename[:-5]+'.EIG')
            print '... eF = %.4f eV'%eF
    else:
        # Set energy reference to SIESTAs internal
        eF = 0.0
    if Normalize and len(usedAtoms)>0:
        pdos = pdos/len(usedAtoms)
        print 'io.siesta.ExtractPDOS: Normalizing PDOS to states/atom/eV'
    if outfile!=None: # Write to file or return lists
        if nspin == 1: # No spin
            print 'io.siesta.ExtractPDOS: Writing', outfile
            f = open(outfile,'w')
            for i in range(len(ev)):
                f.write('%.6f %.9f\n'%(ev[i]-eF,pdos[i]))
            f.close()
        elif nspin == 2: # Spin polarized
            print 'io.siesta.ExtractPDOS: Writing', outfile
            f = open(outfile,'w')
            for i in range(len(ev)):
                f.write('%.6f %.9f %.9f\n'%(ev[i]-eF,pdos[2*i],-pdos[2*i+1]))
            f.close()
    else: 
        if nspin == 2:
            return nspin, ev-eF, [pdos[::2], pdos[1::2]]
        else:
            return nspin, ev-eF, [pdos]

# Functions specific to syslabel.PROJBANDS files
# (k-resolved PDOS is available from a modified SIESTA by D. Sanchez-Portal)

def GetPROJBANDSnbands(dom):
    """ Returns the (integer) number of bands """
    node = dom.getElementsByTagName('nbands')[0] # First (and only) entry
    return int(node.childNodes[0].data)

def GetPROJBANDSkpoint(dom):
    "Returns an array with the (3D) k-points"
    nodes = dom.getElementsByTagName('kpoint')
    kpts = []
    for node in nodes:
        data = node.childNodes[0].data.split()
        for i in range(3):
            data[i] = float(data[i])
        kpts.append(data)
    return N.array(kpts)

def GetPROJBANDSenergies(dom):
    """
    Returns an array with the band energies (in units of eV):
        energies = energies[ k-index, band-index, n-spin ]
    """
    nspin = GetPDOSnspin(dom)
    nbands = GetPROJBANDSnbands(dom)
    nodes = dom.getElementsByTagName('energies')
    energies = []
    for node in nodes:
        data = node.childNodes[0].data.split()
        for i in range(len(data)):
            data[i] = float(data[i])
        data = N.array(data)
        data = N.reshape(data,(nbands,nspin))
        energies.append(data)
    # NB: This conversion is weird and linked to the initial output of the modified SIESTA!!!
    return PC.Rydberg2eV**2*N.array(energies)

def GetPROJBANDSfromOrbitals(dom,index=[],atom_index=[],species=[],nlist=[],llist=[]):
    """
    Returns an array of the DOS projected onto bands, summing up contributions from
    orbitals as specified in the function call. The indexing is the following:
       pdos = pdos[ k-index, band-index, spin-index ]
    """
    nspin = GetPDOSnspin(dom)
    nbands = GetPROJBANDSnbands(dom)
    pdos = 0.0*GetPROJBANDSenergies(dom)
    nodes = dom.getElementsByTagName('orbital')
    k=-1
    for node in nodes:
        ok = True
        i = int(node.attributes['index'].value)
        if i==1: k+=1 # Next k-point
        ai = int(node.attributes['atom_index'].value)
        s = node.attributes['species'].value
        n = int(node.attributes['n'].value)
        l = int(node.attributes['l'].value)
        if i not in index and index!=[]:
            ok = False
        if ai not in atom_index and atom_index!=[]:
            ok = False
        if s not in species and species!=[]:
            ok = False
        if n not in nlist and nlist!=[]:
            ok = False
        if l not in llist and llist!=[]:
            ok = False
        if ok:
            print 'Adding PDOS from (k=%i,i=%i,ai=%i,s=%s,n=%i,l=%i)'%(k,i,ai,s,n,l)
            data = node.getElementsByTagName('bandproj')[0] # First (and only) entry
            data = data.childNodes[0].data.split()
            for i in range(len(data)):
                data[i] = float(data[i])
            data = N.array(data)
            data = N.reshape(data,(nbands,nspin))
            pdos[k] += N.array(data)
    return pdos

def ReadPROJBANDSfile(filename,index=[],atom_index=[],species=[],nlist=[],llist=[]):
    dom = xml.parse(filename)
    nspin = GetPDOSnspin(dom)
    nbands = GetPROJBANDSnbands(dom)
    norb = GetPDOSnorbitals(dom)
    energies = GetPROJBANDSenergies(dom)
    kresPDOS = GetPROJBANDSfromOrbitals(dom,index,atom_index,species,nlist,llist)
    return nspin,nbands,norb,energies,kresPDOS


def ExtractPROJBANDS(filename,outfile,index=[],atom_index=[],species=[],nlist=[],llist=[],emin=-5.0,emax=5.0):
    eF = GetFermiEnergy('RUN.out')
    dom = xml.parse(filename)
    nspin = GetPDOSnspin(dom)
    nbands = GetPROJBANDSnbands(dom)
    norb = GetPDOSnorbitals(dom)
    kpts = GetPROJBANDSkpoint(dom)
    energies = GetPROJBANDSenergies(dom)
    kresPDOS = GetPROJBANDSfromOrbitals(dom,index,atom_index,species,nlist,llist)
    if nspin == 1: # No spin
        points = 0
        f = open(outfile+'.dx.dat','w')
        for i in range(len(kpts)):
            for j in range(len(energies[i])):
                if energies[i][j]-eF >= emin and energies[i][j]-eF <= emax:
  	            f.write('%i %.6f %.9f\n'%(i,energies[i][j][0]-eF,kresPDOS[i][j][0]))
                    points += 1
        f.close()
    else:
        sys.exit('Not yet implemented for spin polarized data.\n')
    # Write DX datafile
    f = open(outfile+'.dx.general','w')
    f.write('file = %s.dx.dat\n'%outfile)
    f.write('points = %i\n'%points)
    f.write('format = ascii\ninterleaving = field\nfield = locations, field0\n')
    f.write('structure = 2-vector, scalar\ntype = float, float\n\nend\n')
    f.close



#--------------------------------------------------------------------------------
# label.ion.nc - files

def ReadIonNCFile(filename,printnorm=False):
    """
    Reads a NetCDF file that describes the basis orbitals of a given species
    """
    class ion:
        pass
    
    file = NC4.Dataset(filename,'r')
    print 'Reading Basis from %s' % filename

    # General attributes
    ion.filename = filename
    ion.element = file.Element
    ion.label = file.Label
    ion.atomnum = file.Atomic_number
    ion.numorb = file.Number_of_orbitals

    # Variables
    ion.L = N.array(file.variables['orbnl_l'][:],N.int)
    ion.N = N.array(file.variables['orbnl_n'][:],N.int)
    ion.Z = N.array(file.variables['orbnl_z'][:],N.int)
    ion.ispol = N.array(file.variables['orbnl_ispol'][:],N.int)
    ion.orb = N.array(file.variables['orb'][:],N.float)
    ion.cutoff = N.array(file.variables['cutoff'][:],N.float)
    ion.delta = N.array(file.variables['delta'][:],N.float)
    
    print '   Element: %s   Atom number: %i,  L-orbs '% (ion.element,ion.atomnum), ion.L
    for i in range(len(ion.L)):
        rr = ion.delta[i] * N.array(range(len(ion.orb[i])),N.float)
        ion.orb[i] = ion.orb[i]*(rr**ion.L[i])/(PC.Bohr2Ang**(3./2.))
        rr = rr*PC.Bohr2Ang
        # check normalization:
        if printnorm:
            print '   orb %i (L=%i),    Norm = %.6f'%(i,ion.L[i],N.sum(rr*rr*(ion.orb[i]**2))* ion.delta[i]*PC.Bohr2Ang)
    ion.delta = PC.Bohr2Ang*ion.delta
    file.close()
    return ion

def BuildBasis(FDFfile,FirstAtom,LastAtom,lasto):
    """
    Builds the information for each basis orbital in the Hamiltonian
    """
    class basis:
        pass
    CSL = GetFDFblock(FDFfile,'ChemicalSpeciesLabel')
    systemlabel = GetFDFlineWithDefault(FDFfile,'SystemLabel', str, 'siesta', 'io.siesta')
    head,tail =  os.path.split(FDFfile)
    if head=='': 
        head = '.'
    XVfile = '%s/%s.XV'%(head,systemlabel)
    try:
        # XV file prevails
        vectors,speciesnumber,atomnumber,xyz = ReadXVFile(XVfile)
    except:
        vectors,xyz,speciesnumber,atomnumber,natoms = ReadFDFFile(FDFfile)
    ions = {}
    for i in range(len(CSL)):
        # Read ion-nc file for each SIESTA species
        ions[int(CSL[i][0])] = ReadIonNCFile(head+'/%s.ion.nc'%CSL[i][2])

    # Determine the basis dimension nn
    nn = 0
    for i in range(FirstAtom-1,LastAtom): # Python counts from zero
        nn += ions[speciesnumber[i]].numorb

    if nn!=lasto[LastAtom]-lasto[FirstAtom-1]:
        print "Length of basis set build: %i"%nn
        print "Size of Hamiltonian: %i"%(lasto[LastAtom]-lasto[FirstAtom-1])
        print "Error: Could not build basis set. Check if all ion.nc files are there!"
        kuk

    # Initiate basis variables
    basis.ii = N.zeros((nn,),N.int)
    basis.L = N.zeros((nn,),N.int)
    basis.M = N.zeros((nn,),N.int)
    basis.N = N.zeros((nn,),N.int)
    basis.atomnum = N.zeros((nn,),N.int)
    basis.xyz = N.zeros((nn,3),N.float)
    basis.delta = N.zeros((nn,),N.float)
    basis.orb, basis.label = [], []
    basis.coff = N.zeros((nn,),N.float)
    
    # Describe each basis orbital
    iorb = 0
    for ii in range(FirstAtom-1,LastAtom):
        an = atomnumber[ii]
        ion = ions[speciesnumber[ii]]
        for jj in range(len(ion.L)):
            for kk in range(-ion.L[jj],ion.L[jj]+1):
                basis.ii[iorb]=ii+1
                basis.atomnum[iorb]= an
                basis.L[iorb] = ion.L[jj]
                basis.M[iorb] = kk
                basis.N[iorb] = ion.N[jj]
                basis.xyz[iorb,:] = xyz[ii]
                basis.delta[iorb] = ion.delta[jj]
                basis.orb.append(ion.orb[jj,:])
                basis.coff[iorb] = (len(ion.orb[jj,:])-1)*ion.delta[jj]
                basis.label.append(ion.label)
                iorb = iorb+1

    print 'io.siesta.BuildBasis: Generated basis'
    print '... First atom      = %i (Siesta numbering)'%FirstAtom
    print '... Last atom       = %i (Siesta numbering)'%LastAtom
    print '... Basis dimension = %i'%len(basis.L)
    return basis

#--------------------------------------------------------------------------------
# Gaussian Cube files

def ReadCubeFile(filename):
    print 'io.siesta.ReadCubeFile: Reading geometry from',filename
    file = SIO_open(filename,'r')
    # Read comments (lines 1-2)
    comm1 = file.readline()
    comm2 = file.readline()

    # Read number of atoms and coordinate origin (line 3)
    data = string.split(file.readline())
    numberOfAtoms = string.atoi(data[0])
    origin = [string.atof(data[j+1])*PC.Bohr2Ang for j in range(3)]
        
    # Read cell vectors (lines 4-6)
    vox,vectors = [],[]
    for i in range(3):
        data = string.split(file.readline())
        vox.append(string.atoi(data[0]))
        vectors.append([string.atoi(data[0])*string.atof(data[j+1])*PC.Bohr2Ang for j in range(3)])

    # Read geometry (lines 7-7+N)
    atomnumber,xyz = [],[]
    for i in range(numberOfAtoms):
        data = string.split(file.readline())
        atomnumber.append(string.atoi(data[0]))
        xyz.append([string.atof(data[j+2])*PC.Bohr2Ang for j in range(3)])

    file.close()
    return vectors,atomnumber,xyz


#--------------------------------------------------------------------------------
# Consistency checks of SIESTA output files (RUN.out)

import os,time

def GetSiestaStarttime(infile):
    stdin,stdout,stderr = os.popen3('head '+infile)
    for line in stdout:
        if '>>' in line.split():
            starttime = time.strptime(line,'>> Start of run:  %d-%b-%Y  %H:%M:%S     ')
    return starttime

def GetSiestaEndtime(infile):
    stdin,stdout,stderr = os.popen3('tail '+infile)
    for line in stdout:
        if '>>' in line.split():
            endtime = time.strptime(line,'>> End of run:  %d-%b-%Y  %H:%M:%S     ')
    return endtime

def GetSiestaWalltime(infile):
    "Returns the walltime (in hours)"
    try:
        seconds = time.mktime(GetSiestaEndtime(infile))-time.mktime(GetSiestaStarttime(infile))
        hours = seconds/60**2
    except:
        hours = 0.0
        print 'io.siesta.GetSiestaWalltime: Walltime could not be extracted from',infile 
    return hours

def CheckTermination(infile):
    try:
        GetSiestaEndtime(infile)
        return True
    except:
        return False

#--------------------------------------------------------------------------------
# New TSHS file format used for FCrun, TSrun and onlyS

#
#  UNITS! Always eV and Angstrom!
#         k-values always given in range [0,1.0] (or [-0.5,0.5])
#         They are not in reciprocal space. Instead they corresponds
#         to the mathematical orthogonal space that is fourier 
#         transformed.
#

class HS:
    """
    Create full HS from TSHS file. Read fn and assemble for specified k-point
    
    External variables:
    N            : Size of matrices
    H[ispin,i,j] : Hamiltonian
    S[i,j]       : Overlap

    Internal variables from TS (index described as fortran def):
    (NOTE all python lists start with index 0!)
    nua          : Number of atoms in unitcell
    nuo          : Number of orbitals in unitcell
    nspin        : Number of spin
    no           : Number of orbitals in supercell
    maxnh        : Size of sparse matrices
    isa(1:nua)   : Atomic number (not in netcdf)
    lasto(0:nua) : Last orbital of atom in unitcell
    xa(1:nua,1:3): Position of atoms (NOTE transpose of TS standard)
    indxuo(1:no) : Index of equivalent orbital in unitcell
    numh(1:nuo)  : Number of non-zero elements in row of H
    listh(1:nuo) : Column index of H element 
    Hsparse(1:nspin,1:maxnh) : Sparse H
    Ssparse(1:nspin,1:maxnh) : Sparse S
    qtot         : ??
    temp         : ??
    ef           : Fermi energy (why spin dependent?)
    cell(1:3,1:3): Unitcell (ixyz,ivec)
    gamma        : Logical Gamma-point
    xij(1:maxnh,1:3) : Vector between orbital centers (NOTE transpose of TS)
                       NOTE: xij=Rj-Ri where i,j correspond to Hij
    
    Derived internal variables:
    listhptr(1:nuo) : Start of row-1 in sparse matrix 
    atomindx(1:nuo) : Atom index corresponding to orbital in unitcell
    rcell(1:3,1:3) : Reciprocal lattice vectors (ivec,ixyz) (rcell . cell = I)

    For onlyS: Hsparse is not avalable and gamma point is assumed
        gamma: xij is set to 0 and indxuo set manually to 1:nou and -1 for nou+1:no to catch errors!
    """
    def __init__(self,fn,BufferAtoms=N.empty((0,)),UseF90helpers=True):
        self.fn = fn
        if UseF90helpers and fn.endswith('.gz'):
            sys.exit('io.siesta.HS.__init__: F90helpers do not support reading of gzipped TSHS-files. Please unzip and try again.\n')
        
        if UseF90helpers and F90imported:
            print 'io.siesta.HS.__init__: Reading',fn
            self.gamma, self.onlyS, self.nuo, self.no, self.nspin, self.maxnh, self.qtot, \
                self.temp, self.nua, self.ef, self.cell, self.ts_kscell, self.ts_kdispl,\
                self.ts_gamma_scf, self.istep, self.ia1 = F90.readtshs.read(fn)
            # Logical 
            self.gamma, self.onlyS = self.gamma!=0, self.onlyS!=0
            na = BufferAtoms.size
            if na > 0:
                # Remove any buffer-atoms
                ns = self.no / self.nuo
                self.nua, self.nuo, self.maxnh = F90.readtshs.remove_atoms(BufferAtoms)
                self.no = self.nuo * ns
                
            # Arrays
            arr = F90.readtshs
            try:
                self.version = arr.version.copy()
            except:
                self.version = 0
            self.lasto    = arr.lasto.copy()
            self.xa       = arr.xa.copy()
            self.numh     = arr.numh.copy()
            self.listh    = arr.listh.copy()
            self.listhptr = arr.listhptr.copy()
            self.xij      = arr.xij.copy()
            self.Ssparse  = arr.s.copy()
            if not self.onlyS: 
                self.Hsparse = arr.h.copy()
            F90.readtshs.deallocate
        else:
            if BufferAtoms.size > 0:
                raise ValueError("Buffer atoms are not allowed when reading binary files from python")
            general, sparse, matrices = self.__ReadTSHSFile(fn)
            self.nua, self.nuo, self.no , self.nspin, self.maxnh, \
                self.gamma, self.onlyS, self.istep, self.ia1, \
                self.qtot, self.temp, self.ef = general
            self.lasto, self.numh, self.listh, self.indxuo = sparse
        
            if not self.onlyS:
                self.xa, self.cell, self.xij, self.Ssparse, self.Hsparse = matrices
            else:
                self.xa, self.cell, self.xij, self.Ssparse = matrices
        # Adjust memory layout
        self.xa  = N.require(self.xa,requirements=['A','F'])
        self.xij = N.require(self.xij,requirements=['A','F'])
        if not self.onlyS: 
            self.Hsparse = N.require(self.Hsparse,requirements=['A','F'])
        print "Found %i atoms, (%i, %i) orbitals in super-, unit-cell"%(self.nua, self.no, self.nuo)
        self.N = self.nuo
        self.makeDerivedQuant()
        if not self.gamma and not self.onlyS and self.version == 0:
            self.removeUnitCellXij(UseF90helpers) # Remove phase change in unitcell
        self.resetkpoint() # save time by not repeating

    def resetkpoint(self):
        """
        Resets the kpoint and H,S.
        The garbage collector cannot tell if H or S will be used subsequently
        """
        self.kpoint = N.array([1e10,1e10,1e10],N.float)
        try:
            del self.H
        except: pass
        try:
            del self.S
        except: pass

    def __ReadTSHSFile(self,filename):
        """
        Python version for reading TSHS files.
        For return see code:
        Note that onlyS -> does not return xij or Hsparse
                  gamma -> sets indxuo by hand to 1:nou, -1 for nou+1:nos! and sets xij=0.0
        xa[atomnr,xyz] : Atom positions
        ucell[nr,xyz]  : Unitcell
        xij[nr,xyz]
        """
        print 'io.siesta.__ReadTSHSFile: Reading',filename
        self.version = 0
        # Open binary Fortran file
        file = SIO_open(filename,'rb')
        nau,nou,nos,nspin,maxnh = ReadFortranBin(file,fortranLong,5)
        xa = N.reshape(N.array(ReadFortranBin(file,'d',3*nau)),(nau,3))*PC.Bohr2Ang
        xa = N.require(xa.T,requirements=['A','F'])
        isa = ReadFortranBin(file,fortranLong,nau) ; del isa
        ucell = N.transpose(N.reshape(N.array(ReadFortranBin(file,'d',9)),(3,3)))*PC.Bohr2Ang
        gamma = ReadFortranBin(file,'L',1)[0]!=0        # Read boolean (works with ifort)
        onlyS = ReadFortranBin(file,'L',1)[0]!=0        # Read boolean (works with ifort)
        ts_gamma_scf = ReadFortranBin(file,'L',1)[0]!=0 # Read boolean (works with ifort)
        ts_kscell = N.reshape(N.array(ReadFortranBin(file,fortranLong,9)),(3,3))
        ts_kdispl = N.array(ReadFortranBin(file,'d',3))
        istep, ia1 = ReadFortranBin(file,fortranLong,2)
        lasto = N.array(ReadFortranBin(file,fortranLong,nau+1))
        if not gamma:
            indxuo = N.array(ReadFortranBin(file,fortranLong,nos))
        else:
            # For gamma point make indxuo such that indexes not pointing to unitcell give error, i.e., -1. 
            tmp1 = N.array(range(1,nou+1),N.int)
            tmp2 = -N.ones((nos-nou),N.int)
            indxuo = N.concatenate((tmp1,tmp2))
        numhg = N.array(ReadFortranBin(file,fortranLong,nou))
        qtot, temp = ReadFortranBin(file,'d',2)
        temp = temp * PC.Rydberg2eV
        ef = ReadFortranBin(file,'d',1)[0]*PC.Rydberg2eV
        listh = []
        for ii in range(nou):
            listh += list( ReadFortranBin(file,fortranLong,numhg[ii]) )
        listh = N.array(listh)
        Ssparse, cnt = N.zeros(maxnh,N.float), 0
        for ii in range(nou):
            Ssparse[cnt:cnt+numhg[ii]]=ReadFortranBin(file,'d',numhg[ii])
            cnt=cnt+numhg[ii]
        if not onlyS:
            Hsparse = N.zeros((nspin,maxnh),N.float,order='F')
            for ispin in range(nspin):
                cnt=0
                for ii in range(nou):
                    Hsparse[ispin,cnt:cnt+numhg[ii]]=ReadFortranBin(file,'d',numhg[ii])
                    cnt=cnt+numhg[ii]
            Hsparse = Hsparse.T*PC.Rydberg2eV
            Hsparse = N.require(Hsparse,requirements=['A','F'])

        if not gamma:
            # Read xij
            xij = N.zeros((maxnh,3),N.float,order='F')
            cnt=0
            for ii in range(nou):
                tmp=ReadFortranBin(file,'d',numhg[ii]*3)
                tmp = N.reshape(tmp,(3,numhg[ii]))
                xij[cnt:cnt+numhg[ii],:] = tmp.T
                cnt=cnt+numhg[ii]
            xij = xij.T*PC.Bohr2Ang
            xij = N.require(xij,requirements=['A','F'])

        else:
            xij = N.zeros((3,maxnh),N.float,order='F')
        file.close()
        
        general = [nau,nou,nos,nspin,maxnh,gamma,onlyS,istep,ia1,qtot,temp,ef]
        sparse = [lasto, numhg, listh, indxuo]
        if not onlyS:
            matrices = [xa, ucell, xij, Ssparse, Hsparse]
        else:
            matrices = [xa, ucell, xij, Ssparse]
        return general, sparse, matrices


    def makeDerivedQuant(self):
        """
        Create derived internal variables:
        listhptr(1:nuo) : Start of row-1 in sparse matrix 
        atomindx(1:nuo) : Atom index corresponding to orbital in unitcell
        rcell(1:3,1:3)  : Reciprocal lattice vectors (ivec,ixyz) (rcell . cell = I)
        """

        # numh(1:nuo)  : Number of non-zero elements in row of H
        if not 'listhptr' in self.__dict__:
            self.listhptr = N.empty(self.nuo,N.int)
            self.listhptr[0]  = 0
            self.listhptr[1:] = N.cumsum(self.numh[:-1])
        
        # lasto(0:nua) : Last orbital of atom in unitcell       
        self.atomindx = N.empty(self.nuo,N.int)
        atom = 0
        for io in range(self.nuo): 
            while io>=self.lasto[atom]: atom=atom+1 
            self.atomindx[io] = atom

        # Reciprocal cell
        self.rcell = LA.inv(self.cell)

    def removeUnitCellXij(self,UseF90helpers=True):
        """
        Remove displacements within unitcell from xij
        NOTE: We remove the in cell difference so xij corresponds to 
              lattice vectors to the relevant part of the supercell.
        NOTE: xij = Rj-Ri where Ri,j corresponds to positions of the orbitals H_{i,j} 
        TODO: Check why some orbitals in sparse matrix reported within cell but have xij!
        """

        if F90imported and UseF90helpers:
            #      subroutine f90removeunitcellxij( maxnh, no, nuo, nua,
            # +     numh, xij, xa, listhptr, listh, atomindx, xijo)
            self.xij = F90.removeunitcellxij(nnzs=self.maxnh,no_u=self.nuo,
                                                na_u=self.nua,
                                                numh=self.numh,
                                                xij=self.xij,
                                                xa=self.xa,
                                                listh=self.listh,
                                                atomindx=self.atomindx)
        else:
            for iuo in range(self.nuo):
                for jnz in range(self.numh[iuo]):
                    jo=self.listh[self.listhptr[iuo]+jnz]-1
                    juo = self.indxuo[self.listh[self.listhptr[iuo]+jnz]-1]-1
                    ia,ja = self.atomindx[iuo]-1, self.atomindx[juo]-1
                    #if juo==jo and N.max(abs(self.xij[self.listhptr[iuo]+jnz,:]))>0.1:
                    #    print self.xij[self.listhptr[iuo]+jnz,:]
                    self.xij[:,self.listhptr[iuo]+jnz] = self.xij[:,self.listhptr[iuo]+jnz]-\
                        (self.xa[:,ja]-self.xa[:,ia])
                    #if juo==jo and N.max(abs(self.xij[self.listhptr[iuo]+jnz,:]))>0.1:
                    #    print self.xij[self.listhptr[iuo]+jnz,:]
                    #if juo!=jo and N.max(abs(self.xij[self.listhptr[iuo]+jnz,:]))<0.1:
                    #    print self.xij[self.listhptr[iuo]+jnz,:]

    def setkpoint(self,kpoint,UseF90helpers=True,atype=N.complex,verbose=True):
        "Make full matrices from sparse for specific k-point"
        kpoint = N.array(kpoint,N.float)
        if self.gamma:
            VC.Check("same-kpoint", abs(kpoint),
                     "Trying to set non-zero k-point for Gamma point calculation.")
        if N.any(N.abs(self.kpoint-kpoint) > VC.GetCheck("same-kpoint")):
            if verbose: print "io.siesta.HS.setkpoint: %s k =" % self.fn,kpoint
            self.kpoint = kpoint
            self.S = self.setkpointhelper(self.Ssparse,kpoint,UseF90helpers,atype=atype)
            if not self.onlyS:
                self.H = N.empty((self.nspin,self.nuo,self.nuo),atype)
                for ispin in range(self.nspin):
                    self.H[ispin,:,:] = self.setkpointhelper(self.Hsparse[:,ispin],kpoint,UseF90helpers,atype=atype) \
                        - self.ef * self.S


    def setkpointhelper(self,Sparse,kpoint,UseF90helpers=True,atype=N.complex):
        """
        Make full matrices from sparse for specific k-point
        NOTE: Assumption for Fourier transform                     
        Psi(i) =           sum_R exp(i k.R) Psi_k(i)               NOTE sign!
        (Full waveunction)                   (Unit cell part of wavefunction)
              which gives:
        i,j part of unitcell gives from the rows of the full H: 
            H_k(i,j) = H_(i,j) + sum_R H(i,j+R) exp(i*k*R)         NOTE sign!
            where R corresponds to R_j-R_0 which is Xij            NOTE sign!

        For efficiency this routine is normally run in Fortran90. Compile with f2py:
        cd F90;source compile.bat
        """
        if UseF90helpers and F90imported:
            Full = F90.setkpointhelper(nnzs=self.maxnh,sparse=Sparse,kpoint=kpoint, 
                                        no_u=self.nuo,numh=self.numh, 
                                        rcell=self.rcell,xij=self.xij, 
                                        listhptr=self.listhptr,
                                        listh=self.listh)
            # Ensure correct memory alignment
            Full = N.require(Full,requirements=['A','C'])
            Full.shape = (self.nuo,self.nuo)
        else:
            Full = N.zeros((self.nuo,self.nuo),atype)
            # Phase factor 
            tmp = N.dot(kpoint,N.dot(self.rcell,self.xij))
            phase = N.exp(2.0j*N.pi*tmp)    # exp(2 pi i k*(Rj-Ri)) where i,j from Hij

            for iuo in range(self.nuo):
                for jz in range(self.numh[iuo]):
                    si = self.listhptr[iuo]+jz
                    juo = self.indxuo[self.listh[si]-1]-1
                #if juo==self.listh[si]-1:
                #    if phase[si]!=1.0+0.0j:
                #        print "hej"
                    Full[iuo,juo] += Sparse[si]*phase[si]
        if not (Full.dtype==atype):
            print 'io.siesta: Forcing array from %s to %s'%(Full.dtype,atype)
        if atype==N.float or atype==N.float32 or atype==N.float64:
            return N.array(Full.real,atype)
        else:
            return N.array(Full,atype)

# Easy method to read in number of atoms in a TSHS file
def ReadTSHS(fn,**kwargs):
    print 'io.siesta.ReadTSHS: Reading TSHS header from',fn
    # return dictionary
    d = {}
    # Read in TSHS header (do not read in everything!)
    nou,nos,nspin,nua,maxnh = F90.readtshs.read_size(fn)
    if 'nua' in kwargs: d['nua'] = nua
    if 'nuo' in kwargs: d['nuo'] = nuo
    if 'nso' in kwargs: d['nso'] = nos
    if 'nspin' in kwargs: d['nspin'] = nspin
    if 'maxnh' in kwargs: d['maxnh'] = maxnh
    return d

def GetBufferAtomsList(fn,fdf):
    d = ReadTSHS(fn,nua=True)
    nua = d.pop('nua')
    del d
    bufL = GetFDFlineWithDefault(fdf,'TS.BufferAtomsLeft', int, 0,'GetBuffer')
    bufR = GetFDFlineWithDefault(fdf,'TS.BufferAtomsRight', int, 0,'GetBuffer')
    try:
        data = GetFDFblock(fdf, 'TS.Atoms.Buffer')
    except:
        data = []
    BufferAtoms = []
    if len(data) > 0:
        for sl in data:
            sl = [s.lower() for s in sl]
            if sl[0] not in ['position','atom']: continue
            # Currently Inelastica only accepts the from <> to <>
            # and from <> plus/minus <>
            if sl[1] == 'from':
                f = int(sl[2])
                if f < 0: f = nua + f + 1
                t = int(sl[4])
                if sl[3] == 'plus':
                    t = f + t - 1
                elif sl[3] == 'minus':
                    t = f - t + 1
                if t < 0: t = nua + t + 1
                s = 1
                if len(sl) > 5:
                    try:
                        s = int(sl[6])
                    except:
                        s = int(sl[5])
                if f <= t:
                    for i in range(f,t+1,s):
                        BufferAtoms.append(int(i))
                else:
                    for i in range(f,t-1,s):
                        BufferAtoms.append(int(i))
            else:
                for s in sl[1:]:
                    BufferAtoms.append(int(s))
    else:
        for i in range(bufL):
            BufferAtoms.append(i+1)
        for i in range(nua-bufR,nua):
            BufferAtoms.append(i+1)
    if len(BufferAtoms) == 0:
        return N.empty((0,)),0,0
    # sort it
    BufferAtoms.sort()
    nbuf = N.array(BufferAtoms)
    s = 0
    e = nua + 1
    # They must be consecutive
    for i in range(nua+1):
        if i+1 not in nbuf:
            s = i
            break
    e = N.amin(nbuf[nbuf > s])
    for i in range(e,nua+1):
        if i not in nbuf:
            raise ValueError('Buffer atoms must be consecutive, please correct')
    if not N.all(nbuf > 0) and not N.all(nbuf <= nua):
        raise ValueError('Buffer atoms does not exist in the TSHS file')
    return nbuf,s,nua-e+1
