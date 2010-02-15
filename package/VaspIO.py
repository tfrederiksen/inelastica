import numpy as N
import numpy.linalg as LA
import string, struct, os.path, sys
import MakeGeom as MG
import gzip
import Scientific.IO.NetCDF as nc


#--------------------------------------------------------------------------------
# Interface with VASP

def ReadCONTCAR(filename):
    "Read CONTCAR file"
    print 'VaspIO.ReadCONTCAR: Reading', filename
    file = open(filename,'r')
    label = file.readline()
    scalefactor = float(file.readline())
    vectors = N.zeros((3,3),N.float)
    for ii in range(3):
        tmp = file.readline().split()
        vectors[ii] = N.array(tmp,N.float)
    speciesnumbers = N.array(file.readline().split(),N.int)
    natoms = N.sum(speciesnumbers)
    # Read 'Selective Dynamics' and 'Direct' lines
    file.readline()
    file.readline()
    # Read coordinates and degrees of freedom
    xyz = N.zeros((natoms,6),N.float)
    for ii in range(natoms):
        line = file.readline()
        line = line.replace('F','0')
        line = line.replace('T','1')
        line = line.split()
        xyz[ii] = N.array(line,N.float)
    # Ignore rest of the file
    file.close()
    # Convert to cartesian coordinates
    for ii in range(natoms):
        xyz[ii][:3] = xyz[ii,0]*vectors[0]+xyz[ii,1]*vectors[1]+xyz[ii,2]*vectors[2]
    return label,scalefactor,vectors,speciesnumbers,xyz

def WritePOSCAR(filename,vectors,speciesnumbers,xyz,label='LABEL',scalefactor=1.0):
    "Write POSCAR file"
    print 'VaspIO.WritePOSCAR: Writing',filename
    file = open(filename,'w')
    file.write(label)
    file.write('  %.12f \n'%scalefactor)
    for ii in range(3):
        for jj in range(3):
            file.write(string.rjust('%.9f'%vectors[ii][jj],16)+' ')
        file.write('\n')
    for ii in range(len(speciesnumbers)):
        file.write('  %i'%speciesnumbers[ii])
    file.write('\n')
    file.write('Selective dynamics\nCartesian\n')
    for ii in range(len(xyz)):
        line  = string.rjust('%.9f'%xyz[ii][0],16)+' '
        line += string.rjust('%.9f'%xyz[ii][1],16)+' '
        line += string.rjust('%.9f'%xyz[ii][2],16)+' '
        for jj in range(3):
            if xyz[ii,3+jj] == 1.0:
                line += '  T'
            else:
                line += '  F'
        file.write(line+'\n')


def GetEnergies(OUTCAR):
    file = open(OUTCAR,'r')
    #
    freeE, Etot, EtotSigma0 = 1e100, 1e100, 1e100
    for line in file:
        if 'TOTEN' in line:
            l = line.split()
            freeE = float(l[4])      # Pick last appearance
        if 'energy  without entropy=' in line:
            l = line.split()
            Etot = float(l[3])       # Pick last appearance
            EtotSigma0 = float(l[6]) # Pick last appearance
    file.close()
    return freeE, Etot, EtotSigma0
