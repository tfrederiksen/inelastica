import numpy as N
import numpy.linalg as LA
import string, struct, os.path, sys
import MakeGeom as MG
import gzip
import Scientific.IO.NetCDF as nc


#--------------------------------------------------------------------------------
# Interface with VASP

def ReadPOSCAR(filename):
    "Read POSCAR file"
    print 'VaspIO.ReadPOSCAR: Reading', filename
    file = open(filename,'r')
    label = file.readline()
    scalefactor = float(file.readline())
    vectors = N.zeros((3,3),N.float)
    for ii in range(3):
        tmp = file.readline().split()
        vectors[ii] = N.array(tmp,N.float)
    speciesnumbers = N.array(file.readline().split,N.int)
    natoms = N.sum(speciesnumbers)
    xyz = N.zeros((natoms,3),N.float)
    for ii in natoms:
        tmp = file.readline().split()
        xyz[ii] = N.array(tmp,N.float)
    return label,scalefactor,vectors,speciesnumbers,xyz

def WritePOSCAR(filename,vectors,speciesnumbers,xyz,label='LABEL'):
    "Write POSCAR file"
    print 'VaspIO.WritePOSCAR: Writing',filename
    file = open(filename,'w')
    file.write(label+'\n')
    file.write('  %.12f \n'%1.0)
    for ii in range(3):
        for jj in range(3):
            file.write(string.rjust('%.9f'%vectors[ii][jj],16)+' ')
        file.write('\n')
    for ii in range(len(speciesnumbers)):
        file.write('  %i'%speciesnumbers[ii])
    file.write('\n')
    file.write('Selective dynamics\nCartesian\n')
    for ii in range(len(xyz)):
        line=string.rjust('%.9f'%xyz[ii][0],16)+' '
        line+=string.rjust('%.9f'%xyz[ii][1],16)+' '
        line+=string.rjust('%.9f'%xyz[ii][2],16)+' '
        line+='  F  F  F\n'
        file.write(line)


