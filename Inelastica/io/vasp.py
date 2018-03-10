"""

:mod:`Inelastica.io.vasp`
=========================

IO interface with VASP

.. currentmodule:: Inelastica.io.vasp

"""

import numpy as N
import string
import gzip


def VIO_open(filename, mode='r'):
    "A io.vasp redefinition of the function open() to handle gzip format"
    try:
        if filename[-3:] == '.gz':
            # filename is explicitly a gzip file
            vfile = gzip.open(filename, mode)
        else:
            # filename is given as a non-zip file
            vfile = open(filename, mode)
    except:
        # if filename is not existing upon read, then try append the '.gz' ending
        vfile = gzip.open(filename+'.gz', mode)
    return vfile

#--------------------------------------------------------------------------------
# Interface with VASP


def ReadCONTCAR(filename):
    "Read CONTCAR file"
    print 'io.vasp.ReadCONTCAR: Reading', filename
    ccarfile = VIO_open(filename, 'r')
    label = ccarfile.readline()
    scalefactor = float(ccarfile.readline())
    vectors = N.zeros((3, 3), N.float)
    for ii in range(3):
        tmp = ccarfile.readline().split()
        vectors[ii] = N.array(tmp, N.float)
    # The species labels are not always included in CONTCAR
    firstline = ccarfile.readline().split()
    line = ccarfile.readline().split()
    try:
        specieslabels = firstline
        speciesnumbers = N.array(line, N.int)
    except:
        speciesnumbers = N.array(firstline, N.int)
        specieslabels = []
        for i, s in enumerate(speciesnumbers):
            specieslabels += [raw_input('Element label for group %i (%i atoms): '%(i+1, s))]
    print 'specieslabels =', specieslabels
    print 'speciesnumbers =', list(speciesnumbers)
    natoms = N.sum(speciesnumbers)
    # Read 'Selective Dynamics' and 'Direct' lines
    dircoor = True # Default is reading direct coordinates
    while line[0].upper() != 'DIRECT' and line[0].upper() != 'CARTESIAN':
        line = ccarfile.readline().split()
        if line[0].upper() == 'CARTESIAN':
            dircoor = False
    # Read coordinates and degrees of freedom
    xyz = N.zeros((natoms, 6), N.float)
    for ii in range(natoms):
        line = ccarfile.readline()
        line = line.replace('F', '0')
        line = line.replace('T', '1')
        line = line.split()
        try:
            xyz[ii] = N.array(line, N.float)
        except:
            # No constraints given
            xyz[ii, :3] = N.array(line, N.float)
    # Ignore rest of the file
    ccarfile.close()
    # Convert to cartesian coordinates
    print 'Read direct coordinates?', dircoor
    for ii in range(natoms):
        if dircoor:
            xyz[ii][:3] = xyz[ii, 0]*vectors[0]+xyz[ii, 1]*vectors[1]+xyz[ii, 2]*vectors[2]
        else:
            xyz[ii][:3] = xyz[ii, :3]
    return label, scalefactor, vectors, specieslabels, speciesnumbers, xyz


def WritePOSCAR(filename, vectors, specieslabels, speciesnumbers, xyz, label='LABEL', scalefactor=1.0, constrained=[]):
    "Write POSCAR file"
    print 'io.vasp.WritePOSCAR: Writing', filename
    pcarfile = open(filename, 'w')
    if label[:-2] != '\n':
        pcarfile.write(label+'\n')
    else:
        pcarfile.write(label)
    pcarfile.write('  %.12f \n'%scalefactor)
    for ii in range(3):
        for jj in range(3):
            pcarfile.write(string.rjust('%.9f'%vectors[ii][jj], 16)+' ')
        pcarfile.write('\n')
    for lbl in specieslabels:
        pcarfile.write('  %s'%lbl)
    pcarfile.write('\n')
    for ii in speciesnumbers:
        pcarfile.write('  %i'%ii)
    pcarfile.write('\n')
    pcarfile.write('Selective dynamics\nCartesian\n')
    for ii, xyzval in enumerate(xyz):
        line = string.rjust('%.9f'%xyzval[0], 16)+' '
        line += string.rjust('%.9f'%xyzval[1], 16)+' '
        line += string.rjust('%.9f'%xyzval[2], 16)+' '
        if len(constrained) > 0:
            for jj in range(3):
                if constrained[ii, jj] > 0:
                    line += ' T'
                else:
                    line += ' F'
            line += '\n'
        else:
            line += ' F F F\n'
        pcarfile.write(line)


def GetEnergies(OUTCAR):
    ocarfile = VIO_open(OUTCAR, 'r')
    print 'io.vasp.GetEnergies: Reading', OUTCAR
    #
    freeE, Etot, EtotSigma0 = 1e100, 1e100, 1e100
    for line in ocarfile:
        if 'free energy    TOTEN' in line:
            l = line.split()
            freeE = float(l[4])      # Pick last appearance
        if 'energy  without entropy=' in line:
            l = line.split()
            Etot = float(l[3])       # Pick last appearance
            EtotSigma0 = float(l[6]) # Pick last appearance
        if 'energy without entropy =' in line:
            l = line.split()
            Etot = float(l[4])       # Pick last appearance
            EtotSigma0 = float(l[7]) # Pick last appearance
    ocarfile.close()
    return freeE, Etot, EtotSigma0


def GetEnergiesFromOszi(OSZICAR):
    oszicarfile = VIO_open(OSZICAR, 'r')
    print 'io.vasp.GetEnergiesFromOszi: Reading', OSZICAR
    #
    f, e0 = 1e100, 1e100
    for line in oszicarfile:
        if 'F=' in line:
            l = line.split()
            f = float(l[2]) # read Free energy
            e0 = float(l[4]) # read energy for sigma->0
    return f, e0


def GetMagnetization(OSZICAR):
    oszicarfile = VIO_open(OSZICAR, 'r')
    print 'io.vasp.GetMagnetization: Reading', OSZICAR
    #
    mag = 1e100
    for line in oszicarfile:
        if 'mag=' in line:
            l = line.split()
            mag = float(l[-1])      # Pick last appearance
    return mag


def GetSpecies(OUTCAR):
    ocarfile = VIO_open(OUTCAR, 'r')
    print 'io.vasp.GetSpecies: Reading', OUTCAR
    atoms = []
    for line in ocarfile:
        if 'TITEL' in line:
            l = line.split()
            print l
            atoms += [l[3]]
    return atoms


def GetVibModesNoScaling(OUTCAR):
    ocarfile = VIO_open(OUTCAR, 'r')
    print 'io.vasp.GetVibrations: Reading', OUTCAR
    freq = []
    modes = []
    v = []
    datablock = False
    for line in ocarfile:
        if 'Eigenvectors and eigenvalues of the dynamical matrix' in line:
            datablock = True # beginning of data block
        if 'Eigenvectors after division by SQRT(mass)' in line \
               or 'Finite differences POTIM=' in line:
            datablock = False # end of data block
        if datablock:
            l = line.split()
            if 'meV' in line:
                # grep frequency
                print line,
                if 'f/i' in line: # imaginary as negative
                    freq.append(-float(l[-2]))
                else:
                    freq.append(float(l[-2]))
                # begin vector array
            if len(l) == 6 and l[0] != 'X':
                v.append([float(l[3]), float(l[4]), float(l[5])])
            if len(l) == 0 and len(v) > 0:
                modes.append(N.array(v))
                v = []
    return N.array(freq), N.array(modes)


def GetVibModesMassScaled(OUTCAR):
    ocarfile = VIO_open(OUTCAR, 'r')
    print 'io.vasp.GetVibrations: Reading', OUTCAR
    freq = []
    modes = []
    v = []
    datablock = False
    for line in ocarfile:
        if 'Finite differences POTIM=' in line:
            datablock = False # end of data block
        if datablock:
            l = line.split()
            if 'meV' in line:
                # grep frequency
                print line,
                if 'f/i' in line: # imaginary as negative
                    freq.append(-float(l[-2]))
                else:
                    freq.append(float(l[-2]))
                # begin vector array
            if len(l) == 6 and l[0] != 'X':
                v.append([float(l[3]), float(l[4]), float(l[5])])
            if len(l) == 0 and len(v) > 0:
                modes.append(N.array(v))
                v = []
        if 'Eigenvectors after division by SQRT(mass)' in line:
            datablock = True
    return N.array(freq), N.array(modes)


def ExtractPDOS(filename, outfile, atom_index=[]):
    "Read DOSCAR file and sum over group of atoms (python numbering)"
    print 'io.vasp.ExtractPDOS: Reading', filename
    f = VIO_open(filename, 'r')
    # Read number of atoms on first line
    s = f.readline()
    s = s.split()
    atoms = int(s[0])
    print 'Atoms =', atoms
    # skip 4 lines
    for i in range(4):
        f.readline()
    # Read Emin,Emax,pts,eF
    s = f.readline()
    s = s.split()
    Emax = float(s[0])
    Emin = float(s[1])
    pts = int(s[2])
    eF = float(s[3])
    print 'Emin,Emax,pts =', Emin, Emax, pts
    print 'eF = ', eF
    # If atom_index not specified take all:
    if atom_index == []:
        atom_index = range(atoms)
    # Loop over atom PDOS
    dat = N.zeros((pts, 19), N.float)
    extrablock = 0
    for j in range(atoms):
        for e in range(pts):
            s = f.readline()
            s = s.split()
            # determine spin deg. freedom
            if e == 0:
                if len(s) == 19:
                    spin = 2
                elif len(s) == 10:
                    spin = 1
                else:
                    # VASP wrote 3-column data...
                    extrablock = 1
            for i, sval in enumerate(s):
                s[i] = float(sval)
            if j-extrablock in atom_index:
                dat[e, 0] = s[0]-eF
                dat[e, 1:1+9*spin] += N.array(s[1:1+9*spin])
                if e == 0:
                    print '  adding %i'%(j-extrablock),
            elif e == 0:
                print '  skipping %i'%(j-extrablock),
        # skip header line
        f.readline()
    if extrablock == 1:
        # Need to read one more block
        j += 1
        for e in range(pts):
            s = f.readline()
            s = s.split()
            # determine spin deg. freedom
            if e == 0:
                if len(s) == 19:
                    spin = 2
                elif len(s) == 10:
                    spin = 1
            for i, sval in enumerate(s):
                s[i] = float(sval)
            if j-extrablock in atom_index:
                dat[e, 0] = s[0]-eF
                dat[e, 1:1+9*spin] += N.array(s[1:1+9*spin])
                if e == 0:
                    print '  adding %i'%(j-extrablock),
            elif e == 0:
                print '  skipping %i'%(j-extrablock),
        # skip header line
        f.readline()
        print
    # Make spin=2 negative
    if spin == 2:
        sgn = N.ones(1+9*spin)
        for i in range(9):
            sgn[2+2*i] = -1.0
        sgn2 = N.array(pts*[sgn])
        dat = dat*sgn2
    # Write output
    print 'io.vasp.ExtractPDOS: Writing', outfile
    fout = open(outfile, 'w')
    for i in range(pts):
        s = ''
        for j in range(1+9*spin):
            s += '%.5e '%dat[i, j]
        fout.write(s+'\n')
    fout.close()
