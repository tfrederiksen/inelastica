"""

:mod:`Inelastica.MakeGeom`
==========================

Routines to read, write and manipulate geometries.

.. currentmodule:: Inelastica.MakeGeom

Classes
-------

.. autosummary::
   :toctree:

   Geom

"""
from __future__ import print_function

import numpy as N
import copy
from . import math
import sys
import netCDF4 as NC4
import ast
import Inelastica.io.siesta as SIO
import Inelastica.io.vasp as VIO
import Inelastica.physics.constants as PC


def interpolateGeom(g0, g1, newlength):
    '''
    Interpolate (extra-) to give new length
        g0 has to have findContacts data
    '''
    try:
        dz = (newlength-g0.ContactSeparation)
    except:
        dz = (newlength-g0.zlength)

    # Find unit cell length
    NN, z0len = 0, 0
    for ii in range(3):
        if g0.pbc[ii][2] > z0len:
            NN, z0len = ii, g0.pbc[ii][2]

    z1len = max(N.array(g1.pbc)[:, 2])

    scale = dz/(z1len-z0len)

    new = Geom()
    for ii in range(g0.natoms):
        r0 = N.array(g0.xyz[ii])
        r1 = N.array(g1.xyz[ii])
        rn = (r1-r0)*scale+r0
        new.addAtom(rn, g0.snr[ii], g0.anr[ii])

    new.pbc = N.array(g0.pbc).copy()
    new.pbc[NN][2] = z0len+dz
    print('MakeGeom.interpolateGeom: Electrode separation interpolated to L = %f' %newlength)
    return new


def CrossProd(A, B):
    "Returns the cross product of two geometric vectors"
    [ax, ay, az] = A
    [bx, by, bz] = B
    return N.array([ay*bz-az*by, az*bx-ax*bz, ax*by-ay*bx])


def GetDist(r1, r2):
    r1 = N.array(r1, N.float)
    r2 = N.array(r2, N.float)
    v12 = r2-r1
    return N.dot(v12, v12)**.5


def GetAngle(r1, r2, r3):
    r1 = N.array(r1, N.float)
    r2 = N.array(r2, N.float)
    r3 = N.array(r3, N.float)
    v12 = r1-r2
    v23 = r3-r2
    angle = math.acos(N.dot(v12, v23)/(N.dot(v12, v12)*N.dot(v23, v23))**.5)
    return 360.0*angle/(2*N.pi)


def GetDihedral(r1, r2, r3, r4):
    r1 = N.array(r1, N.float)
    r2 = N.array(r2, N.float)
    r3 = N.array(r3, N.float)
    r4 = N.array(r4, N.float)
    v12 = r2-r1
    v23 = r3-r2
    v34 = r4-r3
    # atan2-definition of torsion angle with sign
    # see http://en.wikipedia.org/wiki/Dihedral_angle
    y = N.dot(v23, v23)**.5*N.dot(v12, CrossProd(v23, v34))
    x = N.dot(CrossProd(v12, v23), CrossProd(v23, v34))
    angle = math.atan2(y, x)
    return 360.0*angle/(2*N.pi)


class Geom(object):

    '''
    Geom class:
       pbc : repetition vectors
       xyz : list of positions
       snr : list of Siesta number for each atom
       anr : list of atomic numbers
       natoms : number of atoms
    '''

    def __init__(self, fn="", BufferAtoms=N.empty((0,))):
        self.natoms = 0
        self.xyz, self.anr, self.snr = [], [], []

        if len(fn) > 1:
            if fn.endswith('.XV') or fn.endswith('.XV.gz') or fn.endswith('.XV2') or fn.endswith('.XV2.gz'):
                self.readXV(fn)
            elif fn.endswith('.xyz') or fn.endswith('.xyz.gz'):
                self.readXYZ(fn)
            elif fn.endswith('.fdf') or fn.endswith('.fdf.gz'):
                self.readFDF(fn)
            elif fn.endswith('.STRUCT_OUT') or fn.endswith('.STRUCT_OUT.gz'):
                self.readSTRUCT_OUT(fn)
            elif 'CONTCAR' in fn:
                self.readCONTCAR(fn)
            elif 'POSCAR' in fn:
                self.readCONTCAR(fn)
            else:
                sys.exit('ERROR: Input file format unknown\n ... currently .XV, .xyz, .fdf, and POSCAR/CONTCAR are supported')

        try:
            self.constrained
        except:
            self.constrained = 0*N.array(self.xyz)

        # Remove buffer atoms
        if BufferAtoms.size > 0:
            tmp = BufferAtoms.copy() - 1
            self.xyz = N.delete(self.xyz, tmp, axis=0)
            self.anr = N.delete(self.anr, tmp)
            self.snr = N.delete(self.snr, tmp)
            self.constrained = N.delete(self.constrained, tmp, axis=0)

    # BASIC FUNCTIONS

    def rmAtom(self, rmnr):
        "Remove the offending atom"
        if rmnr < 0 or rmnr > self.natoms-1:
            print("ERROR: You tried to remove an atom that isn't there")
            sys.exit(1)
        self.xyz = self.xyz[0:rmnr]+self.xyz[rmnr+1:self.natoms]
        self.snr = self.snr[0:rmnr]+self.snr[rmnr+1:self.natoms]
        self.anr = self.anr[0:rmnr]+self.anr[rmnr+1:self.natoms]
        self.natoms = self.natoms-1

    def addAtom(self, xyz, snr, anr, constrained=[0, 0, 0]):
        "Add atom"
        self.xyz = list(self.xyz)
        self.snr = list(self.snr)
        self.anr = list(self.anr)
        self.xyz.append(xyz[:])
        self.snr.append(snr)
        self.anr.append(anr)
        self.xyz = N.array(self.xyz)
        self.snr = N.array(self.snr)
        self.anr = N.array(self.anr)
        self.constrained = list(self.constrained)+[list(constrained)]
        self.natoms += 1

    def prependAtom(self, xyz, snr, anr, index=0):
        "Add atom"
        self.xyz = list(self.xyz)
        self.snr = list(self.snr)
        self.anr = list(self.anr)
        self.xyz.insert(index, xyz[:])
        self.snr.insert(index, snr)
        self.anr.insert(index, anr)
        self.xyz = N.array(self.xyz)
        self.snr = N.array(self.snr)
        self.anr = N.array(self.anr)
        self.natoms += 1

    def sort(self):
        input('Sure you want to sort the structure?')
        new = Geom()
        tmp = []
        for i in range(self.natoms):
            tmp.append([self.xyz[i][2], self.xyz[i][1], self.xyz[i][0],
                        self.snr[i], self.anr[i], self.constrained[i]])
        tmp.sort()
        for i in range(self.natoms):
            #print i,tmp[i][0],tmp[i][5]
            new.addAtom([tmp[i][2], tmp[i][1], tmp[i][0]], tmp[i][3], tmp[i][4], constrained=tmp[i][5])
        self.xyz = new.xyz
        self.snr = new.snr
        self.anr = new.anr
        self.constrained = N.array(new.constrained)

    def reverse(self):
        # TF/051114
        self.xyz = list(self.xyz)
        self.snr = list(self.snr)
        self.anr = list(self.anr)
        self.xyz.reverse()
        self.snr.reverse()
        self.anr.reverse()
        self.xyz = N.array(self.xyz)
        self.snr = N.array(self.snr)
        self.anr = N.array(self.anr)

    def roundDigits(self, digits=9):
        for i in range(len(self.xyz)):
            for j in range(3):
                self.xyz[i][j] = round(self.xyz[i][j], digits)

    def repeteGeom(self, vec, rep=3):
        # TF 050606: MP 140302 changed repetition order to keep unitcells together
        # TF 140313: Vectorial implementation (much faster for large systems)
        xyz = N.array(self.xyz)
        dr = N.outer(N.ones(len(xyz)), N.array(vec))
        for j in range(1, rep):
            xyz += dr
            for i, xyzval in enumerate(xyz):
                self.addAtom(xyzval, self.snr[i], self.anr[i])

    def addGeom(self, Other):
        #TF/050606
        for i in range(Other.natoms):
            self.addAtom(Other.xyz[i], Other.snr[i], Other.anr[i])

    def prependGeom(self, Other):
        #TF/051114
        for i in range(Other.natoms):
            self.prependAtom(Other.xyz[i], Other.snr[i], Other.anr[i])

    def move(self, dr):
        "Move atoms by dr"
        for ii in range(self.natoms):
            for jj in range(3):
                self.xyz[ii][jj] += dr[jj]

    def move2origo(self):
        self.move([-min(N.array(self.xyz)[:, 0]),
                   -min(N.array(self.xyz)[:, 1]),
                   -min(N.array(self.xyz)[:, 2])])

    def rotate(self, axisvector, angle, RotationCenter=None, RotateSubset=None, Degrees=True,
               RotateLatticeVectors=False):
        # Rotation around an axis specified by some axisvector
        # See the "rotation formula" in Goldstein 2nd ed. p. 165
        vec = N.array(axisvector)/(N.dot(axisvector, axisvector)**0.5) # Normalized unit vector
        if Degrees:
            angle = 2*N.pi/360.0*angle
        if not RotateSubset:
            RotateThese = list(range(self.natoms)) # Rotate all atoms
        else:
            RotateThese = RotateSubset
        if RotationCenter.any():
            RotationCenter = N.array(RotationCenter) # Avoid changes, if reference to an atom
            self.move(-RotationCenter)
        for i in RotateThese:
            r0 = N.array(self.xyz[i])
            r = r0*math.cos(angle) \
                + vec*N.dot(vec, r0)*(1-math.cos(angle)) \
                + CrossProd(r0, vec)*math.sin(angle)
            self.xyz[i] = [r[0], r[1], r[2]]
        if RotationCenter.any():
            # Move back after rotation
            self.move(RotationCenter)
        # Lattice vectors will only be rotated around origo
        if RotateLatticeVectors:
            for i in range(3):
                r0 = N.array(self.pbc[i])
                r = r0*math.cos(angle) \
                    + vec*N.dot(vec, r0)*(1-math.cos(angle)) \
                    + CrossProd(r0, vec)*math.sin(angle)
                self.pbc[i] = [r[0], r[1], r[2]]

    def AlignPlane(self, v1, v2, normal=[0, 0, 1]):
        # Align a plane (specified by v1 and v2) such that "normal" becomes a normal vector
        v1 = N.array(v1, N.float)
        v2 = N.array(v2, N.float)
        normal = N.array(normal, N.float)
        p12 = CrossProd(v1, v2)
        p12n = CrossProd(p12, normal)
        angle = math.acos(N.dot(p12, normal)/(N.dot(p12, p12)*N.dot(normal, normal))**.5)
        self.rotate(p12n, -angle, Degrees=False)

    def PlaceInXYplane(self, atomindices):
        # This function orientates the geometry such that three specified
        # atom indices fall in the same xy-plane
        if len(atomindices) != 3:
            raise ValueError("You need at least 3 atoms here")
        v1 = N.array(self.xyz[atomindices[0]])-N.array(self.xyz[atomindices[1]])
        v2 = N.array(self.xyz[atomindices[0]])-N.array(self.xyz[atomindices[2]])
        self.AlignPlane(v1, v2)

    def GetGeometricCenter(self, Subset=None):
        x0, y0, z0 = 0.0, 0.0, 0.0
        if not Subset:
            sset = self.xyz
        else:
            sset = [self.xyz[i] for i in Subset]
        for coord in set:
            x0 += coord[0]
            y0 += coord[1]
            z0 += coord[2]
        return [x0/len(sset), y0/len(sset), z0/len(sset)]

    def MoveInsideUnitCell(self):
        # Assuming xyz-basis vectors
        #self.move2origo()
        print('MakeGeom.MoveInsideUnitCell: Moving atoms.')
        for i in range(self.natoms):
            self.xyz[i][0] = self.xyz[i][0]%self.pbc[0][0]
            self.xyz[i][1] = self.xyz[i][1]%self.pbc[1][1]
            self.xyz[i][2] = self.xyz[i][2]%self.pbc[2][2]

    def CalcZmatrix(self, first, last):
        'Calculates the Zmatrix for a molecule (SIESTA numbering)'
        zmat = N.zeros((last-first+1, 6), N.float)
        print('MakeGeom.CalcZmatrix: Calculating Zmatrix (from atom %i to %i, SIESTA numbering)...'%(first, last), end=' ')
        f = first-1 # Python numbering
        # 1st - Cartesian coordinates
        if last-first >= 0:
            zmat[0, :3] = N.array([0, 0, 0])
            zmat[0, 3:] = self.xyz[f]

        # 2nd - Spherical coordinates
        if last-first >= 1:
            origo = N.array([0., 0., 0.])
            d = GetDist(self.xyz[f], self.xyz[f+1])
            v01 = N.array(self.xyz[f+1])-N.array(self.xyz[f])
            theta = GetAngle([0., 0., 1.], origo, v01)
            v01[2] = 0.0 # project into xy-plane
            phi = GetAngle([1., 0., 0.], origo, v01)
            zmat[1, :3] = N.array([1, 0, 0])
            zmat[1, 3:] = N.array([d, theta, phi])

        # 3rd - Dihedral angle with pseudoatom coordinate r0
        if last-first >= 2:
            r0 = N.array(self.xyz[f])+N.array([0., 0., 10.])
            zmat[2, :3] = N.array([2, 1, 0])
            zmat[2, 3:] = N.array([GetDist(self.xyz[f+1], self.xyz[f+2]),
                                   GetAngle(self.xyz[f], self.xyz[f+1], self.xyz[f+2]),
                                   GetDihedral(r0, self.xyz[f], self.xyz[f+1], self.xyz[f+2])])
        # Remaining atoms
        if last-first >= 3:
            for i in range(f+3, last):
                zmat[i-f, :3] = N.array([i-f, i-f-1, i-f-2])
                zmat[i-f, 3:] = N.array([GetDist(self.xyz[i-1], self.xyz[i]),
                                         GetAngle(self.xyz[i-2], self.xyz[i-1], self.xyz[i]),
                                         GetDihedral(self.xyz[i-3], self.xyz[i-2], self.xyz[i-1], self.xyz[i])])
        print('Done!')
        return zmat

    # PURPOSE SPECIFIC FUNCTIONS

    def findContactsAndDevice(self, AtomsPerLayer, tol=1e-4):
        '''
        If there exists AtomsPerLayer of atoms with the same z-coordinate
        these are assumed to be belonging to the contacts.
        The remaining atoms are taken to be the device.
        The left contact atoms are assumed to appear before the device,
        and the right contact atoms after the device.

        The function identifies and sets the following variables:
        - leftContactList
        - rightContactList
        - deviceList
        - zLeftContact
        - ContactSeparation
        '''

        self.leftContactList = []
        self.rightContactList = []
        self.deviceList = []
        if AtomsPerLayer > 1:
            for i in range(self.natoms):
                tmp = []
                for j in range(self.natoms):
                    if abs(self.xyz[i][2]-self.xyz[j][2]) < tol: tmp.append(j)
                if len(tmp) == AtomsPerLayer and len(self.deviceList) == 0:
                    self.leftContactList.append(i+1)  # Siesta numbering starts from 1
                elif len(tmp) == AtomsPerLayer and len(self.deviceList) > 0:
                    self.rightContactList.append(i+1) # Siesta numbering starts from 1
                else:
                    self.deviceList.append(i+1) # Siesta numbering starts from 1
        else: # AtomsPerLayer == 1
            lsep = self.xyz[1][2]-self.xyz[0][2]
            # First atom belongs per definition to the left contact:
            self.leftContactList.append(1) # Siesta numbering starts from 1
            for i in range(1, self.natoms-1):
                distL = self.xyz[i][2]-self.xyz[i-1][2]
                distR = self.xyz[i+1][2]-self.xyz[i][2]
                if abs(distL-lsep) < tol and len(self.rightContactList) == 0:
                    self.leftContactList.append(i+1)
                elif abs(distR-lsep) < tol:
                    self.rightContactList.append(i+1)
                else:
                    self.deviceList.append(i+1)
            # Last atom belongs per definition to the right contact:
            self.rightContactList.append(self.natoms) # Siesta numbering starts from 1
        if len(self.leftContactList) > 0 and len(self.rightContactList) > 0:
            self.zLeftContact = self.xyz[self.leftContactList[-1]-1][2]
            self.ContactSeparation = self.xyz[self.rightContactList[0]-1][2] \
                                     - self.zLeftContact
        elif len(self.leftContactList) > 0:
            self.zLeftContact = self.xyz[self.leftContactList[-1]-1][2]
            self.ContactSeparation = max(N.array(self.pbc)[:, 2]) \
                                     - self.xyz[self.leftContactList[-1]-1][2]
        elif len(self.rightContactList) > 0:
            self.zLeftContact = 0.0
            self.ContactSeparation = self.xyz[self.rightContactList[0]-1][2]
        elif len(self.leftContactList) == 0 and len(self.rightContactList) == 0:
            self.zLeftContact = 0.0
            self.ContactSeparation = max([self.pbc[ii][2] for ii in range(3)])
        print('MakeGeom.findContactsAndDevice: Electrode separation detected was L = %f' \
              %self.ContactSeparation)
        counts = len(self.leftContactList), len(self.deviceList), \
                 len(self.rightContactList), len(self.xyz)
        if counts[0]+counts[1]+counts[2] != counts[3]:
            raise ValueError("Non conforming")
        print('   ... Atoms (Left/Device/Right):  %i / %i / %i  =  %i'%counts)
        print('   ... DeviceFirst, DeviceLast = %i, %i  (Siesta numbering)'\
              %(self.deviceList[0], self.deviceList[-1]))

    def stretch2NewContactSeparation(self, NewContactSeparation, AtomsPerLayer,
                                     ListL=None, ListR=None):
        "Stretch system to NewContactSeparation"
        self.findContactsAndDevice(AtomsPerLayer)
        print('MakeGeom.stretch2NewContactSeparation:')
        print('   ... Stretching from L = %f Ang --> L = %f Ang' \
              %(self.ContactSeparation, NewContactSeparation))
        if ListL: print('   ... ListL = [%i:%i] (Siesta numbering)'%(ListL[0], ListL[-1]))
        else: ListL = self.leftContactList
        if ListR: print('   ... ListR = [%i:%i] (Siesta numbering)'%(ListR[0], ListR[-1]))
        else: ListR = self.rightContactList

        # Scale device coordinates
        for ii in self.deviceList:
            if (ii not in ListL) and (ii not in ListR):
                self.xyz[ii-1][2] = (self.xyz[ii-1][2]-self.zLeftContact)*\
                                    (NewContactSeparation/self.ContactSeparation)+self.zLeftContact
            if ii in ListR:
                self.xyz[ii-1][2] = self.xyz[ii-1][2]+NewContactSeparation-self.ContactSeparation
        # Move right contact to new position
        for ii in self.rightContactList:
            self.xyz[ii-1][2] = self.xyz[ii-1][2]+NewContactSeparation-self.ContactSeparation
        # Set new size of unit cell
        NN, zmax = 0, 0
        for ii in range(3):
            # Find cell vector with largest z-component
            if self.pbc[ii][2] > zmax:
                NN, zmax = ii, self.pbc[ii][2]
        self.pbc[NN][2] = zmax+NewContactSeparation-self.ContactSeparation
        self.ContactSeparation = NewContactSeparation

    def PasteElectrodeLayers(self, BlockSize, AtomsPerLayer, LayersLeft, LayersRight):
        """
        Copy atoms from range(BlockSize) as many times as specified by
        BlocksLeft/BlocksRight, and including enlargement of the unit cell.
        """
        # Determine interlayer separation
        LayersPerBlock = BlockSize/AtomsPerLayer
        LayerSep = self.xyz[AtomsPerLayer][2]-self.xyz[0][2]

        print('MakeGeom.PasteElectrodeLayers:')
        print('   ... Initial number of atoms      = %i' %len(self.xyz))
        print('   ... AtomsPerLayer                = %i atoms' %AtomsPerLayer)
        print('   ... Block size                   = %i atoms' %BlockSize)
        print('   ... Block periodicity            = %.4f Ang'%(LayersPerBlock*LayerSep))
        print('   ... Layer separation             = %.4f Ang'%LayerSep)
        print('   ... Layer(s) pasted to left      = %i' %LayersLeft)
        print('   ... Layer(s) pasted to right     = %i' %LayersRight)

        # Make block geometry
        pieceR = Geom()
        pieceL = Geom()
        for i in range(BlockSize):
            pieceR.addAtom(self.xyz[i], self.snr[i], self.anr[i])
            pieceL.addAtom(self.xyz[i], self.snr[i], self.anr[i])

        # Append layers to the right
        BlocksRight = LayersRight/LayersPerBlock
        ExtraLayersRight = LayersRight%LayersPerBlock
        pieceR.move(self.pbc[2])
        for i in range(BlocksRight):
            self.addGeom(pieceR)
            pieceR.move([0, 0, LayersPerBlock*LayerSep])
        for i in range(AtomsPerLayer*ExtraLayersRight):
            self.addAtom(pieceR.xyz[i], pieceR.snr[i], pieceR.anr[i])
        self.pbc[2][2] += (BlocksRight*LayersPerBlock+ExtraLayersRight)*LayerSep

        # Prepend layers to the left
        BlocksLeft = LayersLeft/LayersPerBlock
        ExtraLayersLeft = LayersLeft%LayersPerBlock
        pieceL.reverse()
        for i in range(BlocksLeft):
            pieceL.move([0, 0, -LayersPerBlock*LayerSep])
            self.prependGeom(pieceL)
        pieceL.move([0, 0, -LayersPerBlock*LayerSep])
        for i in range(AtomsPerLayer*ExtraLayersLeft):
            self.prependAtom(pieceL.xyz[i], pieceL.snr[i], pieceL.anr[i])
        self.pbc[2][2] += (BlocksLeft*LayersPerBlock+ExtraLayersLeft)*LayerSep

        print('   ... Final number of atoms        = %i' %len(self.xyz))

    def BuildOnlyS(self, displacement=0.02):
        "Returns a new geometry object with 7 times as many atoms"
        geom = copy.deepcopy(self)
        for dim in range(3):
            for displdir in [-1, 1]:
                new = copy.deepcopy(self)
                displ = [0., 0., 0.]
                displ[dim] = displdir*displacement
                new.move(displ)
                geom.addGeom(new)
        return geom

    def ShiftPeriodicCellAlongZ(self, IndexShift):
        print('MakeGeom.ShiftPeriodicCellAlongZ: Shifting %i atoms from LEFT to RIGHT side of the unit cell'%IndexShift)
        # Find cell vector with largest z-component
        zmax = 0
        for ii in range(3):
            if self.pbc[ii][2] > zmax:
                zmax = self.pbc[ii][2]
        for i in range(IndexShift):
            # add z-periodicity to z-coordinate
            self.xyz[i][2] += zmax
        self.sort()
        # Bring first layer to z=0.0
        self.move([0, 0, -min(N.array(self.xyz)[:, 2])])

    def PasteLists(self, AddLeftList=[], AddRightList=[]):
        # Extend a geometry with coordinates in the specified lists
        # along the z-direction, automatically shifting the z-coordinates
        # to match the lattice structure. The supercell vector
        # in the z-direction is enlarged correspondingly.
        AddLeftList = N.array(AddLeftList)
        AddRightList = N.array(AddRightList)
        if len(AddLeftList) > 0:
            dz = AddLeftList[AtomsPerLayer, 2]-AddLeftList[0, 2]
            tmp = N.array(geom.xyz)
            minz = min(tmp[:, 2])
            maxz = max(AddLeftList[:, 2])
            for ii in reversed(list(range(len(AddLeftList)))):
                geom.prependAtom(AddLeftList[ii, :]+
                                 (-maxz+minz-dz)*N.array([0, 0, 1], N.float),
                                 geom.snr[0], geom.anr[0])
            geom.pbc[2][2] += len(AddLeftList)/AtomsPerLayer*dz
        if len(AddRightList) > 0:
            tmp = N.array(geom.xyz)
            maxz = max(tmp[:, 2])
            minz = min(AddRightList[:, 2])
            for ii in range(len(AddRightList)):
                geom.addAtom(AddRightList[ii, :]+
                             (maxz-minz+dz)*N.array([0, 0, 1], N.float),
                             geom.snr[0], geom.anr[0])
            geom.pbc[2][2] += len(AddRightList)/AtomsPerLayer*dz

    def StretchAlongEigenvector(self, ncfilename, modeindex, displacement=0.04*PC.Bohr2Ang):
        'Displace a geometry along a calculated eigenmode'
        # TF/080527
        ncfile = NC4.Dataset(ncfilename, 'r')
        hw = ncfile.variables['hw'][modeindex]
        U = ncfile.variables['U'][modeindex]
        print('MakeGeom.StretchAlongEigenvector: Stretching %.3e Ang along mode #%i (hw = %.4f eV).'%(displacement, modeindex, hw))
        d = len(U) // 3
        U = N.reshape(U, (d, 3))
        firstdynatom = int(ncfile.variables['DynamicAtoms'][0]-1)
        for i in range(d):
            self.xyz[i+firstdynatom] += displacement*N.array(U[i])

    # GENERAL IO-FUNCTIONS

    def readXV(self, fn):
        "Read XV file"
        self.pbc, self.snr, self.anr, self.xyz = SIO.ReadXVFile(fn)
        self.natoms = len(self.xyz)
        #self.move2origo()

    def writeXV(self, fn, rep=[1, 1, 1]):
        geom = copy.deepcopy(self)
        for i in range(3):
            geom.repeteGeom(self.pbc[i], rep=rep[i])
            geom.pbc[i] = [rep[i]*x for x in self.pbc[i]]
        SIO.WriteXVFile(fn, geom.pbc, geom.snr, geom.anr, geom.xyz)

    def readXYZ(self, fn):
        label, self.anr, self.xyz = SIO.ReadXYZFile(fn)
        self.natoms = len(self.xyz)
        #self.move2origo()
        self.snr = ast.literal_eval(input('Input snr expression:'))
        if len(self.xyz) != len(self.snr):
            print('Error assigning snr!')
        self.pbc = []
        for i in range(3):
            vec = ast.literal_eval(input('Input cell vector (%i):'%i))
            if len(vec) != 3:
                print('Error assigning cell vector!')
            else:
                self.pbc.append(vec)

    def writeXYZ(self, fn, rep=[1, 1, 1], write_ghosts=False):
        "Write .mkl file with unitcell repeated rep times"
        geom = copy.deepcopy(self)
        for i in range(3):
            geom.repeteGeom(self.pbc[i], rep=rep[i])
            geom.pbc[i] = [rep[i]*x for x in self.pbc[i]]
        SIO.WriteXYZFile(fn, geom.anr, geom.xyz, write_ghosts)

    def readFDF(self, fn):
        self.pbc, self.xyz, self.snr, self.anr, self.natoms = SIO.ReadFDFFile(fn)

    def writeFDF(self, fn, rep=[1, 1, 1]):
        "Write STRUCT.fdf file"
        geom = copy.deepcopy(self)
        for i in range(3):
            geom.repeteGeom(self.pbc[i], rep=rep[i])
            geom.pbc[i] = [rep[i]*x for x in self.pbc[i]]
        SIO.WriteFDFFile(fn, geom.pbc, geom.snr, geom.anr, geom.xyz)

    def writeFDFZmat(self, fn, first=0, last=0, rep=[1, 1, 1]):
        "Write STRUCT.fdf file"
        geom = copy.deepcopy(self)
        for i in range(3):
            geom.repeteGeom(self.pbc[i], rep=rep[i])
            geom.pbc[i] = [rep[i]*x for x in self.pbc[i]]
        if first > 0 and last >= first:
            zmat = geom.CalcZmatrix(first, last)
        else:
            zmat = None
        SIO.WriteFDFFileZmat(fn, geom.pbc, geom.snr, geom.anr, geom.xyz, first, last, zmat)

    def readMKL(self, fn):
        "not implemented yet"

    def writeMKL(self, fn, rep=[1, 1, 1]):
        "Write .mkl file with unitcell repeated rep times"
        geom = copy.deepcopy(self)
        for i in range(3):
            geom.repeteGeom(self.pbc[i], rep=rep[i])
            geom.pbc[i] = [rep[i]*x for x in self.pbc[i]]
        SIO.WriteMKLFile(fn, geom.anr, geom.xyz, [], [], 0, 0)

    def readSTRUCT_OUT(self, fn):
        self.pbc, self.snr, self.anr, self.xyz = SIO.ReadSTRUCT_OUTFile(fn)

    def readCONTCAR(self, fn):
        "Read geometry from VASP CONTCAR file"
        label, scalefactor, vectors, specieslabels, speciesnumbers, xyz = VIO.ReadCONTCAR(fn)
        self.pbc = N.array(vectors)*scalefactor
        self.xyz = N.array(xyz[:, :3])*scalefactor
        self.constrained = N.array(xyz[:, 3:])
        self.natoms = len(xyz)
        self.snr = []
        self.anr = []
        for i, nspecie in enumerate(speciesnumbers):
            self.snr += nspecie*[i+1]
            self.anr += nspecie*[PC.PeriodicTable[specieslabels[i]]]

#--------------------------------------------------------------------------------
# Interface with VASP

    def writePOSCAR(self, fn, constrained=[]):
        "Write POSCAR coordinate file for VASP"
        geom = copy.deepcopy(self)
        tmp = []
        anrnum = N.zeros(max(geom.anr)+1, N.int)
        for i, xyz in enumerate(geom.xyz):
            if geom.anr[i] > 0:
                tmp += [[geom.anr[i], i]]
                anrnum[geom.anr[i]] += 1 # count atoms of each element type
        tmp.sort() # sort according to increasing atom number
        xyz = []
        cons = []
        for t in tmp:
            j = t[1] # sorted index
            xyz += [geom.xyz[j]]
            cons += [list(geom.constrained[j])]
        speciesnumbers = []
        specieslabels = []
        for i, nanr in enumerate(anrnum):
            if nanr != 0:
                speciesnumbers += [nanr]
                specieslabels += [PC.PeriodicTable[i]]
        VIO.WritePOSCAR(fn, geom.pbc, specieslabels, speciesnumbers, xyz, constrained=N.array(cons))

#--------------------------------------------------------------------------------
# Conversion functions:
#


def convert(infile, outfile, rep=[1, 1, 1], MoveInsideUnitCell=False, RoundOff=False):
    """Reads and writes according to the file extensions specified."""
    geom = Geom(infile)
    if MoveInsideUnitCell: geom.MoveInsideUnitCell()
    if RoundOff:
        print('MakeGeom.convert: Performing rounding...')
        for i in range(geom.natoms):
            for j in range(3):
                geom.xyz[i][j] += 1e5
                geom.xyz[i][j] = round(geom.xyz[i][j], 4)
                geom.xyz[i][j] -= 1e5
    if outfile[-4:] == '.XYZ':
        geom.writeXYZ(outfile, rep, write_ghosts=True)
    elif outfile[-4:] == '.xyz':
        geom.writeXYZ(outfile, rep, write_ghosts=False)
    elif outfile[-4:] == '.fdf':
        geom.writeFDF(outfile, rep)
    elif outfile[-3:] == '.XV':
        geom.writeXV(outfile, rep)
    elif outfile[-4:] == '.mkl':
        geom.writeMKL(outfile, rep)
    elif outfile[-6:] == 'POSCAR':
        print('Repetition not implemented with the POSCAR conversion')
        geom.writePOSCAR(outfile, constrained=geom.constrained)
    else:
        print('Error - did not recognize output format in %s' %outfile)


def repeteANI(ANIinfile, XVinfile, outfile, rep=[1, 1, 1]):
    """Generates a new ANI-file with cell repetitions"""
    GeomList, Energy = SIO.ReadANIFile(ANIinfile)
    initGeom = Geom(XVinfile)
    newGeomList = []
    for i in range(len(GeomList)):
        geom = copy.deepcopy(GeomList[i])
        geom.pbc = N.zeros((3, 3), N.float)
        for i in range(3):
            geom.repeteGeom(initGeom.pbc[i], rep=rep[i])
            geom.pbc[i] = [rep[i]*x for x in initGeom.pbc[i]]
        newGeomList.append(geom)
    SIO.WriteANIFile(outfile, newGeomList, Energy)
