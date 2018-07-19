"""

:mod:`Inelastica.NEB`
=====================

Nudged elastic band

.. currentmodule:: Inelastica.NEB

Classes
-------

.. autosummary::
   :toctree:

   general
   savedData
   step

=======================
https://doi.org/10.1142/9789812839664_0016
"""
from __future__ import print_function

import numpy as N
import numpy.linalg as LA
import sys, glob, os, copy, time, pickle, string
import Inelastica.io.siesta as SIO
import Inelastica.SetupRuns as SUR
import Inelastica.MakeGeom as MG
import Inelastica.math as MM
import Inelastica.physics.constants as PC
import Inelastica.misc.valuecheck as VC


##################### Global variabels #################
# steps, opts
class savedData(object):
    pass

########################################################
##################### Main routine #####################
########################################################


def main(options):
    global opts
    opts = options
    runNEB()


########################################################
def runNEB():
    global opts, steps, savedData

    # Restarting?
    fns = glob.glob('NEB_*/CGrun')
    if len(fns) == 0:
        restart = False
    elif len(fns) != opts.NNEB:
        sys.exit('Number of NEB_??? directories %i does not match number of steps %i'%(len(fns), opts.NNEB))
    else:
        if opts.onlyInit:
            sys.exit('NEB_* directories already exist, exititng.')
        else:
            restart = True

    fns = ["NEB_%i/CGrun"%ii for ii in range(opts.NNEB)]
    fns = [opts.initial+"/CGrun"]+fns+[opts.final+"/CGrun"]

    i, f = step(fns[0], False, 0), step(fns[opts.NNEB+1], False, opts.NNEB+1)
    if restart:
        try:
            savedData.E, savedData.F, savedData.Fmax, \
                savedData.FmaxWOspring, savedData.geom = pickle.load(open('NEB_0/savedData.pickle', 'r'))
        except:
            savedData.E = []
            savedData.F = []
            savedData.geom = []
            savedData.Fmax = []
    else:
        savedData.E = []
        savedData.F = []
        savedData.geom = []
        savedData.Fmax = []

    steps = [step(fn, restart, ii, initial=i, final=f) for ii, fn in enumerate(fns)]

    checkConst(i, f)
    if restart:
        for ii in range(len(steps)-1):
            checkConst(steps[ii],steps[ii+1])

    done = opts.onlyInit

    while not done:
        # Start calculations
        for ii in steps:
            ii.run()

        nextstep = False
        while not nextstep:
            for ii in steps:
                if not ii.checkDone(): # Reads forces!
                    break
            else:
                nextstep = True
            if not nextstep:
                time.sleep(10)

        savedData.E += [[SIO.GetTotalEnergy(ii.dirr+'/RUN.out') for ii in steps]]
        savedData.geom += [[copy.deepcopy(ii.XVgeom) for ii in steps]]

        overShoot = False
        for ii, jj in enumerate(steps[:-1]):
            if ii != 0: # Skip first and last
                overShoot = overShoot or jj.shootOver(savedData.geom[-1][ii-1], savedData.geom[-1][ii+1])
        if overShoot: # Kill speed if one image overshoot ... otherwise kinks will appear
            for ii in steps:
                ii.v = ii.v*0

        for ii, jj in enumerate(steps[:-1]):
            if ii != 0: # Skip first and last
                jj.update(savedData.geom[-1][ii-1], savedData.geom[-1][ii+1])

        # Write status
        savedData.F += [[ii.F for ii in steps]]
        savedData.Fmax += [[ii.Fmax for ii in steps]]
        for ii in range(1, len(steps)-1):
            geoms = [savedData.geom[jj][ii] for jj in range(len(savedData.geom))]
            Fmax = [savedData.Fmax[jj][ii] for jj in range(len(savedData.geom))]
            Fs = [savedData.F[jj][ii] for jj in range(len(savedData.geom))]
            SIO.WriteANIFile('NEB_%i/Steps.ANI'%(ii-1), geoms, Fmax)
            SIO.WriteAXSFFiles('NEB_%i/Steps.XASF'%(ii-1), geoms, forces=Fs)

        geoms = [ii.XVgeom for ii in steps]
        Fmax = [ii.Fmax for ii in steps]
        F = [ii.F for ii in steps]
        SIO.WriteANIFile('NextStep.ANI', geoms, \
                         [ii-savedData.E[-1][0] for ii in savedData.E[-1]])
        SIO.WriteAXSFFiles('NextStep.XASF', geoms, forces=F)
        done = True
        for ii in steps:
            done = done and ii.converged
        pickle.dump((savedData.E, savedData.F, savedData.Fmax, savedData.geom),\
                    open('NEB_%i/savedData.pickle'%0, 'w'))
        f = open('Convergence','a')
        f.write(('####### Iteration %i #######\n#Fmax '+('%2.3f '*(opts.NNEB+2))+'\n')%\
                tuple([len(savedData.Fmax),]+savedData.Fmax[-1]))
        f.write(('#step length '+(('%2.4f ')*(opts.NNEB+2))+'\n')%\
                tuple([LA.norm(ii.v) for ii in steps]))
        f.write('# Barrier [meV]:\n')
        f.write((('%4.1f '*(opts.NNEB+2))+'\n')%tuple([1000*(savedData.E[-1][ii]-savedData.E[-1][0])\
                                                       for ii in range(opts.NNEB+2)]))     
        f.write('# delta E compared to start/restart [meV]:\n')
        f.write((('%4.1f '*(opts.NNEB+2))+'\n')%tuple([1000*(savedData.E[0][ii]-savedData.E[-1][ii])\
                                                       for ii in range(opts.NNEB+2)]))
        f.close()

#################### Class for each step ###############


class step(object):
    global steps, opts

    def __init__(self, dirr, restart, iistep, initial=None, final=None):
        self.dirr = dirr
        self.converged, self.Fmax = False, 0.0
        self.ii = iistep
        self.fixed = (iistep == 0) or (iistep == opts.NNEB+1)

        if restart or self.fixed:
            self.FDFgeom = MG.Geom(self.dirr+"/"+opts.fn)
            self.XVgeom = readxv(self.dirr)
            self.v = N.array(self.FDFgeom.xyz)*0

        if not restart and not self.fixed:
            os.makedirs(dirr)
            SUR.CopyInputFiles(opts.initial+"/CGrun/", dirr,\
                                   ['.fdf', '.vps', '.psf'])
            # Interpolate
            ixyz, fxyz = N.array(initial.XVgeom.xyz), N.array(final.XVgeom.xyz)
            mix = float(iistep)/(opts.NNEB+1.0)
            xyz = (1-mix)*ixyz+mix*fxyz
            self.FDFgeom = copy.copy(initial.XVgeom)
            self.FDFgeom.xyz = [xyz[ii, :] for ii in range(len(xyz))]
            self.FDFgeom.writeFDF(self.dirr+"/STRUCT.fdf")
            self.v = N.array(self.FDFgeom.xyz)*0

            # Append lines to RUN.fdf
            elm = dirr+"/"+opts.fn
            f = open(elm, 'r')
            lines = f.readlines()
            f.close()
            f = open(elm, 'w')
            f.write('### Lines appended %s \n' %time.ctime())
            f.write("MD.TypeOfRun         CG\n")
            f.write("MD.NumCGsteps        0\n")
            f.write("MD.UseSaveXV         False\n")
            f.write("MD.UseSaveCG         False\n")
            f.write('# end of lines appended %s \n' %time.ctime())
            for line in lines:
                f.write(line)
            f.close()

        self.done = self.checkDone()
        const = SIO.GetFDFblock(dirr+"/"+opts.fn, "GeometryConstraints")
        if opts.const2 != None:
            self.const = [[opts.const2[0],opts.const2[0]]]
        else:
            self.const = []
        for ii in const:
            self.const += [[int(ii[2])-1, int(ii[4])-1]]

    def checkDone(self):
        if self.fixed == True:
            self.FDFgeom = MG.Geom(self.dirr+"/"+opts.fn)
            self.XVgeom = readxv(self.dirr)
            self.forces = SIO.ReadForces(self.dirr+"/RUN.out")
            self.forces = self.forces[-len(self.XVgeom.xyz):]
            self.F = N.array(self.forces)*0
            self.v = N.array(self.forces)*0
            self.converged = True
            return True
        else:
            done = SIO.CheckTermination(self.dirr+"/RUN.out")
            if done:
                try:
                    self.FDFgeom = MG.Geom(self.dirr+"/"+opts.fn)
                    self.XVgeom = readxv(self.dirr)
                    self.forces = SIO.ReadForces(self.dirr+"/RUN.out")
                    if N.allclose(self.XVgeom.xyz, self.FDFgeom.xyz, 1e-6):
                        done = True
                except:
                    done = False
            return done

    def run(self):
        if (not self.done) and (not self.converged):
            try:
                os.remove(self.dirr+"/RUN.out")
            except Exception as e:
                print(e)
            fns = glob.glob(self.dirr+'/*.XV')
            for fn in fns:
                os.remove(fn)
            fns = glob.glob(self.dirr+'/*.ANI')
            for fn in fns:
                os.remove(fn)
            self.FDFgeom.writeFDF(self.dirr+"/STRUCT.fdf")
            SUR.MakePBS(None, self.dirr+"/RUN.pbs",\
                        [['$NODES$', '1:ppn=%i'%opts.proc]],\
                        True, rtype = 'TS')

    def shootOver(self, left, right):
        if self.ii == 0 or self.ii == opts.NNEB+1: return False

        xyz, F = N.array(self.XVgeom.xyz), N.array(self.forces)
        F = F[:, 1:4]
        lxyz, rxyz = N.array(left.xyz), N.array(right.xyz)

        tangent = rxyz-lxyz
        tangent = tangent/LA.norm(tangent) # Normalize

        F = F-N.sum(F*tangent)*tangent  # orthogonal to tangent

        # Apply constraints
        for s, f in self.const:
            for ii in range(s, f+1):
                F[ii, :] = 0
        if opts.const2 != None: # constraints 2
            indx, vec = opts.const2[1], opts.const2[2]
            F[indx,:] = N.dot(F[indx,:],vec)*vec # Allow along vec 
            indx, vec = opts.const2[3], opts.const2[4]
            F[indx,:] = F[indx,:]-N.dot(F[indx,:],vec)*vec # Plane perp to vec      

        return N.sum(F*self.v) < 0

    def update(self, left, right):
        if self.ii == 0 or self.ii == opts.NNEB+1: return
        # calculate new geometry
        xyz, F = N.array(self.XVgeom.xyz), N.array(self.forces)
        F = F[:, 1:4]
        lxyz, rxyz = N.array(left.xyz), N.array(right.xyz)

        tangent = rxyz-lxyz
        tangent = tangent/LA.norm(tangent) # Normalize

        # Calculate distance to move along tangent to reach midpoint
        dxyz = (rxyz+lxyz-2*xyz)/2.0*opts.tMix # Move to find mid-point /3 to converge
        dxyz = N.sum(dxyz*tangent)*tangent # Constraint to tangent
        steplen = N.sqrt(N.sum(dxyz*dxyz))
        if steplen > opts.maxDist:
            dxyz = dxyz/steplen*opts.maxDist
        
        F = F-N.sum(F*tangent)*tangent  # orthogonal to tangent

        # Apply constraints
        for s, f in self.const:
            for ii in range(s, f+1):
                F[ii, :] = 0
        if opts.const2 != None: # constraints 2
            indx, vec = opts.const2[1], opts.const2[2]
            F[indx,:] = N.dot(F[indx,:],vec)*vec # Allow along vec 
            indx, vec = opts.const2[3], opts.const2[4]
            F[indx,:] = F[indx,:]-N.dot(F[indx,:],vec)*vec # Plane perp to vec      
            
        self.F = F

        self.Fmax = N.max(N.sqrt(N.sum(F*F, 1)))
        self.converged = self.Fmax < opts.convCrit

        v = self.v + F*opts.moveK
        steplen = N.sqrt(N.sum(v*v))
        if steplen > opts.maxDist:
            v = v/steplen*opts.maxDist
        xyz += v

        if not self.converged:
            self.v = v
            self.FDFgeom.xyz = xyz
            self.XVgeom.xyz = xyz

        self.done = False

def checkConst(a,b):
    if not N.allclose(a.const, b.const):
        print("Error: NEB: constraints on initial and final states not the same")
        sys.exit(1)
    for ii in a.const:
        if not N.allclose(a.FDFgeom.xyz[ii[0]:ii[1]+1],b.FDFgeom.xyz[ii[0]:ii[1]+1]):
            sys.exit('Constrained atoms '+str([ii[0]+1,ii[1]+1])+'does not match for all geometries')
    if opts.const2 != None:
        # constraints 2
        indx, vec = opts.const2[1], opts.const2[2]
        dxyz = N.array(a.FDFgeom.xyz[indx])-N.array(b.FDFgeom.xyz[indx])
        if LA.norm(dxyz-N.dot(dxyz,vec)*vec) > 1e-6:
            sys.exit('Constraints not fullfilled for Atom %i along [%f,%f,%f]'%(indx+1,vec[0],vec[1],vec[2]))
        indx, vec = opts.const2[3], opts.const2[4]
        dxyz = N.array(a.FDFgeom.xyz[indx])-N.array(b.FDFgeom.xyz[indx])
        if N.abs(N.dot(dxyz,vec)) > 1e-6:
            sys.exit('Constraints not fullfilled for Atom %i in plane with tangent [%f,%f,%f]'%(indx+1,vec[0],vec[1],vec[2]))
    return

########################################################


def readxv(dirpath):
    global geom
    # Read geometry from first .XV file found in dirpath
    fns = glob.glob(dirpath+'/*.XV')

    if len(fns) > 1:
        print("ERROR: NEB: More than one .XV file in dir:%s"%dirpath)
        sys.exit(1)
    elif len(fns) < 1:
        return None

    print('Reading geometry from "%s" file' % fns[0])
    geom = MG.Geom(fns[0])
    return geom

###########################################################


def GetOptions(argv):
    usage = "usage: %prog [options] Initial Final"
    description = """
Nudged elastic band calculation script. The script will generate the intermediate (linear interpolation) steps in a NEB calculation, submit them to the que system and repeat the calculations until finding the approximate NEB solution.
Initial and Final states are read from the directories: Initial/CGrun, Final/CGrun. 
The directories should contain RUN.fdf (main fdf file), RUN.out (forces) and one .XV file containing the geometry. Constraints can be given as '%block GeometryConstraints' and/or using the Constraint option for small molecules. If the script has died, it will try to restart the calculations. It is usefull to use '-s' to make a better starting interpolation. [https://doi.org/10.1142/9789812839664_0016], note that the implementation along the trajectory novel ... if you obtain buckeling of the path, try to decrease 'acceleration' (-d), if the distance between points oscillate, decrease mixing (-t). Convergence does not demand equal spacing.

Intermediate steps will be written in:
NEB_.../CGrun directories

For help use --help!
"""

    import argparse

    # if text string is specified, convert to list
    print(argv)
    if isinstance(argv, VC.string_types):
        argv = argv.split()
    print(argv)
    p = argparse.ArgumentParser(description=description)
    p.add_argument('initial', help='Initial geometry')
    p.add_argument('final', help='Final geometry')

    p.add_argument("-n", "--NumNEB", dest="NNEB",
                   help="Number of intermediate steps [%(default)i]",
                   type=int, default=9)
    p.add_argument("-f", "--fdf", dest='fn', default='RUN.fdf', type=str,
                   help='Input fdf-file for SIESTA calculations [%(default)s]')
    p.add_argument("-p", "--Proc", dest="proc",
                   help="Number of processors [%(default)i]",
                   type=int, default=1)
    p.add_argument("-d", "--MoveK", dest="moveK",
                   help="Acceleration velocity/force/step [A/(eV/A)/step] [%(default)f]",
                   type=float, default=0.006)
    p.add_argument("-t", "--tangent", dest="tMix",
                   help="Mixing along the trajectory, i.e., move points to midpoint along tangent [%(default)f]",
                   type=float, default=.8)
    p.add_argument("-m", "--MaxDist", dest="maxDist",
                   help="Maximum distance moved (separately enforced along/perpendicular to tangent) [Ang] [%(default)f]",
                   type=float, default=0.04)
    p.add_argument("-e", "--ConvCrit", dest="convCrit",
                   help="Convergence criteria, max force on atom [eV/A] [%(default)f]",
                   type=float, default=0.08)
    p.add_argument("-s", "--Start", dest="onlyInit",
                   help="Only make the initial structures [%(default)s]",
                   default=False, action='store_true')
    p.add_argument("-c","--Constraint", dest="const2",
                   help="Constraints string 'n1, n2, x,y,z, n3, x,y,z': where atom n1 is fully constrained, n2 can move along the vector (x,y,z), and n3 in the plane perpendicular to the vector. [%(default)s]",
                   type=str, default=None)

    opts = p.parse_args()
    opts.module = 'NEB'

    # Parse constraints
    if opts.const2 != None:
        strings = string.split(opts.const2,',')
        opts.const2 = [int(strings[0]), int(strings[1]), N.array([float(strings[2]),float(strings[3]),float(strings[4])]),
                       int(strings[5]), N.array([float(strings[6]),float(strings[7]),float(strings[8])])]
        def normalize(x): 
            return x/LA.norm(x)
        opts.const2[2] = normalize(opts.const2[2])
        opts.const2[4] = normalize(opts.const2[4])

    print('##################################################################################')
    print('## NEB options')
    print(('Number of intermediate steps  : %i'%opts.NNEB))
    print(('Tangent mixing                : %f'%opts.tMix))
    print(('Processors                    : %i'%opts.proc))
    print(('Acceleration [A/(eV/A)/step]  : %f'%opts.moveK))
    print(('Max move distance/coord [A]   : %f'%opts.maxDist))
    print(('Convergence criteria [eV/A]   : %f'%opts.convCrit))
    print(('Constraints                   : ',opts.const2))
    print('##################################################################################')

    if opts.const2 != None: # Internal numbering
        opts.const2[0] -= 1
        opts.const2[1] -= 1
        opts.const2[3] -= 1

    return opts

################# Math helpers ################################
mm = MM.mm
outerAdd = MM.outerAdd
dist = MM.dist
mysqrt = MM.mysqrt
dagger = MM.dagger

