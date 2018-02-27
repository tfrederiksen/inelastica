"""

NEB (:mod:`Inelastica.NEB`)
===========================

Nudged elastic band

.. currentmodule:: Inelastica.NEB

classes
-------

.. autosummary::
   :toctree:

   general
   savedData
   step

"""

import Inelastica.io.siesta as SIO
import Inelastica.SetupRuns as SUR
import Inelastica.MakeGeom as MG
import Inelastica.MiscMath as MM
import Inelastica.physics.constants as PC
import Inelastica.io.xmgrace as XMGR

import numpy as N
import numpy.linalg as LA
import sys
import string
import struct
import glob
import os
import copy
import time
import pickle
from optparse import OptionParser, OptionGroup

##################### Global variabels #################


class general:
    pass
steps=[]


class savedData:
    pass

########################################################
##################### Main routine #####################
########################################################


def main():
    setupParameters()
    runNEB()


########################################################
def runNEB():
    global general, steps, savedData

    # Restarting?
    fns=glob.glob('NEB_*/CGrun')
    if len(fns)==0:
        restart = False
    elif len(fns)!=general.NNEB:
        sys.exit('Number of NEB_??? directories %i does not match number of steps %i'%(len(fns), general.NNEB))
    else:
        restart=True

    fns = ["NEB_%i/CGrun"%ii for ii in range(general.NNEB)]
    fns = [general.initial+"/CGrun"]+fns+[general.final+"/CGrun"]

    i, f = step(fns[0], False, 0), step(fns[general.NNEB+1], False, general.NNEB+1)
    checkConst(i, f)

    if restart:
        try:
            savedData.E, savedData.F, savedData.Fmax, savedData.geom, savedData.v=pickle.load(open('NEB_%i/savedData.pickle'%0, 'r'))
        except:
            savedData.E    = []
            savedData.F    = []
            savedData.geom = []
            savedData.Fmax = []
            savedData.v    = N.zeros((general.NNEB+2, len(i.XVgeom.xyz), 3), N.float)
    else:
        savedData.E    = []
        savedData.F    = []
        savedData.geom = []
        savedData.Fmax = []
        savedData.v    = N.zeros((general.NNEB+2, len(i.XVgeom.xyz), 3), N.float)

    steps = [step(fn, restart, ii, initial=i, final=f) for ii, fn in enumerate(fns)]

    done=False

    while not done:
        # Start  calculations
        for ii in steps:
            ii.run()

        nextstep=False
        while not nextstep:
            for ii in steps:
                if not ii.checkDone():
                    break
            else:
                nextstep=True
            if not nextstep:
                time.sleep(30)

        savedData.E += [[SIO.GetTotalEnergy(ii.dir+'/RUN.out') for ii in steps]]
        savedData.geom += [[copy.copy(ii.XVgeom) for ii in steps]]

        oldgeom = [copy.deepcopy(ii.XVgeom) for ii in steps]
        for ii, jj in enumerate(steps[1:-1]):
            jj.update(oldgeom[ii], oldgeom[ii+2])

        # Write status
        savedData.F += [[ii.Ftot for ii in steps]]
        savedData.Fmax += [[ii.Fmax for ii in steps]]
        for ii in range(1, len(steps)-1):
            geoms = [savedData.geom[jj][ii] for jj in range(len(savedData.geom))]
            Ftots = [savedData.Fmax[jj][ii] for jj in range(len(savedData.geom))]
            Fs    = [savedData.F[jj][ii] for jj in range(len(savedData.geom))]
            SIO.WriteANIFile('NEB_%i/Steps.ANI'%(ii-1), geoms, Ftots)
            SIO.WriteAXSFFiles('NEB_%i/Steps.XASF'%(ii-1), geoms, forces=Fs)

        geoms = [ii.XVgeom for ii in steps]
        Fmax  = [ii.Fmax for ii in steps]
        Ftot  = [ii.Ftot for ii in steps]
        SIO.WriteANIFile('NEB_%i/NextStep.ANI'%0, geoms, Fmax)
        SIO.WriteAXSFFiles('NEB_%i/NextStep.XASF'%0, geoms, forces=Ftot)
        done = True
        for ii in steps:
            done = done and ii.converged
        pickle.dump((savedData.E, savedData.F, savedData.Fmax, savedData.geom, savedData.v),\
                    open('NEB_%i/savedData.pickle'%0, 'w'))

#################### Class for each step ###############


class step:
    global steps, general

    def __init__(self, dir, restart, iistep, initial=None, final=None):
        self.dir=dir
        self.converged, self.Fmax = False, 0.0
        self.ii=iistep
        self.fixed = (iistep==0) or (iistep==general.NNEB+1)

        if not restart and not self.fixed:
            os.makedirs(dir)
            SUR.CopyInputFiles(general.initial+"/CGrun/", dir,\
                                   ['.fdf', '.vps', '.psf'])
            # Interpolate
            ixyz, fxyz = N.array(initial.XVgeom.xyz), N.array(final.XVgeom.xyz)
            mix = float(iistep)/(general.NNEB+1.0)
            xyz = (1-mix)*ixyz+mix*fxyz
            self.FDFgeom = copy.copy(initial.XVgeom)
            self.FDFgeom.xyz = [xyz[ii, :] for ii in range(len(xyz))]
            self.FDFgeom.writeFDF(self.dir+"/STRUCT.fdf")

            # Append lines to RUN.fdf
            elm = dir+"/RUN.fdf"
            f = open(elm, 'r')
            lines = f.readlines()
            f.close()
            f = open(elm, 'w')
            f.write('### Lines appended %s \n' %time.ctime())
            f.write("SolutionMethod       diagon\n")
            f.write("MD.TypeOfRun         CG\n")
            f.write("MD.NumCGsteps        0\n")
            f.write("TS.SaveHS            True\n")
            f.write("MD.UseSaveXV         False\n")
            f.write("MD.UseSaveCG         False\n")
            f.write('# end of lines appended %s \n' %time.ctime())
            for line in lines: f.write(line)
            f.close()

        self.done = self.checkDone()
        self.const = SIO.GetFDFblock(dir+"/RUN.fdf", "GeometryConstraints")

    def checkDone(self):
        if self.fixed==True:
            self.FDFgeom = MG.Geom(self.dir+"/RUN.fdf")
            self.XVgeom  = readxv(self.dir)
            self.forces  = SIO.ReadForces(self.dir+"/RUN.out")
            self.forces  = self.forces[-len(self.XVgeom.xyz):]
            self.Ftot    = self.forces
            self.converged = True
            return True
        else:
            SIO.CheckTermination(self.dir+"/RUN.out")
            self.FDFgeom = MG.Geom(self.dir+"/RUN.fdf")
            done = False
            try:
                self.XVgeom  = readxv(self.dir)
                self.forces  = SIO.ReadForces(self.dir+"/RUN.out")
                print len(self.forces[0])
                if N.allclose(self.XVgeom.xyz, self.FDFgeom.xyz, 1e-6):
                    done=True
            except:
                done=False
            return done

    def run(self):
        if (not self.done) and (not self.converged):
            try:
                os.remove(self.dir+"/RUN.out")
            except:
                pass
            fns=glob.glob(self.dir+'/*.XV')
            for fn in fns: os.remove(fn)
            fns=glob.glob(self.dir+'/*.ANI')
            for fn in fns: os.remove(fn)
            self.FDFgeom.writeFDF(self.dir+"/STRUCT.fdf")
            SUR.MakePBS(None, self.dir+"/RUN.pbs",\
                            [['$NODES$', '1:ppn=%i'%general.proc]],\
                            True, type = 'TS')

    def update(self, left, right):
        # calculate new geometry
        xyz, F = N.array(self.XVgeom.xyz), N.array(self.forces)
        if len(F[0])!=4:
            print F
            print F.shape
        F = F[:, 1:4]
        lxyz, rxyz = N.array(left.xyz), N.array(right.xyz)

        tangent = rxyz-lxyz
        tangent = tangent/N.sqrt(N.sum(tangent*tangent)) # Normalize
        FS = general.SK*(rxyz+lxyz-2*xyz) # Spring forces
        FS = N.sum(FS*tangent)*tangent    # Allong tangent
        F  = F-N.sum(F*tangent)*tangent  # orthogonal to tangent
        Ftot = F+FS

        # Apply constraints
        for s, f in general.const:
            for ii in range(s, f+1):
                Ftot[ii, :]=0

        self.Ftot=Ftot
        savedData.v[self.ii, :, :]=savedData.v[self.ii, :, :]-N.sum(savedData.v[self.ii, :, :]*tangent)*tangent
        if N.sum(Ftot*savedData.v[self.ii, :, :])<0:
            savedData.v[self.ii, :, :]=0
        savedData.v[self.ii, :, :]=savedData.v[self.ii, :, :]+Ftot*general.moveK
        steplen=N.sqrt(N.sum(savedData.v[self.ii, :, :]*savedData.v[self.ii, :, :]))
        if steplen>general.maxDist:
            savedData.v[self.ii, :, :]=savedData.v[self.ii, :, :]/steplen*general.maxDist
        xyz = xyz+savedData.v[self.ii, :, :]
        xyz=[xyz[ii, :] for ii in range(len(xyz))]
        self.FDFgeom.xyz = xyz
        self.XVgeom.xyz = xyz

        self.Fmax = N.max(N.sqrt(N.sum(Ftot*Ftot, 1)))
        self.converged = self.Fmax<general.convCrit
        self.done = False


def checkConst(initial, final):
    if initial.const!=final.const:
        print "Error: NEB: constraints on initial and final states not the same"
        sys.exit(1)
    general.const = []
    for ii in initial.const:
        general.const+=[[int(ii[2])-1, int(ii[4])-1]]
    print "Constraints on atoms:"
    for ii in general.const:
        print "from %i to %i (Siesta numbering)"%(ii[0]+1, ii[1]+1)

########################################################


def readxv(dir):
    global geom
    # Read geometry from first .XV file found in dir
    fns=glob.glob(dir+'/*.XV')

    if len(fns)>1:
        print "ERROR: NEB: More than one .XV file in dir:%s"%dir
        sys.exit(1)
    elif len(fns)<1:
        return None

    print('Reading geometry from "%s" file' % fns[0])
    geom = MG.Geom(fns[0])
    return geom

###########################################################


def setupParameters():
    global general

    usage = "usage: %prog [options] Initial Final"
    description = """
Nudged elastic band calculation.
Initial and Final and final states are read from the directories:
Initial/CGrun
Final/CGrun
The directories should contain RUN.fdf (main fdf file) and STRUCT.fdf containing the structure of the the system (read from RUN.fdf with %include). RUN.out should contain the output of the Siesta calculations.

Intermediate steps will be written in:
NEB_.../CGrun directories

For help use --help!
"""
    parser = OptionParser(usage, description=description)
    EC = OptionGroup(parser, "Options for NEB")
    EC.add_option("-n", "--NumNEB", dest="NNEB",\
                      help="Number of intermediate steps [%default]",\
                      type='int', default=10)
    EC.add_option("-k", "--SpringK", dest="SK",\
                      help="Spring constant [eV/A] [%default]",\
                      type='float', default=100.0)
    EC.add_option("-p", "--Proc", dest="proc",\
                      help="Number of processors [%default]",\
                      type='int', default=1)
    EC.add_option("-d", "--MoveK", dest="moveK",\
                      help="Distance moved/force [A/(eV/A)] [%default]",\
                      type='float', default=0.004)
    EC.add_option("-m", "--MaxDist", dest="maxDist",\
                      help="Maximum distance moved [A/direction] [%default]",\
                      type='float', default=0.04)
    EC.add_option("-c", "--ConvCrit", dest="convCrit",\
                      help="Convergence criteria, max force on atom [eV/A] [%default]",\
                      type='float', default=0.08)

    parser.add_option_group(EC)

    (general, args) = parser.parse_args()
    print description

    print args
    if len(args)<2:
        parser.error("ERROR: You need to specify initial and final geometries")
    if len(args)>2:
        parser.error("ERROR: you have more than 2 geometries")

    general.initial = args[0]
    general.final = args[1]

    class myopen:
        # Double stdout to RUN.out and stdout
        def write(self, x):
            self.stdout.write(x)
            self.file.write(x)

    fo = myopen()
    fo.stdout, fo.file = sys.stdout, open('NEB.out', 'w', 0)
    sys.stdout = fo

    argv=""
    for ii in sys.argv: argv+=" "+ii
    print(argv)
    print('##################################################################################')
    print('## NEB options')
    print('Number of intermediate steps  : %i'%general.NNEB)
    print('Spring constant [eV/A^2]      : %f'%general.SK)
    print('Processor                     : %i'%general.proc)
    print('Move constant [A/(eV/A)]      : %f'%general.moveK)
    print('Max move distance/coord [A]   : %f'%general.maxDist)
    print('Convergence criteria [eV/A]   : %f'%general.convCrit)
    print('##################################################################################')

################# Math helpers ################################
mm = MM.mm
outerAdd = MM.outerAdd
dist = MM.dist
mysqrt = MM.mysqrt
dagger = MM.dagger

##################### Start main routine #####################

if __name__ == '__main__':
    main()
