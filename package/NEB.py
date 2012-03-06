print "SVN $Id$"

"""
Nudged elastic band
"""
try: 
    import Inelastica.SiestaIO as SIO
    import Inelastica.SetupRuns as SUR
    import Inelastica.MakeGeom as MG
    import Inelastica.MiscMath as MM
    import Inelastica.PhysicalConstants as PC
except:
    import SiestaIO as SIO
    import SetupRuns as SUR
    import MakeGeom as MG
    import MiscMath as MM
    import PhysicalConstants as PC

import numpy as N
import numpy.linalg as LA
import sys, string, struct, glob, os, copy
from optparse import OptionParser, OptionGroup

##################### Global variabels #################
class general:
    pass
steps=[]


########################################################
##################### Main routine #####################
########################################################
def main():
    setupParameters()
    runNEB()


########################################################
def runNEB():
    global general, steps
    
    # Restarting?
    fns=glob.glob('NEB_*/CGrun')
    if len(fns)==0:
        restart = False
    elif len(fns)!=general.NNEB:
        sys.exit('Number of NEB_??? directories %i does not match number of steps %i'%(len(fns),general.NNEB))
    else:
        restart=True

    fns = ["NEB_%i/CGrun"%ii for ii in range(general.NNEB)]
    fns = [general.initial+"/CGrun"]+fns+[general.final+"/CGrun"]
    
    i, f = step(fns[0],False,0), step(fns[general.NNEB+1],False,general.NNEB+1)
    checkConst(i,f)

    steps = [step(fn,restart,ii,initial=i,final=f) for ii,fn in enumerate(fns)]

    # Start initial calculation
    for ii in steps:
        ii.run()

#################### Class for each step ###############
class step:
    global steps, general
    def __init__(self,dir,restart,iistep,initial=None,final=None):
        self.dir=dir
        self.ii=iistep
        self.fixed = (iistep==0) or (iistep==general.NNEB+1)
        
        if not restart and not self.fixed:
            os.makedirs(dir)
            SUR.CopyInputFiles(general.initial+"/CGrun/",dir,\
                                   ['.fdf','.vps','.psf'])
            # Interpolate
            ixyz, fxyz = N.array(initial.XVgeom.xyz), N.array(final.XVgeom.xyz)
            mix = float(iistep)/(general.NNEB+1.0)
            xyz = mix*ixyz+(1.0-mix)*fxyz
            self.FDFgeom = copy.deepcopy(initial.XVgeom)
            self.FDFgeom.xyz = [xyz[ii,:] for ii in range(len(xyz))]
            self.FDFgeom.writeFDF(dir+"/STRUCT.fdf")

        self.done=False or self.fixed
        
        self.FDFgeom = MG.Geom(dir+"/RUN.fdf")
        try:
            self.XVgeom  = readxv(dir)
            self.forces  = SIO.ReadForces(dir+"/RUN.out")
            if N.allclose(self.XVgeom.xyz,self.FDFgeom.xyz,1e-6):
                self.done=True
        except:
            pass
        self.const = SIO.GetFDFblock(dir+"/RUN.fdf","GeometryConstraints")

    def run(self):
        if not self.done:
            SUR.MakePBS(None, self.dir+"/RUN.pbs",\
                            [['$NODES$','1:ppn=%i'%general.proc]],\
                            True, type = 'TS')

def checkConst(initial,final):
    if initial.const!=final.const:
        print "Error: NEB: constraints on initial and final states not the same"
        sys.exit(1)
    general.const = []
    for ii in initial.const:
        general.const+=[[int(ii[2])-1,int(ii[4])-1]]
    print "Constraints on atoms:"
    for ii in general.const:
        print "from %i to %i (Siesta numbering)"%(ii[0]+1,ii[1]+1)

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
The directories should contain RUN.fdf (main fdf file) and STRUCT.fdf containing the structure of the the system (read from RUN.fdf with %include).

Intermediate steps will be written in:
NEB_.../CGrun directories

For help use --help!
"""
    parser = OptionParser(usage,description=description)
    EC = OptionGroup(parser, "Options for NEB")
    EC.add_option("-n", "--NumNEB", dest="NNEB",\
                      help="Number of intermediate steps [%default]",\
                      type='int', default=10)
    EC.add_option("-k", "--SpringK", dest="SK",\
                      help="Spring constant [eV/A] [%default]",\
                      type='float', default=1.0)
    EC.add_option("-p", "--Proc", dest="proc",\
                      help="Number of processors [%default]",\
                      type='int', default=1)
    parser.add_option_group(EC)
    
    (general, args) = parser.parse_args()
    print description

    if len(sys.argv)<3:
        parser.error("ERROR: You need to specify initial and final geometries")
    if len(sys.argv)>3:
        parser.error("ERROR: you have more than 2 geometries")

    general.initial = sys.argv[1]
    general.final = sys.argv[2]

    class myopen:
        # Double stdout to RUN.out and stdout
        def write(self,x):
            self.stdout.write(x)
            self.file.write(x)

    fo = myopen()
    fo.stdout, fo.file = sys.stdout, open('NEB.out','w',0)
    sys.stdout = fo

    argv=""
    for ii in sys.argv: argv+=" "+ii
    print(argv)
    print('##################################################################################')
    print('## NEB options')
    print('Number of intermediate steps  : %i'%general.NNEB)
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

