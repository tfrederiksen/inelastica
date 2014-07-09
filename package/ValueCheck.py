"""
Contains general quantities that are used to do checks in the 
Inelastica package.

It allows users to edit the criterias of crashes without editing
the code.

Currently if this is to be used with Inelastica or pyTBT, it should
be in a wrapper call which emulates the options.
"""

import sys as s
import numpy as _np

# Place holder for checks
_check = {}

def _s(name,value): _check[name] = value

# The lower magnitude of the Hamiltonian/overlap elements between 
# the device region and electrode
_s("Device-Elec-overlap", 1e-7)

# The tolerance of difference in onlyS and FC displacements
_s("displacement-tolerance", 0.05)

# The tolerance for which the k-points are the same
_s("same-kpoint",1e-9)

# The tolerance when ensuring that the imaginary part is zero
_s("zero-imaginary-part",1e-8)
_s("trans-imaginary-part",1e-10)

# Convergence criteria for Lopez-Sancho algorithm
_s("Lopez-Sancho",1e-5)
_s("Lopez-Sancho-warning",1e-8)

del _s

def GetCheck(name):
    global _check
    return _check[name]

def EditCheck(name,value):
    """
    Sets the check of "name" to the value "value".
    It will notify the user so that any log files will retain
    this information.
    """
    global _check
    ov = "None"
    if name in _check:
        ov = _check[name]
    _check[name] = value
    print("WARNING: Overriding variable '{0}'".format(name))
    print("         Old value: {0}".format(ov))
    print("         New value: {0}".format(value))
    
def Check(name,val,*msgs):
    """
    Checks the value and exits if the value "val" is
    above the one specified
    """
    global _check
    if not name in _check:
        # the name hasn't been set
        # hence, it will always go through
        # This allows the coder to insets fully user tests 
        # Say if the code should generally never stop, we can make 
        # some checks where users might have certain criterias.
        return
    
    if _np.any(_np.array(val) > _check[name]):
        print("ERROR:")
        for msg in msgs:
            print(msg)
        print("       Check against '{0}' failed...".format(name))
        print("       Maximum value: {0}".format(_np.max(val)))
        print("       Error value  : {0}".format(_check[name]))
        raise ArithmeticError("Criteria not met. Please check output...")

def OptionsCheck(opts,exe):
    """
    Generic routine for adjusting most used options for routines.
    I.e. Inelastica/EigenChannels/pyTBT.
    """
    import SiestaIO as SIO
    import os, os.path as osp
        
    # Destination directory
    if not osp.isdir(opts.DestDir):
        print('\n'+exe+': Creating folder {0}'.format(opts.DestDir))
        os.mkdir(opts.DestDir)

    if not osp.isfile(opts.fn):
        raise IOError("FDF-file not found: "+opts.fn)

    # Read SIESTA files
    opts.head,tail = osp.split(opts.fn)
    print(exe+": Reading keywords from {0} \n".format(opts.fn))

    opts.systemlabel = SIO.GetFDFlineWithDefault(opts.fn,'SystemLabel', str, 'Systemlabel', exe) 
    opts.TSHS = '%s/%s.TSHS'%(opts.head,opts.systemlabel)

    # Electrodes
    opts.fnL  = opts.head+'/'+SIO.GetFDFlineWithDefault(opts.fn,'TS.HSFileLeft', str, None, exe)
    opts.NA1L = SIO.GetFDFlineWithDefault(opts.fn,'TS.ReplicateA1Left', int, 1, exe)
    opts.NA2L = SIO.GetFDFlineWithDefault(opts.fn,'TS.ReplicateA2Left', int, 1, exe)
    opts.fnR  = opts.head+'/'+SIO.GetFDFlineWithDefault(opts.fn,'TS.HSFileRight', str, None, exe)
    opts.NA1R = SIO.GetFDFlineWithDefault(opts.fn,'TS.ReplicateA1Right', int, 1, exe)
    opts.NA2R = SIO.GetFDFlineWithDefault(opts.fn,'TS.ReplicateA2Right', int, 1, exe)

    if opts.UseBulk < 0:
        opts.UseBulk = SIO.GetFDFlineWithDefault(opts.fn,'TS.UseBulkInElectrodes', bool, True, exe)

    # Read in number of buffer atoms
    opts.buffer,L,R = SIO.GetBufferAtomsList(opts.TSHS,opts.fn)
    opts.bufferL = L
    opts.bufferR = R

    if 'maxBias' in opts.__dict__:
        # Bias range
        opts.maxBias = abs(opts.maxBias)
        opts.minBias = -opts.maxBias

    # Device region
    if opts.DeviceFirst<=0:
        opts.DeviceFirst = SIO.GetFDFlineWithDefault(opts.fn,'TS.TBT.PDOSFrom',int,1,exe)
    opts.DeviceFirst -= L
    if opts.DeviceLast<=0:
        opts.DeviceLast = SIO.GetFDFlineWithDefault(opts.fn,'TS.TBT.PDOSTo',int,1e10,exe)
    opts.DeviceLast -= L
    opts.NumberOfAtoms = SIO.GetFDFlineWithDefault(opts.fn,'NumberOfAtoms',int,1e10,exe)
    opts.NumberOfAtoms -= L + R
    if opts.DeviceLast<opts.DeviceFirst:
        print(exe+' error: DeviceLast<DeviceFirst not allowed. Setting DeviceLast=DeviceFirst')
        opts.DeviceLast = opts.DeviceFirst
    opts.DeviceAtoms = [max(opts.DeviceFirst,1),min(opts.DeviceLast,opts.NumberOfAtoms)]

    # Voltage
    opts.voltage = SIO.GetFDFlineWithDefault(opts.fn,'TS.Voltage', float, 0.0, exe)

    #############
    # Here comes some specifics related to different executables:
    #############

    if "k1" in opts.__dict__:
        # Check k-point input
        opts.kPoint = _np.array([opts.k1,opts.k2,0.0],_np.float)

    if "VfracL" in opts.__dict__:
        if opts.VfracL < 0.0 or opts.VfracL > 1.0:
            raise RuntimeError('Option VfracL must be a value in the range [0,1].')

    if "PhononNetCDF" in opts.__dict__:
        if not osp.isfile(opts.PhononNetCDF):
            raise IOError("Electron-phonon coupling NetCDF file not found: "+opts.PhononNetCDF)
    if "eta" in opts.__dict__:
        if opts.eta < 0:
            raise RuntimeError("eta must be a posivite number")
    if "etaLead" in opts.__dict__:
        if opts.etaLead < 0:
            raise RuntimeError("etaLead must be a posivite number")
    if "PhExtDamp" in opts.__dict__:
        if opts.PhExtDamp < 0:
            raise RuntimeError("PhExtDamp must be a posivite number")
    if "biasPoints" in opts.__dict__:
        if opts.biasPoints < 6:
            raise AssertionError("BiasPoints must be larger than 5")
    if "iSpin" in opts.__dict__:
        if not opts.iSpin in [0,1]:
            raise AssertionError("Spin must be either 0 or 1")

    if "Emin" in opts.__dict__:
        if opts.Emin == 1e10:
            opts.Emin = SIO.GetFDFlineWithDefault(opts.fn,'TS.TBT.Emin', float, 0.0, 'pyTBT')
    if "Emax" in opts.__dict__:
        if opts.Emax == 1e10:
            opts.Emax = SIO.GetFDFlineWithDefault(opts.fn,'TS.TBT.Emax', float, 1.0, 'pyTBT')
    if "NPoints" in opts.__dict__:
        if opts.NPoints <= 0:
            opts.NPoints = SIO.GetFDFlineWithDefault(opts.fn,'TS.TBT.NPoints', int, 1, 'pyTBT')

        # Create list of energies
        opts.dE = (opts.Emax-opts.Emin)/float(opts.NPoints-1)
        opts.Elist = _np.arange(opts.Emin,opts.Emax,opts.dE)
        if opts.NPoints == 1:
            opts.dE = 0.
            opts.Elist = _np.array((opts.Emin,),_np.float)


def GetPositional(args,msg="You have not specified any positional argument"):
    if len(args) < 1:
        raise ValueError(msg)
    pos = args.pop(0)
    return pos

# Save the stdout pipes
_default_stdout = s.stdout
_default_stderr = s.stderr
def CreatePipeOutput(f):
    global _default_stdout, _default_stderr
    import subprocess, os, os.path as osp, errno

    # First ensure that the path to the file exists
    # In case one wishes to create a log folder this should
    # not be limited.
    head, tail = osp.split(f)
    try:
        # Create directory tree
        os.makedirs(head)
    except OSError as exc:
        if exc.errno == errno.EEXIST and osp.isdir(head):
            pass
        else: raise # forward error...

    class TeeLog(object):
        def __init__(self,f,term):
            self.term = term
            self.log = open(f,'w') # Consider doing this optionally appending?
        def write(self,message):
            self.term.write(message)
            self.log.write(message)
        def flush(self):
            self.term.flush()
            self.log.flush()

    # Overwrite the std-out and std-err
    s.stdout = TeeLog(f,_default_stdout)
    s.stderr = TeeLog(f,_default_stderr)

