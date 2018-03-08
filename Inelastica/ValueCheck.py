"""

ValueCheck (:mod:`Inelastica.ValueCheck`)
=========================================

Contains general quantities that are used to do checks in the
`Inelastica`_ package.

It allows users to edit the criterias of crashes without editing
the code.

Currently if this is to be used with `Inelastica`_ or `pyTBT`, it should
be in a wrapper call which emulates the options.

.. currentmodule:: Inelastica.ValueCheck


"""

# Simple case for allowing the newer print function (Python3 compliant)
from __future__ import print_function

import sys as s
import numpy as _np


# Python 3-version
# Check of strings
if s.version_info[0] == 3:
    string_types = str,
else:
    string_types = basestring,

# Place holder for checks
_check = {}


def _s(name, value): _check[name] = value

# The lower magnitude of the Hamiltonian/overlap elements between
# the device region and electrode
_s("Device-Elec-overlap", 1e-7)

# The tolerance of difference in onlyS and FC displacements
_s("displacement-tolerance", 0.05)

# The tolerance for which the k-points are the same
_s("same-kpoint", 1e-9)

# The tolerance when ensuring that the imaginary part is zero
_s("zero-imaginary-part", 1e-8)
_s("trans-imaginary-part", 1e-10)

# Convergence criteria for Lopez-Sancho algorithm
_s("Lopez-Sancho", 1e-5)
_s("Lopez-Sancho-warning", 1e-8)

del _s


def GetCheck(name):
    global _check
    return _check[name]


def EditCheck(name, value):
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


def Check(name, val, *msgs):
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


def OptionsCheck(opts, exe):
    """
    Generic routine for adjusting most used options for routines.
    I.e. Inelastica/EigenChannels/pyTBT.
    """
    import Inelastica.io.siesta as SIO
    import os
    import os.path as osp

    # Destination directory
    if not osp.isdir(opts.DestDir):
        print('\n'+exe+': Creating folder {0}'.format(opts.DestDir))
        os.mkdir(opts.DestDir)

    if not osp.isfile(opts.fn):
        raise IOError("FDF-file not found: "+opts.fn)

    # Read SIESTA files
    opts.head = osp.split(opts.fn)[0]
    if opts.head == '': # set filepath if missing
        opts.head = '.'
    print(exe+": Reading keywords from {0} \n".format(opts.fn))

    opts.systemlabel = SIO.GetFDFlineWithDefault(opts.fn, 'SystemLabel', str, 'siesta', exe)
    opts.TSHS = '%s/%s.TSHS'%(opts.head, opts.systemlabel)

    # These first keys can be used, but they are superseeded by keys in the TS.Elec.<> block
    # Hence if they are read in first it will do it in correct order.

    if opts.UseBulk < 0:
        # Note NRP:
        #  in principle this is now a per-electrode setting which
        #  may be useful for certain systems...
        opts.UseBulk = SIO.GetFDFlineWithDefault(opts.fn, 'TS.UseBulkInElectrodes', bool, True, exe)
        opts.UseBulk = SIO.GetFDFlineWithDefault(opts.fn, 'TS.Elecs.Bulk', bool, opts.UseBulk, exe)

    def get_elec_vars(lr):

        # Look up old format first
        TSHS = SIO.GetFDFlineWithDefault(opts.fn, 'TS.HSFile'+lr, str, '', exe)
        NA1 = SIO.GetFDFlineWithDefault(opts.fn, 'TS.ReplicateA1'+lr, int, 1, exe)
        NA2 = SIO.GetFDFlineWithDefault(opts.fn, 'TS.ReplicateA2'+lr, int, 1, exe)

        # default semi-inf direction
        semiinf = 2

        # Proceed looking up new format, which precedes
        belec = 'TS.Elec.' + lr
        print('Looking for new electrode format in: %%block {}'.format(belec))

        # Default replication stuff
        TSHS = SIO.GetFDFlineWithDefault(opts.fn, belec+'.TSHS', str, TSHS, exe)
        NA1 = SIO.GetFDFlineWithDefault(opts.fn, belec+'.Bloch.A1', int, NA1, exe)
        NA2 = SIO.GetFDFlineWithDefault(opts.fn, belec+'.Bloch.A2', int, NA2, exe)
        NA3 = SIO.GetFDFlineWithDefault(opts.fn, belec+'.Bloch.A3', int, 1, exe)

        # Overwrite block
        block = SIO.GetFDFblock(opts.fn, KeyWord=belec)

        for line in block:
            print(line)
            # Lower-case, FDF is case-insensitive
            key = line[0].lower()
            if key in ['tshs', 'tshs-file', 'hs', 'hs-file']:
                TSHS = line[1]
            elif key in ['replicate-a', 'rep-a', 'replicate-a1', 'rep-a1', 'bloch-a1']:
                NA1 = int(line[1])
            elif key in ['replicate-b', 'rep-b', 'replicate-a2', 'rep-a2', 'bloch-a2']:
                NA2 = int(line[1])
            elif key in ['replicate-c', 'rep-c', 'replicate-a3', 'rep-a3', 'bloch-a3']:
                NA3 = int(line[1])
            elif key in ['replicate', 'rep', 'bloch']:
                # We have *at least* 2 integers
                NA1 = int(line[1])
                NA2 = int(line[2])
                NA3 = int(line[3])

            elif key in ['semi-inf-direction', 'semi-inf-dir', 'semi-inf']:
                # This is lower-case checked
                axis = line[1][1:].lower()
                if 'a' == axis or 'a1' == axis:
                    semiinf = 0
                elif 'b' == axis or 'a2' == axis:
                    semiinf = 1
                elif 'c' == axis or 'a3' == axis:
                    semiinf = 2

                # Simple check of input, this may be overwritten later
                if semiinf != 2 and NA1 * NA2 * NA3 > 1:
                    raise ValueError(("Cannot provide semi-infinite directions "
                                      "other than A3-direction with repetition."))
        if TSHS[0] != '/':
            # path is relative
            TSHS = opts.head+'/'+TSHS
        return TSHS, NA1, NA2, semiinf

    # Look up electrode block
    block = SIO.GetFDFblock(opts.fn, KeyWord='TS.Elecs')
    if len(block) == 0:
        # Did not find the electrode block, defaults to old naming scheme
        opts.fnL, opts.NA1L, opts.NA2L, opts.semiinfL = get_elec_vars('Left')
        opts.fnR, opts.NA1R, opts.NA2R, opts.semiinfR = get_elec_vars('Right')
    elif len(block) == 2:
        # NB: The following assumes that the left electrode is the first in the block!
        opts.fnL, opts.NA1L, opts.NA2L, opts.semiinfL = get_elec_vars(block[0][0])
        opts.fnR, opts.NA1R, opts.NA2R, opts.semiinfR = get_elec_vars(block[1][0])
    else:
        print(block)
        raise IOError('Currently only two electrodes are supported')

    # Read in number of buffer atoms
    opts.buffer, L, R = SIO.GetBufferAtomsList(opts.TSHS, opts.fn)
    opts.bufferL = L
    opts.bufferR = R

    if 'maxBias' in opts.__dict__:
        # Bias range
        opts.maxBias = abs(opts.maxBias)
        opts.minBias = -abs(opts.maxBias)

    # Device region
    if opts.DeviceFirst <= 0:
        opts.DeviceFirst = SIO.GetFDFlineWithDefault(opts.fn, 'TS.TBT.PDOSFrom', int, 1, exe)
    opts.DeviceFirst -= L
    if opts.DeviceLast <= 0:
        opts.DeviceLast = SIO.GetFDFlineWithDefault(opts.fn, 'TS.TBT.PDOSTo', int, 1e10, exe)
    opts.DeviceLast -= L
    opts.NumberOfAtoms = SIO.GetFDFlineWithDefault(opts.fn, 'NumberOfAtoms', int, 1e10, exe)
    opts.NumberOfAtoms -= L + R
    if opts.DeviceLast < opts.DeviceFirst:
        print(exe+' error: DeviceLast<DeviceFirst not allowed. Setting DeviceLast=DeviceFirst')
        opts.DeviceLast = opts.DeviceFirst
    opts.DeviceAtoms = [max(opts.DeviceFirst, 1), min(opts.DeviceLast, opts.NumberOfAtoms)]

    # Voltage
    opts.voltage = SIO.GetFDFlineWithDefault(opts.fn, 'TS.Voltage', float, 0.0, exe)

    #############
    # Here comes some specifics related to different executables:
    #############

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
        if not opts.iSpin in [0, 1]:
            raise AssertionError("Spin must be either 0 or 1")

    if "Emin" in opts.__dict__:
        if opts.Emin == 1e10:
            opts.Emin = SIO.GetFDFlineWithDefault(opts.fn, 'TS.TBT.Emin', float, 0.0, 'pyTBT')
    if "Emax" in opts.__dict__:
        if opts.Emax == 1e10:
            opts.Emax = SIO.GetFDFlineWithDefault(opts.fn, 'TS.TBT.Emax', float, 1.0, 'pyTBT')
    if "NPoints" in opts.__dict__:
        if opts.NPoints <= 0:
            opts.NPoints = SIO.GetFDFlineWithDefault(opts.fn, 'TS.TBT.NPoints', int, 1, 'pyTBT')

        # Create list of energies
        try:
            opts.Elist
            # Do not overwrite if some Elist is already specified
        except:
            if opts.NPoints == 1:
                opts.Elist = _np.array((opts.Emin,), _np.float)
            else:
                # Linspace is just what we need
                opts.Elist = _np.linspace(opts.Emin, opts.Emax, opts.NPoints)


def GetPositional(args, msg="You have not specified any positional argument"):
    if len(args) < 1:
        raise ValueError(msg)
    pos = args.pop(0)
    return pos
