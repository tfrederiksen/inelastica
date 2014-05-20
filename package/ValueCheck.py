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
        raise ValueError("Criteria not met. Please check output...")
