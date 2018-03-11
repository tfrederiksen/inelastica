"""

:mod:`Inelastica.io.log`
========================

This module provides the functions to generate logfiles in the `Inelastica` package.

.. currentmodule:: Inelastica.io.log

"""

from __future__ import print_function, absolute_import

import sys
import time

# Save the stdout pipes
_default_stdout = sys.stdout
_default_stderr = sys.stderr


def CreatePipeOutput(f):
    global _default_stdout, _default_stderr
    import os
    import os.path as osp
    import errno

    # First ensure that the path to the file exists
    # In case one wishes to create a log folder this should
    # not be limited.
    head = osp.split(f)[0]
    try:
        # Create directory tree
        os.makedirs(head)
    except OSError as exc:
        if exc.errno == errno.EEXIST and osp.isdir(head):
            pass
        else: raise # forward error...

    class TeeLog(object):

        def __init__(self, f, term):
            self.term = term
            self.log = open(f, 'w') # Consider doing this optionally appending?

        def write(self, message):
            self.term.write(message)
            self.log.write(message)

        def flush(self):
            self.term.flush()
            self.log.flush()

    # Overwrite the std-out and std-err
    sys.stdout = TeeLog(f, _default_stdout)
    sys.stderr = TeeLog(f, _default_stderr)


def PrintMainHeader(name, options):
    import Inelastica.info as info
    print('=======================================================================')
    print('INELASTICA VERSION : %s'%(info.version))
    if info.git_count > 0:
        print('               GIT : %s'%(info.git_revision))
    print('RUNNING %s : %s'%(name.upper(), time.ctime()))
    if options:
        print('\nOPTIONS :')
        opts_dict = vars(options)
        keys = sorted(opts_dict)
        for i in keys:
            print('    ', i, '-->', opts_dict[i])
    print('=======================================================================')


def PrintMainFooter(name):
    print('=======================================================================')
    print('FINISHED %s : %s'%(name.upper(), time.ctime()))
    print('=======================================================================')


def PrintScriptSummary(argv, dT):
    print('SCRIPT SUMMARY:')

    # Write function call
    print('Call:', ' '.join(argv))

    # Timing
    hours = dT.days/24.+(dT.seconds+dT.microseconds*1.e-6)/60.**2
    minutes = hours*60.
    seconds = minutes*60.
    print('Program finished:  %s '%time.ctime())
    print('Walltime: %.2f hrs = %.2f min = %.2f sec'%(hours, minutes, seconds))
    print('=======================================================================')