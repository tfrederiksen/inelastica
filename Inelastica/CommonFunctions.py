"""

CommonFunctions (:mod:`Inelastica.CommonFunctions`)
===================================================

.. currentmodule:: Inelastica.CommonFunctions

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
    try:
        import Inelastica.info as info
        infover = info.label
    except Exception as e:
        print('Could not import info:')
        print(str(e))
        infover = ''

    print('=======================================================================')
    print('INELASTICA VERSION : %s'%(infover))
    print('RUNNING %s : %s'%(name.upper(), time.ctime()))
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

######################################################################
#
# Multiprocessing


def runParallel(function, argList, nCPU=None):
    # Run in parallel the function with arguments given in the list
    # return list of results. You have to wrap the normal function with:
    # def myFuncPar(resQue, ii, *args):
    #     resQue.put( (ii,)+(myFunc(*args),))
    # Which returns the results of the arguments

    import multiprocessing as MP
    import os

    try: # Remove interfering OMP threading
        OMP = os.environ['OMP_NUM_THREADS']
    except:
        OMP = None
    try:
        OBLAS = os.environ['OPENBLAS_NUM_THREADS']
    except:
        OBLAS = None

    os.environ['OMP_NUM_THREADS']='1'
    os.environ['OPENBLAS_NUM_THREADS']='1'
    if nCPU==None:
        nCPU=MP.cpu_count()
    print("Running on %i CPUS"%(nCPU))

    resQue = MP.Queue() # return que
    chunks = [argList[ii*nCPU:(ii+1)*nCPU] for ii, jj in enumerate(argList[::nCPU])]
    res = [None]*len(argList)
    for ii, chunk in enumerate(chunks):
        threads=[]
        for jj, args in enumerate(chunk):
            t= MP.Process(target=function, args =(resQue, ii*nCPU+jj,)+args)
            t.start()
            threads += [t]
        for jj in range(len(threads)):
            #print('Joining')
            out = resQue.get()
            #print('Joined',out[0])
            res[out[0]]=out[1]
            threads[out[0]-ii*nCPU].join()
            if threads[out[0]-ii*nCPU].exitcode>0:
                sys.exit('Something wrong inside process ....')

    if OMP==None: # Reset threading
        del os.environ['OMP_NUM_THREADS']
    else:
        os.environ['OMP_NUM_THREADS']=OMP
    if OBLAS==None:
        del os.environ['OPENBLAS_NUM_THREADS']
    else:
        os.environ['OPENBLAS_NUM_THREADS']=OBLAS
    return res
