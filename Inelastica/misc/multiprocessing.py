"""

:mod:`Inelastica.misc.multiprocessing`
======================================

.. currentmodule:: Inelastica.misc.multiprocessing

"""

from __future__ import print_function, absolute_import

import sys
import multiprocessing as MP
import os


def runParallel(function, argList, nCPU=None):
    # Run in parallel the function with arguments given in the list
    # return list of results. You have to wrap the normal function with:
    # def myFuncPar(resQue, ii, *args):
    #     resQue.put( (ii,)+(myFunc(*args),))
    # Which returns the results of the arguments

    try: # Remove interfering OMP threading
        OMP = os.environ['OMP_NUM_THREADS']
    except:
        OMP = None
    try:
        OBLAS = os.environ['OPENBLAS_NUM_THREADS']
    except:
        OBLAS = None

    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    if nCPU == None:
        nCPU = MP.cpu_count()
    print("Running on %i CPUS"%(nCPU))

    resQue = MP.Queue() # return que
    chunks = [argList[ii*nCPU:(ii+1)*nCPU] for ii, jj in enumerate(argList[::nCPU])]
    res = [None]*len(argList)
    for ii, chunk in enumerate(chunks):
        threads = []
        for jj, args in enumerate(chunk):
            t = MP.Process(target=function, args=(resQue, ii*nCPU+jj,)+args)
            t.start()
            threads += [t]
        for jj in range(len(threads)):
            #print('Joining')
            out = resQue.get()
            #print('Joined',out[0])
            res[out[0]] = out[1]
            threads[out[0]-ii*nCPU].join()
            if threads[out[0]-ii*nCPU].exitcode > 0:
                sys.exit('Something wrong inside process ....')

    if OMP == None: # Reset threading
        del os.environ['OMP_NUM_THREADS']
    else:
        os.environ['OMP_NUM_THREADS'] = OMP
    if OBLAS == None:
        del os.environ['OPENBLAS_NUM_THREADS']
    else:
        os.environ['OPENBLAS_NUM_THREADS'] = OBLAS
    return res
