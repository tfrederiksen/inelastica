import sys, profile
sys.path+=['../..']

import Inelastica.SiestaIO as SIO
import numpy as N
import numpy.random as RA

def main():
    buf = SIO.GetBufferAtomsList('sample.onlyS','test_buf1.fdf')
    a = N.arange(1,4)
    a = N.append(a,N.arange(41,45))
    if N.any(a-buf != 0):
        print buf
        print a
        raise ValueError('Error on reading in buffer atoms from test_buf1.fdf')
    buf = SIO.GetBufferAtomsList('sample.onlyS','test_buf2.fdf')
    a = N.arange(1,4)
    a = N.append(a,N.arange(4,6))
    a = N.append(a,N.arange(42,45))
    if N.any(a-buf != 0):
        print buf
        print a
        raise ValueError('Error on reading in buffer atoms from test_buf2.fdf')

    try:
        buf = SIO.GetBufferAtomsList('sample.onlyS','test_buf3.fdf')
        raise ValueError('Error in algorithm. This list is non-consecutive')
    except:
        pass # Test passed


main()
