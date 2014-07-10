import sys, profile
sys.path+=['../..']

import Inelastica.SiestaIO as SIO
import numpy as N
import numpy.random as RA

def main():
    # test 1
    buf = SIO.GetBufferAtomsList('sample.onlyS','test_buf1.fdf')
    a = N.arange(1,4)
    a = N.append(a,N.arange(41,45))
    print 'buf =',buf
    print 'a =  ',a
    if N.any(a-buf[0] != 0):
        raise ValueError('Error on reading in buffer atoms from test_buf1.fdf')
    else:
        print 'Test 1 passed'
    # test 2
    buf = SIO.GetBufferAtomsList('sample.onlyS','test_buf2.fdf')
    a = N.arange(1,4)
    a = N.append(a,N.arange(4,6))
    a = N.append(a,N.arange(42,45))
    print 'buf =',buf
    print 'a =  ',a
    if N.any(a-buf[0] != 0):
        raise ValueError('Error on reading in buffer atoms from test_buf2.fdf')
    else:
        print 'Test 2 passed'
    # test 3
    try:
        buf = SIO.GetBufferAtomsList('sample.onlyS','test_buf3.fdf')
        raise ValueError('Error in algorithm. This list is non-consecutive')
    except:
        print 'Test 3 passed'

main()
