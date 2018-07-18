"""

:mod:`Inelastica.io.netcdf`
===========================

IO interface for handling netCDF4 files

.. currentmodule:: Inelastica.io.netcdf

Classes
-------

.. autosummary::
   :toctree:

   NCfile

"""
from __future__ import print_function

import numpy as N
import netCDF4 as NC4


def write(fn, array, label, vartype='d'):
    """
    The simplest possible way to write an array to
    a new or existing NetCDF file
    """
    nc = NCfile(fn)
    nc.write(array, label, vartype)
    nc.close()


class NCfile(object):

    def __init__(self, fn):
        "Returning instance of class."
        try:
            self.file = NC4.Dataset(fn, 'a')
        except:
            self.file = NC4.Dataset(fn, 'w')
        self.fn = fn
        self.dimensions = self.file.dimensions
        self.variables = self.file.variables
        self.invdim = {}
        for d in self.dimensions:
            val = self.dimensions[d]
            self.invdim[val] = d

    def close(self):
        "Closes the file instance"
        self.file.close()

    def write(self, A, label, vartype='d'):
        "Writes numpy array to file"
        dim = self.__checkDimensions(A)
        print('io.netcdf: Writing variable %s(%s) to file %s'%(label, vartype, self.fn))
        try:
            self.file.createVariable(label, vartype, dim)
        except:
            print('  ...variable %s already exist. Overwriting!!!'%label)
        self.variables[label][:] = N.array(A, dtype=vartype)

    def __checkDimensions(self, A):
        shape = N.shape(A)
        dim = []
        for i in shape:
            try:
                d = self.invdim[i]
            except:
                d = 'd%.2i'%len(self.dimensions)
                'io.netcdf: Generating dimension %s'%d
                self.file.createDimension(d, i)
                self.invdim[i] = d
            dim.append(d)
        return tuple(dim)


def test():
    # sample arrays
    a = N.ones((3, 3))
    b = N.diag([1, 2, 3])
    c = a-b
    d = N.ones((3, 3, 3))/10
    e = N.arange(10, dtype=N.float)/100.
    #
    fn = 'test.nc'
    write(fn, a, 'a')
    nc = NCfile(fn)
    nc.write(b, 'b')
    write(fn, c, 'c')
    write(fn, d, 'mytest')
    write(fn, e, 'energy')
    write(fn, c, 'He_ph', 'd')
    nc.close()


if __name__ == '__main__':
    test()
