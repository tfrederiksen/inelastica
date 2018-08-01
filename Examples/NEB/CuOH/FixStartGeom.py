from __future__ import print_function

import glob
import string
import os 
import Inelastica.MakeGeom as MG
import numpy as N
import numpy.linalg as LA

# Find geometries
fns = glob.glob('NEB*/CGrun/RUN.fdf')

# Fix order
indx = [int(string.split(string.split(ii, '_')[1], '/')[0]) for ii in fns]
fns = [fns[ii] for ii in N.argsort(indx)]


geoms = [MG.Geom(ii) for ii in fns]
xyz = [N.array(ii.xyz) for ii in geoms]

# Find bond length
d = LA.norm(xyz[0][37]-xyz[0][36])

# Move in z to keep bond length
for ii in xyz:
    ii[37, 2] = ii[36, 2]+N.sqrt(d**2-(ii[36, 0]-ii[37, 0])**2-(ii[36, 1]-ii[37, 1])**2)

# Redistribute points

xyz, nxyz = N.array(xyz), N.array(xyz) # Two copies

def redist(L,C,R):
    t = (R-L)/LA.norm(L-R)
    dx = ((R-C)-(C-L))/2 # Shift to midpoint 
    dx = N.sum(dx*t)*t # Project along tangent
    return dx

print([LA.norm(xyz[ii]-xyz[ii+1]) for ii in range(len(fns)-1)])

for ii in range(100):
    for jj in range(1, len(fns)-1):
        nxyz[jj] += 0.8*redist(xyz[jj-1], xyz[jj], xyz[jj+1])
    xyz = 1.0*nxyz

# Write output
for ii, fn in enumerate(fns):
    geoms[ii].xyz = xyz[ii]
    geoms[ii].writeFDF(os.path.dirname(fn)+'/STRUCT.fdf')
