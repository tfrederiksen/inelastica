import numpy as N
a= 4.23
a1=N.array([-1,1,1],N.float)*a/2.0
a2=N.array([1,-1,1],N.float)*a/2.0
a3=N.array([1,1,-1],N.float)*a/2.0

print a1*7
print a2*7
print a3*7

for zz in range(7):
    for yy in range(7):
        for xx in range(7):
            tmp=xx*a1+yy*a2+zz*a3
            print "%f %f %f 1"%(tmp[0],tmp[1],tmp[2])

