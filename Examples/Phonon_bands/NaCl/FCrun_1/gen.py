import numpy as N
a= 5.64
a1=N.array([0,1,1],N.float)*a/2.0
a2=N.array([1,0,1],N.float)*a/2.0
a3=N.array([1,1,0],N.float)*a/2.0

for zz in range(6):
    for yy in range(6):
        for xx in range(6):
            tmp=xx*a1+yy*a2+zz*a3
            print "%f %f %f 1"%(tmp[0],tmp[1],tmp[2])
            print "%f %f %f 2"%(tmp[0],tmp[1],tmp[2]+a/2.0)

print a1*6
print a2*6
print a3*6
