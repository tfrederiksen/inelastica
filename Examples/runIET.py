#!/usr/bin/env python

import sys,os
from Inelastica import ValueCheck as VC
from Inelastica import Inelastica as I
from Inelastica import EigenChannels as E
from Inelastica import pyTBT as T

# Change the tolerance level of the imaginary part
VC.EditCheck("zero-imaginary-part",1e-7)

############################
#       Inelastica         #
############################
tmp_argv = ['-p','./PHrun/Dev_09-14.nc','-f','./TSrun/RUN.fdf']
# Copy sys argv to override
if len(sys.argv) > 1:
    tmp_argv.extend([a for a in sys.argv[1:]])

# Get option parser
opts = I.GetOptions(tmp_argv)

I.main(opts)


############################
#      EigenChannels       #
############################
tmp_argv = ['-w','cube','-f','./TSrun/RUN.fdf','-e','0.']
# Copy sys argv to override
if len(sys.argv) > 1:
    tmp_argv.extend([a for a in sys.argv[1:]])

# Get option parser
opts = E.GetOptions(tmp_argv)

for e in [-.1,0.,.1]:
    # Change energy...
    opts.energy = e
    opts.DestDir = 'Eig_'+str(e)
    # Redefine the piping output
    VC.CreatePipeOutput(opts.DestDir+'/Eig.log')
    E.main(opts)


############################
#          pyTBT           #
############################
tmp_argv = ['-f','./TSrun/RUN.fdf']
# Copy sys argv to override
if len(sys.argv) > 1:
    tmp_argv.extend([a for a in sys.argv[1:]])

# Get option parser
opts = T.GetOptions(tmp_argv,log='MyLog.log')

T.main(opts)

