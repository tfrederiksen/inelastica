"""
Example: How to use the SetupRuns functions to handle:
  Geometry optimizations   : CGrun
  Force constants (phonons): FCrun
  Overlap calculations     : OSrun
  TranSIESTA               : TSrun
  Phonons                  : PHrun
  Inelastica               : INrun
"""


from Inelastica.SetupRuns import *                     
from Inelastica.Phonons import *                       
import os                                   
import numpy as N

submitJob = False           # Automatically submit?

# Substitution rules for generating PBS scripts
TSPBSsubs = [['$NODES$','1:ppn=4'],['$MEM$','4gb'],['$WALLTIME$','100:00:00']]
PYPBSsubs = [['$NODES$','1:ppn=1'],['$MEM$','1gb'],['$WALLTIME$',  '5:00:00']]
OSPBSsubs = [['$NODES$','1:ppn=1'],['$MEM$','1gb'],['$WALLTIME$',  '1:00:00']]

# Choose which type of runs. For the example:
# 1: Run CG, wait for geometry to relax.
# 2: Start FC, OS, TS
# 3: Start PH when FC and OS are finnished
# 4: Try Inelastica and Eigenchannels when TS has finnished

T, F = True, False

# Step 1
CG = T

# Step 2
FC = F
OS = F
TS = F

# Step 3
PH = F

# Step 4
IN = F

if CG:
    # Stretch device and start geometry relaxation
    newL = 10.0  # New electrode separation
    SetupCGrun('./L9.68/CGrun','./L%.2f/CGrun'%newL,newL,
               AtomsPerLayer=1,
               overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs)                    

# For the remaining runs we will explore the L = 10.00 Ang case
geom = './L10.00'
head,tail = os.path.split(geom)   

if FC:
    for i in range(2):
        # Dynamic region (force constants)
        F,L = 11+i,11+i
        SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
                   FCfirst=F,FClast=L,  
                   overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs)                        

if OS:
    SetupOSrun(geom+'/CGrun',geom+'/OSrun',
               overwrite=False,submitJob=submitJob,PBSsubs=OSPBSsubs)                        

if TS:
    SetupTSrun(geom+'/CGrun','./L9.68/TSrun',geom+'/TSrun',
               AtomsPerLayer=1,BlockSize=2,
               LayersLeft=1,LayersRight=1,
               AddLeftList=[],AddRightList=[],
               overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs)

if PH:
    print "Try:"
    print "Phonons -c -F 9 -L 14 --FCfirst=11 --FClast=12 --EPHfirst=11 --EPHlast=12 PHrun"

if IN:
    print "Try:"
    print "Inelastica -F 10 -L 15 -p ../PHrun/Output.nc --LOEscale=0.0 INrun"

