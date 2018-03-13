"""
Example: How to use the SetupRuns functions to handle:
  Force constants (phonons): FCrun
  Overlap calculations     : OSrun
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

FC = T
OS = T



if FC:
    # Dynamic region (force constants) SiC part1
    F,L = 46,51
    SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
               FCfirst=F,FClast=L,  
               overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs)     
    # Dynamic region (force constants) SiC part2
    F,L = 52,57
    SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
               FCfirst=F,FClast=L,  
               overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs) 
    # Dynamic region (force constants) Buffer 
    F,L = 58,65  
    SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
               FCfirst=F,FClast=L,  
               overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs)   
    
    # Dynamic region (force constants) graphene
    F,L = 66, 73
    SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
               FCfirst=F,FClast=L,  
               overwrite=False,submitJob=submitJob,PBSsubs=TSPBSsubs)                   

if OS:
    SetupOSrun(geom+'/CGrun',geom+'/OSrun',
               overwrite=True,submitJob=submitJob,PBSsubs=OSPBSsubs)                        


