"""
Example: How to use the SetupRuns functions to handle:
  Geometry optimizations   : CGrun
  Overlap calculations     : OSrun
  Force constants (phonons): FCrun
  TranSIESTA               : TSrun
  Phonons                  : PHrun
  Inelastica               : INrun
"""


from Inelastica.SetupRuns import *                     
from Inelastica.Phonons import *                       
import os                                   
import numpy as N

submitJob=True             # Automatically submitt?
nodes='1:ppn=4'            # Nodes for CGrun, FCrun and TSrun
TSPBSsubs = [['$NODES$',nodes]]
PYPBSsubs = [['$NODES$','1:ppn=1']]
OSPBSsubs = [['$NODES$','1:ppn=1']]


CG = False
FC = False
OS = False
TS = False
PH = True
IN = False

FCFL = [[45,45],[46,46]]    # Split FCrun into these ranges
FCF = N.min(N.array(FCFL))  # (vibrational subspace) 
FCL = N.max(N.array(FCFL))
DevF, DevL = FCF-12, FCL+1  # Device size (e-ph coupling subspace)

geom = 'L4.4'
head,tail = os.path.split(geom)   


if CG:                            
    SetupCGrun(geom+'/CGrun','../L5.4/CGrun',14.044300-3.4+5.4,
               AtomsPerLayer=8,
               IndexShift=0,
               overwrite=False,submitJob=submitJob,                      
               PBSsubs=TSPBSsubs)                    
    
if FC:
    for F,L in FCFL:
        SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
                   FCfirst=F,FClast=L,                     
                   overwrite=False,submitJob=submitJob,          
                   PBSsubs=TSPBSsubs)                        

if OS:
    SetupOSrun(geom+'/CGrun',geom+'/OS',
               overwrite=False,submitJob=submitJob,
               PBSsubs=OSPBSsubs)                        

a1=N.array([
    [0.000000000     , 0.000000000    ,  5.911410000 ],
    [2.090000000     , 1.477850000    ,  7.389270000 ],
    [2.090000000     , 4.433560000    ,  7.389270000 ],
    [4.180000000     , 0.000000000    ,  5.911410000 ],
    [4.180000000     , 2.955710000    ,  5.911410000 ],
    [6.270000000     , 4.433560000    ,  7.389270000 ],
    [6.270000000     , 1.477850000    ,  7.389270000 ],
    [0.000000000     , 2.955710000    ,  5.911410000 ],
    [0.000000000     , 5.911410000    ,  5.911410000 ],
    [2.090000000     , 7.389270000    ,  7.389270000 ],
    [2.090000000     ,10.345000000    ,  7.389270000 ],
    [4.180000000     , 5.911410000    ,  5.911410000 ],
    [4.180000000     , 8.867120000    ,  5.911410000 ],
    [6.270000000     ,10.345000000    ,  7.389270000 ],
    [6.270000000     , 7.389270000    ,  7.389270000 ],
    [0.000000000     , 8.867120000    ,  5.911410000 ]],N.float)
dz =  5.911410000-19.955710000+17
addL = N.zeros((64,3),N.float)
for iz in range(4):
    for jj in range(16):
        addL[iz*16+jj,:]=a1[jj,:]+N.array([0,0,dz*iz])
addR=addL.copy()
            
if TS:
    SetupTSrun(geom+'/CGrun','../TempTSrun',geom+'/TSrun',
               AtomsPerLayer=16,BlockSize=32,
               LayersLeft=0,LayersRight=0,
               AddLeftList=addL,AddRightList=addR,
               overwrite=False,submitJob=submitJob,
               PBSsubs=TSPBSsubs)                        
    
if PH:
    SetupPHrun(geom+'/PHrun','FCrun*',
               onlySdir='../OS',
               DeviceFirst=DevF,DeviceLast=DevL,
               FCfirst=(FCF),FClast=(FCL),
               outlabel=(tail+'Dev+1'),
               overwrite=False,
               submitJob=submitJob,
               PBSsubs=PYPBSsubs)

if IN:
    SetupInelastica('../Inelastica/',geom+'/Inelastica_Dev+1',
                    geom+'/TSrun','./TemplateInelasticaInput.py',
                    tail+'2_Dev+1.py',
                    '../PhononCalcDev+1/%sDev+1.nc'%tail,
                    clustername='e4200',
                    submitJob=submitJob,nodes='1:ppn=1')


