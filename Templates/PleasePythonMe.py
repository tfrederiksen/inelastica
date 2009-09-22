from SetupRuns import *                     
from Phonons import *                       
import os                                   
import numpy as N

clustername='bionano'
submitjob=True

if False:
    CG = False
    FC = True
    OS = True
    TS = True
    Ph = False
    In = False
else:
    CG = False
    FC = False
    OS = False
    TS = True
    Ph = False
    In = False

FCFL = [[28,28],[29,29]]       # Split FCrun into these ranges
FCF = N.min(N.array(FCFL))
FCL = N.max(N.array(FCFL))
DevF, DevL = 28-9, 29+9        # Device size

for nnn in [8]: 
    nodes='1:ppn=%i'%nnn       # CGrun, FCrun and TSrun, others serial
    geom = '../%i/'%nnn
    head,tail = os.path.split(geom)   
                                      
                                      
    if CG:                            
        SetupCGrun(geom+'../8/CGrun','../L0.5/CGrun',13.817176-0.5,
                   AtomsPerLayer=9,clustername=clustername,           
                   IndexShift=0,interpolCGrun=None,                     
                   overwrite=False,submitJob=submitjob,                      
                   nodes=nodes,additionalPBS='')                    

    if FC:
        for F,L in FCFL:
            SetupFCrun(geom+'/CGrun',geom+'/FCrun_%i-%i'%(F,L),
                       FCfirst=F,FClast=L,                     
                       clustername=clustername,              
                       overwrite=False,submitJob=submitjob,          
                       nodes=nodes)                        

    if OS:
        SetupOnlyS(geom+'/CGrun',geom+'/onlyS',
                   clustername=clustername,
                   overwrite=False,submitJob=submitjob,          
                   nodes='1:ppn=1')

    if TS:
        AL= N.array([[    0.000000000  ,    0.000000000  ,   -6.408588004],
            [2.616295010  ,    0.000000000  ,   -6.408588004],
            [5.232590019  ,    0.000000000  ,   -6.408588004],
            [1.308147769  ,    2.265778219  ,   -6.408588004],
            [3.924442779  ,    2.265778219  ,   -6.408588004],
            [6.540737789  ,    2.265778219  ,   -6.408588004],
            [2.616295010  ,    4.531555909  ,   -6.408588004],
            [5.232590019  ,    4.531555909  ,   -6.408588004],
            [7.848885029  ,    4.531555909  ,   -6.408588004],
            [0.000000000  ,    1.510518813  ,   -4.272392003],
            [2.616295010  ,    1.510518813  ,   -4.272392003],
            [5.232590019  ,    1.510518813  ,   -4.272392003],
            [1.308147769  ,    3.776296503  ,   -4.272392003],
            [3.924442779  ,    3.776296503  ,   -4.272392003],
            [6.540737789  ,    3.776296503  ,   -4.272392003],
            [2.616295010  ,    6.042074722  ,   -4.272392003],
            [5.232590019  ,    6.042074722  ,   -4.272392003],
            [7.848885029  ,    6.042074722  ,   -4.272392003],
            [0.000000000  ,    3.021037096  ,   -2.136196001],
            [2.616295010  ,    3.021037096  ,   -2.136196001],
            [5.232590019  ,    3.021037096  ,   -2.136196001],
            [1.308147769  ,    5.286815316  ,   -2.136196001],
            [3.924442779  ,    5.286815316  ,   -2.136196001],
            [6.540737789  ,    5.286815316  ,   -2.136196001],
            [2.616295010  ,    7.552593535  ,   -2.136196001],
            [5.232590019  ,    7.552593535  ,   -2.136196001],
            [7.848885029  ,    7.552593535  ,   -2.136196001]])
        AR=AL.copy()
        SetupTSrun(geom+'/CGrun','../TSrun_temp',geom+'/TSrun',
                   AtomsPerLayer=9,BlockSize=27,
                   LayersLeft=0,LayersRight=0,
                   AddLeftList=AL,
                   AddRightList=AR,
                   clustername=clustername,
                   overwrite=False,submitJob=submitjob,
                   nodes=nodes)

    if Ph:
        SetupPhononCalc(geom+'/PhononCalcDev+1','FCrun*',
                        clustername=clustername,
                        DeviceFirst=DevF,DeviceLast=DevL,
                        FCfirst=(FCF),FClast=(FCL),
                        outlabel=(tail+'Dev+1'),
                        append2existingJobs=False,overwrite=False,
                        submitJob=submitjob,nodes='1:ppn=1')

    if In:
        SetupInelastica('../Inelastica/',geom+'/Inelastica_Dev+1',
                        geom+'/TSrun','./TemplateInelasticaInput.py',
                        tail+'2_Dev+1.py',
                        '../PhononCalcDev+1/%sDev+1.nc'%tail,
                        clustername='e4200',
                        submitJob=submitjob,nodes='1:ppn=1')


