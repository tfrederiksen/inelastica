from Inelastica.Phonons import *


Analyze('../FCrun_*',
        onlySdir='../OS',
        DeviceFirst=0,DeviceLast=0,
        FCfirst=259,FClast=260,
        outlabel='NaCl',
        CalcCoupl=False,
        AuxNCfile=False,
        PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
        PrintSOrbitals=True,
        PhBandStruct='FCC')
