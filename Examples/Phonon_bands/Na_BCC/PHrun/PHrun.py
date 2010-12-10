from Inelastica.Phonons import *

Analyze('../FCrun*',
        onlySdir='../OS',
        DeviceFirst=0,DeviceLast=0,
        FCfirst=171,FClast=171,
        outlabel='NaBCC',
        CalcCoupl=False,
        AuxNCfile=False,
        PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
        PrintSOrbitals=True,
        PhBandStruct="BCC",
        PhBandRadie=7.0)
