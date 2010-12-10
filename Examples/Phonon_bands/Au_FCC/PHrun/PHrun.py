from Inelastica.Phonons import *


Analyze('../FCrun*',
        onlySdir='../OS',
        DeviceFirst=1,DeviceLast=20,
        FCfirst=130,FClast=130,
        outlabel='Au6x6x6',
        CalcCoupl=False,
        AuxNCfile=False,
        PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
        PrintSOrbitals=True,
        PhBandStruct='FCC')
