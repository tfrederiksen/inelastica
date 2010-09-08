from Inelastica.Phonons import *


Analyze('..','FCrun*',
        onlySdir='../OS',
        newPHrun='PHrun',
        DeviceFirst=1,DeviceLast=20,
        FCfirst=130,FClast=130,
        output2file=False,outlabel='Au6x6x6',
        CorrPotentialShift=True,
        TryMirrorSymmetry=False,
        CalcCoupl=False,
        AuxNCfile=False,
        PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
        PrintSOrbitals=True,
        PhBandStruct='FCC')
