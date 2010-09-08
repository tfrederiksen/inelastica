from Inelastica.Phonons import *


Analyze('..','FCrun_*',
        onlySdir='../OS',
        newPHrun='PHrun',
        DeviceFirst=0,DeviceLast=0,
        FCfirst=259,FClast=260,
        output2file=False,outlabel='NaCl',
        CorrPotentialShift=True,
        TryMirrorSymmetry=False,
        CalcCoupl=False,
        AuxNCfile=False,
        PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
        PrintSOrbitals=True,
        PhBandStruct='FCC')
