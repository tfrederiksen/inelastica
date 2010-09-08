from Inelastica.Phonons import *

Analyze('..','FCrun*',
        onlySdir='../OS',
        newPHrun='PHrun',
        DeviceFirst=0,DeviceLast=0,
        FCfirst=171,FClast=171,
        output2file=False,outlabel='NaBCC',
        CorrPotentialShift=True,
        TryMirrorSymmetry=False,
        CalcCoupl=False,
        AuxNCfile=False,
        PerBoundCorrFirst=-1,PerBoundCorrLast=-1,
        PrintSOrbitals=True,
        PhBandStruct="BCC")
