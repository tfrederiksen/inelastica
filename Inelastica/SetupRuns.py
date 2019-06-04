"""

:mod:`Inelastica.SetupRuns`
===========================

Functions to setup (Tran)SIESTA calculation and submit pbs jobs.

* CG: Optimization including scaling of electrode distances
* FC: For phonon and e-ph calculations
* OS: only overlap (S) calculation used to obtain e-ph coupling
* PH: Phonon frequencies and e-ph coupling
* TS: Transiesta adding electrode layers to CG geometry
* IN: (not working at the moment)

The basic idea is to use an existing CG (geometry optimization) setup
to create all the other calculation directories needed.

.. currentmodule:: Inelastica.SetupRuns


"""
from __future__ import print_function

import glob
import time
import sys
import shutil
import os
import netCDF4 as NC4
import numpy as N
import Inelastica.MakeGeom as MG
import Inelastica.physics.constants as PC
from Inelastica.io.siesta import copy_chemical_info


# -----------------------------------------------------------------------------------------------------
#                                          SetupCGrun
# -----------------------------------------------------------------------------------------------------

def SetupCGrun(templateCGrun, newCGrun, NewContactSeparation, AtomsPerLayer,
               IndexShift=0, ListL=None, ListR=None, ZmatRanges=None,
               RotationAngle=None, RotationCenter=None, RotationAxis=None, RotationSubset=None,
               overwrite=False, PBStemplate=None, PBSsubs=None, submitJob=False):
    """
    templateCGrun        : Path+foldername to a template CGrun folder which contain
                              the necessary SIESTA files and a structure which is to
                              be modified for a new electrode separation
                              (this function writes a new STRUCT.fdf file)
    newCGrun             : Path+foldername for the CGrun to be created
    NewContactSeparation : The new electrode separation in Ang (defined as the z-distance
                              over the device between two electrode layers which has not
                              been relaxed)
    AtomsPerLayer        : Number of atoms in one electrode layer in a cross-section wrt.
                              the transport direction
    IndexShift           : Number of atoms which is periodically translated from left to
                              right during the formation of the new STRUCT.fdf from the
                              template *.XV file
    ListL/ListR          : These atom indices (SIESTA numbering) are forced to be a part of
                              the electrodes and hence not affected by stretching
    RotationAngle/       : Rotate whole CG geometry (or only RotationSubset if specified)
    RotationCenter/           an angle RotationAngle (in degrees) around RotationAxis vector through
    RotationAxis/             atom index RotationCenter (SIESTA numbering). These manipulations
    RotationSubset            are implemented BEFORE electrode separation is adjusted.
    overwrite            : (True/False) Specifies if one is allowed to write in existing
                              directories
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...]' will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """

    # Make new directories
    head = os.path.split(newCGrun)[0]
    if not os.path.isdir(head):
        print('\nSetupRuns.SetupCGrun: Creating folder %s' %head)
        os.mkdir(head)
    if not os.path.isdir(newCGrun):
        print('\nSetupRuns.SetupCGrun: Creating folder %s' %newCGrun)
        os.mkdir(newCGrun)
    elif not overwrite:
        print('\nSetupRuns.SetupCGrun: %s already exists. No further actions taken.'\
              %newCGrun)
        return '' # Quit script
    else:
        print('\nSetupRuns.SetupCGrun: %s already exists. OVERWRITING FILES!!!'\
              %newCGrun)
    # Copy template files
    CopyInputFiles(templateCGrun, newCGrun, ['.fdf', '.vps', '.psf', 'pbs', '.TSHS'])
    # Read relaxed geometry
    XVfiles = glob.glob(templateCGrun+'/*.XV*')
    if len(XVfiles) == 1:
        geom = MG.Geom(XVfiles[0])
    elif len(XVfiles) > 1:
        print('More than one XV file was found in folder %s:'%templateCGrun)
        for i, xvfile in enumerate(XVfiles):
            print('   No. %i :'%i, xvfile)
        select = input('   ... select file:')
        geom = MG.Geom(XVfiles[int(select)])
    else:
        print('No XV file was found in folder %s:'%templateCGrun)
        input('   ... Continue reading geometry from RUN.fdf?')
        geom = MG.Geom(templateCGrun+'/RUN.fdf')
    # Rotate via indexshift?
    if IndexShift > 0:
        print('SetupRuns.SetupCGrun: Applying IndexShift =', IndexShift)
        for ii in range(IndexShift):
            geom.xyz[0][2] += geom.pbc[2][2]
            geom.addAtom(geom.xyz[0], geom.snr[0], geom.anr[0])
            geom.rmAtom(0)
    # Rotation?
    if RotationAngle and RotationCenter and RotationAxis:
        print('SetupRuns.SetupCGrun: Rotation applied:')
        print('   ... Rotation angle  =', RotationAngle, ' deg')
        print('   ... Rotation center = atom %i (SIESTA numbering)'%RotationCenter)
        print('   ... Rotation axis   =', RotationAxis)
        if RotationSubset:
            print('   ... RotationSubset  =', RotationSubset, ' (SIESTA numbering)')
            # Change from SIESTA numbering to Python indices
            RotationSubset = [x-1 for x in RotationSubset]
        center = N.array(geom.xyz[RotationCenter-1])
        geom.rotate(RotationAxis, RotationAngle, center, RotateSubset=RotationSubset)
    # Adjust electrode separation
    if NewContactSeparation:
        geom.stretch2NewContactSeparation(NewContactSeparation, AtomsPerLayer, ListL=ListL, ListR=ListR)
    if not ZmatRanges:
        geom.writeFDF(newCGrun+'/STRUCT.fdf')
    else:
        geom.writeFDFZmat(newCGrun+'/STRUCT.fdf', ZmatRanges)
    geom.writeXYZ(newCGrun+'/STRUCT.xyz')
    geom.writeXYZ(newCGrun+'/STRUCT2.xyz', rep=[2, 2, 2])
    # PBS files
    MakePBS(PBStemplate, newCGrun+'/RUN.pbs', PBSsubs, submitJob, rtype='TS')


# -----------------------------------------------------------------------------------------------------
#                                          SetupFCrun
# -----------------------------------------------------------------------------------------------------

def SetupFCrun(CGrun, newFCrun, FCfirst, FClast, displacement=0.02,
               overwrite=False, PBStemplate=None, PBSsubs=None, submitJob=False):
    """
    CGrun                : Path+foldername to a relaxed structure CGrun folder on
                              which to perform a FCrun calculation. This folder must
                              contain an *.XV file
    newFCrun             : Path+foldername for the FCrun to be created
    FCfirst              : Number of the first atom to displace (SIESTA numbering)
    FClast               : Number of the last atom to displace (SIESTA numbering)
    displacement         : Displacement
    overwrite            : (True/False) Specifies if one is allowed to write in existing
                              directories
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...]' will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """
    # Make new directory
    if not os.path.isdir(newFCrun):
        print('\nSetupRuns.SetupFCrun: Creating folder %s' %newFCrun)
        os.mkdir(newFCrun)
    elif not overwrite:
        print('\nSetupRuns.SetupFCrun: %s already exists. No further actions taken.'\
              %newFCrun)
        return '' # Quit script
    else:
        print('\nSetupRuns.SetupFCrun: %s already exists. OVERWRITING FILES!!!'\
              %newFCrun)
    # Copy template files
    CopyInputFiles(CGrun, newFCrun, ['.fdf', '.vps', '.psf', '.DM', '.XV', '.pbs', '.TSDE', '.TSHS'])
    # Read relaxed geometry and overwrite STRUCT files
    XVfiles = glob.glob(CGrun+'/*.XV*')
    if len(XVfiles) == 1:
        geom = MG.Geom(XVfiles[0])
    elif len(XVfiles) > 1:
        print('More than one XV file was found in folder %s:'%CGrun)
        for i, xvfile in enumerate(XVfiles):
            print('   No. %i :'%i, xvfile)
        select = input('   ... select file:')
        geom = MG.Geom(XVfiles[int(select)])
    else:
        print('No XV file was found in folder %s:'%CGrun)
        input('   ... Continue reading geometry from RUN.fdf?')
        geom = MG.Geom(CGrun+'/RUN.fdf')
    geom.writeFDF(newFCrun+'/STRUCT.fdf')
    geom.writeXYZ(newFCrun+'/STRUCT.xyz')
    geom.writeXYZ(newFCrun+'/STRUCT2.xyz', rep=[2, 2, 2])
    if os.path.isfile(CGrun+"/STRUCT.fdf"):
        copy_chemical_info(CGrun+"/STRUCT.fdf", newFCrun+"/STRUCT.fdf")

    # Prepend lines to RUN.fdf
    for elm in glob.glob(newFCrun+'/RUN.fdf'):
        if os.path.isfile(elm):
            print('SetupRuns.SetupFCrun: Modifying %s' %elm)
            f = open(elm, 'r')
            lines = f.readlines()
            f.close()
            f = open(elm, 'w')
            f.write('### Lines appended %s \n' %time.ctime())
            f.write('MD.TypeOfRun  FC\n')
            f.write('MD.FCfirst   %i\n' %FCfirst)
            f.write('MD.FClast    %i\n' %FClast)
            f.write('TS.HS.Save True\n')
            f.write('TS.SaveHS    True\n')
            f.write('MD.FCDispl   %.8f Ang\n'%displacement)
            for line in lines: f.write(line)
            f.close()
    # PBS files
    MakePBS(PBStemplate, newFCrun+'/RUN.pbs', PBSsubs, submitJob, rtype='TS')


# -----------------------------------------------------------------------------------------------------
#                                          SetupOSrun
# -----------------------------------------------------------------------------------------------------

def SetupOSrun(CGrun, newOSrun, displacement=0.02,
               overwrite=False, PBStemplate=None, PBSsubs=None, submitJob=False):
    """
    CGrun                : Path+foldername to a relaxed structure CGrun folder on
                              which to run an onlyS calculation. This folder must
                              contain an *.XV file
    newOSrun             : Path+foldername for the onlyS folder to be created
    displacement         : Displacement used in the FCrun displacements
    overwrite            : (True/False) Specifies if one is allowed to write in existing
                              directories
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...]' will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """
    # Make new directory
    if not os.path.isdir(newOSrun):
        print('\nSetupRuns.SetupOSrun: Creating folder %s' %newOSrun)
        os.mkdir(newOSrun)
    elif not overwrite:
        print('\nSetupRuns.SetupOSrun: %s already exists. No further actions taken.'\
              %newOSrun)
        return '' # Quit script
    else:
        print('\nSetupRuns.SetupOSrun: %s already exists. OVERWRITING FILES!!!'\
              %newOSrun)
    # Copy files from CGrun
    CopyInputFiles(CGrun, newOSrun, ['.fdf', '.vps', '.psf'])
    # Read original RUN.fdf file
    f = open(newOSrun+'/RUN.fdf', 'r')
    lines = f.readlines()
    f.close()
    # Read relaxed geometry
    XVfiles = glob.glob(CGrun+'/*.XV*')
    if len(XVfiles) == 1:
        infile = XVfiles[0]
    elif len(XVfiles) > 1:
        print('More than one XV file was found in folder %s:'%CGrun)
        for i, xvfile in enumerate(XVfiles):
            print('   No. %i :'%i, xvfile)
        select = input('   ... select file:')
        infile = XVfiles[int(select)]
    else:
        print('No XV file was found in folder %s:'%CGrun)
        input('   ... Continue reading geometry from RUN.fdf?')
        infile = CGrun+'/RUN.fdf'
    # Multiply structure
    BuildOSstruct(infile, newOSrun+'/STRUCT_1.fdf', axes=[0], direction=[-1], displacement=displacement)
    BuildOSstruct(infile, newOSrun+'/STRUCT_2.fdf', axes=[0], direction=[1], displacement=displacement)
    BuildOSstruct(infile, newOSrun+'/STRUCT_3.fdf', axes=[1], direction=[-1], displacement=displacement)
    BuildOSstruct(infile, newOSrun+'/STRUCT_4.fdf', axes=[1], direction=[1], displacement=displacement)
    BuildOSstruct(infile, newOSrun+'/STRUCT_5.fdf', axes=[2], direction=[-1], displacement=displacement)
    BuildOSstruct(infile, newOSrun+'/STRUCT_6.fdf', axes=[2], direction=[1], displacement=displacement)
    structfiles = ['STRUCT_1.fdf', 'STRUCT_2.fdf', 'STRUCT_3.fdf',
                   'STRUCT_4.fdf', 'STRUCT_5.fdf', 'STRUCT_6.fdf']
    if os.path.isfile(CGrun+"/STRUCT.fdf"):
        for struct in structfiles:
            copy_chemical_info(CGrun+"/STRUCT.fdf", newOSrun+"/"+struct)
    inputfiles = ['RUN_1.fdf', 'RUN_2.fdf', 'RUN_3.fdf', 'RUN_4.fdf', 'RUN_5.fdf', 'RUN_6.fdf']
    # Write input files
    for i, inputfile in enumerate(inputfiles):
        print('SetupRuns.SetupOSrun: Writing %s' %(newOSrun+'/'+inputfile))
        f = open((newOSrun+'/'+inputfile), 'w')
        f.write('### Lines written %s \n' %time.ctime())
        #f.write('MD.NumCGSteps 0\n')
        #f.write('DM.NumberPulay 0 \n')
        #f.write('MPN.onlyS .true.\n\n')
        f.write('TS.onlyS .true.\n')
        f.write('SystemName STRUCT_%i\n'%(i+1))
        f.write('SystemLabel STRUCT_%i\n\n'%(i+1))
        f.write('%include '+'./STRUCT_%i.fdf\n\n'%(i+1))
        f.write('### Lines from RUN.fdf \n')
        for line in lines:
            if line.strip().endswith('STRUCT.fdf'):
                f.write('%include '+structfiles[i]+'\n')
            else: f.write(line)
        f.close()
    # PBS files
    MakePBS(PBStemplate, newOSrun+'/RUN.pbs', PBSsubs, submitJob, rtype='OS')


# -----------------------------------------------------------------------------------------------------
#                                          SetupTSrun
# -----------------------------------------------------------------------------------------------------

def SetupTSrun(CGrun, templateTSrun, newTSrun,
               AtomsPerLayer, BlockSize, LayersLeft, LayersRight,
               DeviceLayerInclLeft=0, DeviceLayerInclRight=0, IndexShift=0,
               AddLeftList=[], AddRightList=[], ListL=None, ListR=None,
               NewContactSeparation=None,
               RotationAngle=None, RotationCenter=None, RotationAxis=None, RotationSubset=None,
               overwrite=False, PBStemplate=None, PBSsubs=None, submitJob=False):
    """

    CGrun                : Path+foldername to a relaxed structure CGrun folder on
                              which to run a TranSIESTA calculation. This folder must
                              contain an *.XV file
    templateTSrun        : Path+foldername to an existing TSrun folder which contains
                              all relevant electrode calculations and *.fdf files
                              except for STRUCT.fdf
    newTSrun             : Path+foldername for the TSrun folder to be created
    AtomsPerLayer        : Number of atoms in the electrodes in a cross section along
                              the transport direction
    BlockSize            : Number of atoms in the electrode block from CGrun/*.XV file
    LayersLeft/Right     : Number of blocks to be pasted to the new enlarged STRUCT.fdf.
    DeviceLayerInclLeft  : As default TBTRANS takes all relaxed atoms to be the device
                              If DeviceLayerInclLeft is larger than 0 then TBTRANS
                              includes the specified number of electode atomic planes into
                              the device region.
    DeviceLayerInclRight : See DeviceLayerInclLeft
    IndexShift           : Number of atoms which is periodically translated from left to
                              right contact BEFORE eventual electrode layers are pasted.
    AddLeftList          : Add list of atoms to the left side
    AddRightList         : Add list of atoms to the right side
    ListL/ListR          : These atom indices (SIESTA numbering) are forced to be a part of
                              the electrodes and hence not affected by stretching
    NewContactSeparation : (Optional) value to set a different electrode separation
                              than in the CG geometry
    RotationAngle/       : Rotate whole CG geometry (or only RotationSubset if specified)
    RotationCenter/           an angle RotationAngle (in degrees) around RotationAxis vector through
    RotationAxis/             atom index RotationCenter (SIESTA numbering). These manipulations
    RotationSubset            are implemented BEFORE electrode separation is adjusted.
    overwrite            : (True/False) Specifies if one is allowed to write in existing
                              directories
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...]' will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """
    # Make new directory
    if not os.path.isdir(newTSrun):
        print('\nSetupRuns.SetupTSrun: Creating folder %s' %newTSrun)
        os.mkdir(newTSrun)
    elif not overwrite:
        print('\nSetupRuns.SetupTSrun: %s already exists. No further actions taken.'\
              %newTSrun)
        return '' # Quit script
    else:
        print('\nSetupRuns.SetupTSrun: %s already exists. OVERWRITING FILES!!!'\
              %newTSrun)
    # Copy all sub-directories
    for elm in glob.glob(templateTSrun+'/*'):
        tail = os.path.split(elm)[1]
        if os.path.isdir(elm):
            CopyTree(elm, newTSrun+'/'+tail, overwrite=overwrite)
    # Copy template files
    CopyInputFiles(templateTSrun, newTSrun, ['.fdf', '.vps', '.psf', '.pbs'])
    # Read relaxed geometry
    XVfiles = glob.glob(CGrun+'/*.XV*')
    if len(XVfiles) == 1:
        geom = MG.Geom(XVfiles[0])
    elif len(XVfiles) > 1:
        print('More than one XV file was found in folder %s:'%CGrun)
        for i, xvfile in enumerate(XVfiles):
            print('   No. %i :'%i, xvfile)
        select = input('   ... select file:')
        geom = MG.Geom(XVfiles[int(select)])
    else:
        print('No XV file was found in folder %s:'%CGrun)
        input('   ... Continue reading geometry from RUN.fdf?')
        geom = MG.Geom(CGrun+'/RUN.fdf')
    # Rotate via indexshift?
    if IndexShift > 0:
        print('SetupRuns.SetupTSrun: Applying IndexShift =', IndexShift)
        for ii in range(IndexShift):
            geom.xyz[0][2] += geom.pbc[2][2]
            geom.addAtom(geom.xyz[0], geom.snr[0], geom.anr[0])
            geom.rmAtom(0)
    # Overwrite STRUCT files
    if RotationAngle and RotationCenter and RotationAxis:
        print('SetupRuns.SetupTSrun: Rotation applied:')
        print('   ... Rotation angle  =', RotationAngle, ' deg')
        print('   ... Rotation center = atom %i (SIESTA numbering)'%RotationCenter)
        print('   ... Rotation axis   =', RotationAxis)
        if RotationSubset:
            print('   ... RotationSubset  =', RotationSubset, ' (SIESTA numbering)')
            # Change from SIESTA numbering to Python indices
            RotationSubset = [x-1 for x in RotationSubset]
        center = N.array(geom.xyz[RotationCenter-1])
        geom.rotate(RotationAxis, RotationAngle, center, RotateSubset=RotationSubset)
    # Paste electrode layers via BlockSize specifications
    geom.PasteElectrodeLayers(BlockSize, AtomsPerLayer, LayersLeft, LayersRight)
    # Change contact separation?
    if NewContactSeparation:
        geom.stretch2NewContactSeparation(NewContactSeparation, AtomsPerLayer, ListL=ListL, ListR=ListR)
    # Add electrode atoms to the left
    if len(AddLeftList) > 0:
        dz = AddLeftList[AtomsPerLayer, 2]-AddLeftList[0, 2]
        tmp = N.array(geom.xyz)
        minz = min(tmp[:, 2])
        maxz = max(AddLeftList[:, 2])
        for ii in reversed(list(range(len(AddLeftList)))):
            tmp = list(AddLeftList[ii, :]+(-maxz+minz-dz)*N.array([0.0, 0.0, 1.0], N.float))
            geom.prependAtom(tmp,
                         geom.snr[0], geom.anr[0])
        geom.pbc[2][2] += len(AddLeftList)/AtomsPerLayer*dz
    # Add electrode atoms to the right
    if len(AddRightList) > 0:
        dz = AddRightList[AtomsPerLayer, 2]-AddRightList[0, 2]
        tmp = N.array(geom.xyz)
        maxz = tmp[0, 2]-dz+geom.pbc[2][2]
        minz = min(AddRightList[:, 2])
        for ii in range(len(AddRightList)):
            geom.addAtom(list(AddRightList[ii, :]+
                              (maxz-minz+dz)*N.array([0, 0, 1], N.float)),
                         geom.snr[0], geom.anr[0])
        geom.pbc[2][2] += len(AddRightList)/AtomsPerLayer*dz
    # Write structure to files
    geom.writeFDF(newTSrun+'/STRUCT.fdf')
    geom.writeXYZ(newTSrun+'/STRUCT.xyz')
    geom.writeXYZ(newTSrun+'/STRUCT2.xyz', rep=[2, 2, 2])
    # Find device atoms
    geom.findContactsAndDevice(AtomsPerLayer)
    PDOSfirst = min(geom.deviceList)-DeviceLayerInclLeft*AtomsPerLayer
    PDOSlast = max(geom.deviceList)+DeviceLayerInclRight*AtomsPerLayer
    # Prepend lines to TBTRANS.fdf
    for elm in glob.glob(newTSrun+'/TBTRANS.fdf'):
        if os.path.isfile(elm):
            print('SetupRuns.SetupTSrun: Prepending PDOS = [%i,%i] to %s' \
                  %(PDOSfirst, PDOSlast, elm))
            f = open(elm, 'r')
            lines = f.readlines()
            f.close()
            f = open(elm, 'w')
            f.write('### Lines appended %s \n' %time.ctime())
            f.write('TS.TBT.PDOSFrom   %i\n' %PDOSfirst)
            f.write('TS.TBT.PDOSTo     %i\n\n' %PDOSlast)
            for line in lines: f.write(line)
            f.close()
    # PBS files
    MakePBS(PBStemplate, newTSrun+'/RUN.pbs', PBSsubs, submitJob, rtype='TS')


# -----------------------------------------------------------------------------------------------------
#                                          RunTBT
# -----------------------------------------------------------------------------------------------------

def RunTBT(TSrun, Emin, Emax, NPoints, NumKxy_A1=1, NumKxy_A2=1,
           AtomsPerLayer=0, DeviceLayerInclLeft=0, DeviceLayerInclRight=0,
           PBStemplate=None, PBSsubs=None, submitJob=False):
    """
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...] will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """
    # Remove old Green's functions
    for elm in glob.glob(TSrun+'/*GF'):
        print('SetupRuns.RunTBT: Removing *.GF')
        if os.path.isfile(elm):
            print('   Deleting %s'%elm)
            os.remove(elm)
    # Determine Device PDOS?
    if AtomsPerLayer > 0:
        for elm in glob.glob(TSrun+'/*XV'):
            if os.path.isfile(elm):
                geom = MG.Geom(elm)
                geom.findContactsAndDevice(AtomsPerLayer)
                PDOSfirst = min(geom.deviceList)-DeviceLayerInclLeft*AtomsPerLayer
                PDOSlast = max(geom.deviceList)+DeviceLayerInclRight*AtomsPerLayer
    else:
        PDOSfirst, PDOSlast = 0, 0
    # Prepend lines to TBTRANS.fdf
    tbtfile = TSrun+'/TBTRANS.fdf'
    f = open(tbtfile, 'r')
    lines = f.readlines()
    f.close()
    f = open(tbtfile, 'w')
    f.write('### Lines appended %s \n' %time.ctime())
    if PDOSfirst > 0 and PDOSlast > 0:
        print('SetupRuns.RunTBT: Prepending PDOS = [%i,%i] to %s' \
            %(PDOSfirst, PDOSlast, tbtfile))
        f.write('TS.TBT.PDOSFrom   %i\n' %PDOSfirst)
        f.write('TS.TBT.PDOSTo     %i\n' %PDOSlast)
    print('SetupRuns.RunTBT: Writing Emin=%f, Emax=%f, NPoints=%i to %s' \
        %(Emin, Emax, NPoints, tbtfile))
    f.write('TS.TBT.Emin       %f  eV\n' %Emin)
    f.write('TS.TBT.Emax       %f  eV\n' %Emax)
    f.write('TS.TBT.NPoints    %i\n' %NPoints)
    f.write('TS.NumKxy_A1      %i\n' %NumKxy_A1)
    f.write('TS.NumKxy_A2      %i\n' %NumKxy_A2)
    f.write('pyTBT.K_A1        %i\n'%NumKxy_A1)
    f.write('pyTBT.K_A2        %i\n\n'%NumKxy_A2)
    for line in lines: f.write(line)
    f.close()
    # PBS files
    MakePBS(PBStemplate, newTSrun+'/RUN.pyTBT.pbs', PBSsubs, submitJob, rtype='PY')


# -----------------------------------------------------------------------------------------------------
#                                          SetupPHrun
# -----------------------------------------------------------------------------------------------------

def SetupPHrun(newPHrun, wildcard, onlySdir='../OSrun',
               DeviceFirst=1, DeviceLast=1e3, FCfirst=1, FClast=1e3,
               CalcCoupl=True, outlabel='Out',
               PerBoundCorrFirst=-1, PerBoundCorrLast=-1,
               PrintSOrbitals=True, AuxNCfile=False,
               overwrite=False, PBStemplate=None, PBSsubs=None, submitJob=False):
    """
    wildcard             : Path+wildcard used in the search for FCrun folders
    onlySdir             : Path+foldername to the relevant onlyS directory
    DeviceFirst/Last     : SIESTA atom numbers for the first/last atom to be included
                              in the device region (subspace of the Hamiltonian etc.)
    FCfirst/last         : SIESTA atom numbers for the first/last atom to be allowed
                              to move in the phonon calculation
    PerBoundCorrFirst/   : Optional argument to be passed to Phonons.py script
    PerBoundCorrLast
    PrintSOrbitals       : Optional argument to be passed to Phonons.py script
    AuxNCfile            : Optional argument to be passed to Phonons.py script
    overwrite            : (True/False) Specifies if one is allowed to write in existing
                              directories
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...]' will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """

    ### Make directory for output files etc.
    if not os.path.isdir(newPHrun):
        print('SetupRuns.SetupPHrun: Creating', newPHrun)
        os.mkdir(newPHrun)
    else:
        if not overwrite:
            print("Error: directory already exist ", newPHrun)
            sys.exit(1)

    # find device?
    print('SetupRuns.SetupPHrun: Writing', newPHrun+'/PHrun.py')
    phfile = open(newPHrun+'/PHrun.py', 'w')
    phfile.write('from Inelastica.Phonons import *\n\n')
    phfile.write('\nAnalyze(FCwildcard=\'%s\',onlySdir=\'%s\',\n' %(wildcard, onlySdir))
    phfile.write('        DeviceFirst=%s,DeviceLast=%s,\n' %(DeviceFirst, DeviceLast))
    phfile.write('        FCfirst=%s,FClast=%s,\n' %(FCfirst, FClast))
    phfile.write('        outlabel=\'%s\',\n'%outlabel)
    phfile.write('        CalcCoupl=%s,\n' %CalcCoupl)
    if AuxNCfile:
        phfile.write('        AuxNCfile=\'%s\',\n' %AuxNCfile)
    else:
        phfile.write('        AuxNCfile=False,\n')
    phfile.write('        PerBoundCorrFirst=%i,PerBoundCorrLast=%i,\n'
               %(PerBoundCorrFirst, PerBoundCorrLast))
    phfile.write('        PrintSOrbitals=%s)' %PrintSOrbitals)
    phfile.close()

    # write WritePythonPBS(...)
    # PBS files
    if PBSsubs == None:
        PBSsubs = []
    PBSsubs = [['$PYTHONSCRIPT$', 'PHrun.py']] + PBSsubs
    MakePBS(PBStemplate, newPHrun+'/RUN.pbs', PBSsubs, submitJob, rtype='PY')


# -----------------------------------------------------------------------------------------------------
#                                          SetupInelastica
# -----------------------------------------------------------------------------------------------------

def SetupInelastica(templateInelastica, newInelastica, TSrun,
                    templateInputFile, newInputFilename, PhononNCfile,
                    kPoints=1,
                    overwrite=False, PBStemplate=None, PBSsubs=None, submitJob=False):
    """
    templateInelastica   : Path+foldername to the latest version of Inelastica
    newInelastica        : Path+foldername for the Inelastica folder to be created
    TSrun                : Path+foldername to the TSrun where TBTRANS outputs 'HSSigmaLR.nc'
    templateInputFile    : Path+filename to an Inelastica input file which holds all
                              input variables except for the references to
                              RunThis.ncfile_e and RunThis.ncfile_ph (which are appended
                              automatically)
    newInputFilename     : Filename for the Inelastica input file to be created in folder
                              newInelastica with the above references appended
    PhononNCfile         : Path+filename to the NetCDF-file with phonon modes and couplings
                              (RunThis.ncfile_ph)
    kPoints              : Number of k-points written by TBTrans
    overwrite            : (True/False) Specifies if one is allowed to write in existing
                              directories
    PBStemplate          : Path+foldername to a template RUN.pbs file for PBS queueing
    PBSsubs              : A list of string substitutions to be applied to the template
                              PBS script in order to generate a new PBS script
                              (e.g., PBSsubs=[['JOBNAME','newjobname'],...]' will replace
                              any JOBNAME string with newjobname)
    submitJob            : (True/False) Submit to batch queue via qsub command?
    """
    print('SetupRuns.SetupInelastica: Creating', newInelastica)
    CopyTree(templateInelastica, newInelastica, overwrite=False)

    name, ext = os.path.splitext(newInputFilename)
    cmmd = []
    for i in range(kPoints):
        if kPoints == 1:
            inputfile = newInputFilename
            TBT_NCfile = '/HSSigmaLR.nc'
            Mod_NCfile = '/HSSigmaLR_mod.nc'
        else:
            inputfile = name+'_%i'%(i+1)+ext
            TBT_NCfile = '/HSSigmaLR_%.3i.nc'%(i+1)
            Mod_NCfile = '/HSSigmaLR_mod_%i.nc'%(i+1)
        # Convert TBT nc-files to Inelastica folder
        if not os.path.isfile(newInelastica+Mod_NCfile) or overwrite:
            UnitConvertTBOutput(TSrun+TBT_NCfile, newInelastica+Mod_NCfile)
        print('SetupRuns.SetupInelastica: Creating', newInelastica+'/'+inputfile)
        shutil.copy(templateInputFile, newInelastica+'/'+inputfile)
        infile = open(newInelastica+'/'+inputfile, 'a')
        infile.write('\n# AUTOMATICALLY APPENDED LINES: \n')
        infile.write('RunThis.ncfile_e    = \'.%s\'\n'%Mod_NCfile)
        infile.write('RunThis.ncfile_ph    = \'%s\'\n' %PhononNCfile)
        infile.close()
        cmmd.append('inelastica.py '+inputfile+' -R')
    # PBS files
    MakePBS(PBStemplate, newInelastica+'/RUN.pbs', PBSsubs, submitJob, rtype='PY')


# -----------------------------------------------------------------------------------------------------
#                                          Auxiliary functions
# -----------------------------------------------------------------------------------------------------

def RemoveOutputFiles(folder, suffixes=[], DryRun=True, FullCleanup=False):
    print('SetupRuns.RemoveOutputFiles: Cleaning up %s'%folder)
    # These suffixes can always be deleted
    suffixes += ['.ion', '.ion.xml', '.POT.CONF', 'INPUT_TMP', 'out.fdf',
                 '.alloc', '.Au_pbr', 'WALLTIME', 'fort.66',
                 'FDF_STDIN', '.GF', '.MD', '.MDE', '.VH', '.VL', '.VT']
    # ['.HS','.ion.nc','.TSHS','.TSHS','.TSDE','.ANI','.DM',
    #  '.EIG','.TEIG','.IEIG','.FC','.RHO','.DRHO','.nc']
    for elm in glob.glob(folder+'/*'):
        for s in suffixes:
            if elm.endswith(s) and os.path.isfile(elm):
                print('   Deleting %s'%elm)
                if not DryRun: os.remove(elm)

    if FullCleanup:
        # Remove all files except the input files
        keepSuffixes = ['.fdf', '.pbs', '.log', '.xyz', '.XV', '.vps', '.psf', '.py', '.dat', '.dat.nc', '.freq']
        for elm in glob.glob(folder+'/*'):
            for s in keepSuffixes:
                if elm.endswith(s) and os.path.isfile(elm):
                    print('   Keeping %s'%elm)
                else:
                    print('   Deleting %s'%elm)
                    if not DryRun: os.remove(elm)
                    break


def ZipFolder(folder):
    print('SetupRuns: Zipping ', folder)
    head, tail = os.path.split(folder)
    cwd = os.getcwd()
    os.chdir(head)
    os.system('tar -cf tmp.tar %s'%tail)
    os.system('gzip -c tmp.tar > %s.tar.gz'%tail)
    os.system('rm tmp.tar')
    os.system('rm -r %s'%tail)
    os.chdir(cwd)


def UnzipFolder(zfile):
    print('SetupRuns: Unzipping', zfile)
    os.system('gunzip %s'%zfile)
    head, tail = os.path.split(zfile)
    cwd = os.getcwd()
    os.chdir(head)
    os.system('tar -xf %s'%tail[:-3])
    os.system('rm %s'%tail[:-3])
    os.chdir(cwd)


def CopyInputFiles(infolder, outfolder, suffixes):
    print('SetupRuns.CopyInputFiles:')
    for elm in glob.glob(infolder+'/*'):
        for s in suffixes:
            if elm.endswith(s) and os.path.isfile(elm):
                print('   %s  -->  %s'%(elm, outfolder))
                shutil.copy(elm, outfolder)


def BuildOSstruct(infile, outfile, axes=[0, 1, 2], direction=[-1, 1], displacement=0.04*PC.Bohr2Ang):
    print('BuildOSstruct: displacement = %.6f Ang'%displacement)
    geom = MG.Geom(infile)
    for i in axes:
        for displdir in direction:
            new = MG.Geom(infile)
            displ = [0., 0., 0.]
            displ[i] = displdir*displacement
            new.move(displ)
            geom.addGeom(new)
    geom.writeFDF(outfile)


def CopyTree(templateFolder, newFolder, overwrite=False):
    # Make new directory
    if os.path.isdir(newFolder):
        if overwrite:
            print('SetupRuns.CopyTree: %s already exists. REMOVING EXISTING FOLDER!!!' %newFolder)
            shutil.rmtree(newFolder, ignore_errors=True)
        else:
            print('SetupRuns.CopyTree: %s already exists. No further actions taken.' %newFolder)
            return
    # Copy tree
    print('SetupRuns.CopyTree: Copying')
    print('   %s  -->  %s'%(templateFolder, newFolder))
    shutil.copytree(templateFolder, newFolder)


def UnitConvertTBOutput(TBncfile, newTBncfile):
    # This function is copied from "netcdf_transiestaInterpol.py"
    print('SetupRuns.UnitConvertTBOutput: %s  -->  %s' %(TBncfile, newTBncfile))

    # Read TBncfile (infile)
    infile = NC4.Dataset(TBncfile, 'r')
    En = PC.Rydberg2eV*N.array(infile.variables['En'][:, 0]) \
        + 1j*PC.Rydberg2eV*N.array(infile.variables['En'][:, 1])
    H = PC.Rydberg2eV*N.array(infile.variables['H'][:])
    S = N.array(infile.variables['S'][:])
    try:
        ImH = PC.Rydberg2eV*N.array(infile.variables['ImH'][:])
        ImS = N.array(infile.variables['ImS'][:])
        ImHSexist = True
    except:
        ImHSexist = False
        print('UnitConvertTBOutput: No imaginary parts of H and S found.')
    x = len(En)
    dim = len(H)

    # Write variables to newTBncfile (outfile)
    outfile = NC4.Dataset(newTBncfile, 'w')
    outfile.title = "TBtrans unit-converted output"
    outfile.version = 1
    outfile.createDimension('x', x) # Energy grid
    outfile.createDimension('flt', 1) # Float value
    outfile.createDimension('dim', dim) # Problem dimension
    H2 = outfile.createVariable('H', N.float, ('dim', 'dim'))
    H2[:] = H
    H2.units = 'eV'
    S2 = outfile.createVariable('S', N.float, ('dim', 'dim'))
    S2[:] = S
    S2.units = 'none'
    if ImHSexist:
        ImH2 = outfile.createVariable('ImH', N.float, ('dim', 'dim'))
        ImH2[:] = ImH
        ImH2.units = 'eV'
        ImS2 = outfile.createVariable('ImS', N.float, ('dim', 'dim'))
        ImS2[:] = ImS
        ImS2.units = 'none'
    tmp = outfile.createVariable('ReEn', N.float, ('x',))
    tmp[:] = En.real
    tmp.units = 'eV'
    tmp2 = outfile.createVariable('ImEn', N.float, ('x',))
    tmp2[:] = En.imag
    tmp2.units = 'eV'
    ReSigmaL = PC.Rydberg2eV*N.array(infile.variables['ReSigmaL'][:])
    ReSigmaL2 = outfile.createVariable('ReSigmaL', N.float, ('x', 'dim', 'dim'))
    ReSigmaL2[:, :, :] = ReSigmaL[:, :, :]
    ReSigmaL2.units = 'eV'
    ReSigmaL = []
    ImSigmaL = PC.Rydberg2eV*N.array(infile.variables['ImSigmaL'][:])
    ImSigmaL2 = outfile.createVariable('ImSigmaL', N.float, ('x', 'dim', 'dim'))
    ImSigmaL2[:, :, :] = ImSigmaL[:, :, :]
    ImSigmaL2.units = 'eV'
    ImSigmaL = []
    ReSigmaR = PC.Rydberg2eV*N.array(infile.variables['ReSigmaR'][:])
    ReSigmaR2 = outfile.createVariable('ReSigmaR', N.float, ('x', 'dim', 'dim'))
    ReSigmaR2[:, :, :] = ReSigmaR[:, :, :]
    ReSigmaR2.units = 'eV'
    ReSigmaR = []
    ImSigmaR = PC.Rydberg2eV*N.array(infile.variables['ImSigmaR'][:])
    ImSigmaR2 = outfile.createVariable('ImSigmaR', N.float, ('x', 'dim', 'dim'))
    ImSigmaR2[:, :, :] = ImSigmaR[:, :, :]
    ImSigmaR2.units = 'eV'
    ImSigmaR = []
    # Close files
    infile.close()
    outfile.close()


def CheckIfFinished(outfile):
    print()
    try:
        f = open(outfile, 'r')
        for line in f.readlines(): pass
        if line[:6] == '>> End' or line[:6] == '======':
            print('CheckIfFinished: %s is done.'%outfile)
            return True
        else:
            print('CheckIfFinished: %s is NOT done.'%outfile)
            return False
    except:
        print('CheckIfFinished: %s does NOT exist!'%outfile)
        return False


def FindElectrodeSep(directory, AtomsPerLayer):
    for xvfile in glob.glob(directory+'/*.XV*'):
        print(xvfile)
        g = MG.Geom(xvfile)
        g.findContactsAndDevice(AtomsPerLayer)
        DeviceFirst, DeviceLast = g.deviceList[0], g.deviceList[-1]
    return g.ContactSeparation, DeviceFirst, DeviceLast


def MakePBS(PBStemplate, PBSout, PBSsubs, submitJob, rtype='TS'):
    if PBStemplate == None:
        rtypes = {'TS': 'RUN.TS.pbs', 'OS': 'RUN.OS.pbs', 'PY': 'RUN.py.pbs'}
        PBStemplate = rtypes[rtype]
        if os.path.exists(os.path.expanduser('~/.Inelastica/'+PBStemplate)):
            PBStemplate = os.path.expanduser('~/.Inelastica/'+PBStemplate)
        else:
            InelasticaDir = os.path.split(__file__)[0]
            PBStemplate = os.path.abspath(InelasticaDir+'/PBS/'+PBStemplate)

    if os.path.exists(PBStemplate):
        WritePBS(PBStemplate, PBSout, PBSsubs)
        if submitJob:
            print(PBStemplate)
            workingFolder, PBSfile = os.path.split(os.path.abspath(PBSout))
            SubmitPBS(workingFolder, PBSfile)
    else:
        print("WARNING: Could not find PBS template file", PBStemplate)


def WritePBS(PBStemplate, PBSout, PBSsubs):
    print('SetupRuns.WritePBS: Reading', PBStemplate)
    print('SetupRuns.WritePBS: Writing', PBSout)

    # Make default job name
    fullPath = os.path.split(os.path.abspath(PBSout))[0]
    last2dir = fullPath.split('/')[-2:]
    try: # Check for numbers at start ... not liked by PBS
        int(last2dir[0][0])+1
        last2dir[0] = 'a'+last2dir[0]
    except Exception as e:
        print(e)
    if not PBSsubs: PBSsubs = []
    newPBSsub = PBSsubs+[['$DEFJOBNAME$', last2dir[0]+'-'+last2dir[1]]]
    infile = open(PBStemplate)
    outfile = open(PBSout, 'w')
    for line in infile:
        for sub in newPBSsub:
            line = line.replace(str(sub[0]), str(sub[1]))
        outfile.write(line)
    infile.close()
    outfile.close()


def SubmitPBS(workingfolder, pbsfile):
    cwd = os.getcwd()
    os.chdir(workingfolder)
    os.system('qsub '+pbsfile)
    os.chdir(cwd)
    print('   ...and submitted!!!')
