"""
==============================
Inelastica (:mod:`Inelastica`)
==============================

.. module:: Inelastica

All classes defined in Inelastica.

All classes
===========

.. autosummary::
   :toctree:

   Geom
   SpectralMatrix
   savedData
   step
   SigDir
   SavedSigClass
   ElectrodeSelfEnergy
   GF
   FCrun
   OTSrun
   OSrun
   DynamicalMatrix
   Supercell_DynamicalMatrix
   Symmetry

"""

# Import version string and the major, minor, micro as well
from . import info
from .info import git_revision as __git_revision__
from .info import version as __version__
from .info import major as __major__
from .info import minor as __minor__
from .info import micro as __micro__

from .BandStruct import *
from .CommonFunctions import *
from .EigenChannels import *
from .iets import *
from .MakeGeom import *
from .MiscMath import *
from .NEB import *
from .NEGF import *
from .Phonons import *
from .pyTBT import *
from .SetupRuns import *
from .STM import *
from .STMFD import *
from .SupercellPhonons import *
from .Symmetry import *
from .ValueCheck import *

__all__ = [s for s in dir() if not s.startswith('_')]
__all__ += ['__{}__'.format(s) for s in ['version', 'major', 'minor', 'micro', 'git_revision']]
