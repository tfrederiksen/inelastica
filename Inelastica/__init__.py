"""
============================
Inelastica (:mod:`MakeGeom`)
============================

.. module:: MakeGeom

classes
=======

.. autosummary::
   :toctree:

   Geom

============================
Inelastica (:mod:`MiscMath`)
============================

.. module:: MiscMath

classes
=======

.. autosummary::
   :toctree:

   SpectralMatrix

=======================
Inelastica (:mod:`NEB`)
=======================

.. module:: NEB

classes
=======

.. autosummary::
   :toctree:

   savedData
   step
   SigDir
   SavedSigClass

========================
Inelastica (:mod:`NEGF`)
========================

.. module:: NEGF

classes
=======

.. autosummary::
   :toctree:

   ElectrodeSelfEnergy
   GF

===========================
Inelastica (:mod:`Phonons`)
===========================

.. module:: Phonons

classes
=======

.. autosummary::
   :toctree:

   FCrun
   OTSrun
   OSrun
   DynamicalMatrix

====================================
Inelastica (:mod:`SupercellPhonons`)
====================================

.. module:: SupercellPhonons

classes
=======

.. autosummary::
   :toctree:

   Supercell_DynamicalMatrix

============================
Inelastica (:mod:`Symmetry`)
============================

.. module:: Symmetry

classes
=======

.. autosummary::
   :toctree:

   Symmetry

"""

# Import version string and the major, minor, micro as well
from . import info
from .info import git_revision as __git_revision__
from .info import version as __version__
from .info import major as __major__
from .info import minor as __minor__
from .info import micro as __micro__

__all__ = [s for s in dir() if not s.startswith('_')]
__all__ += ['__{}__'.format(s) for s in ['version', 'major', 'minor', 'micro', 'git_revision']]
