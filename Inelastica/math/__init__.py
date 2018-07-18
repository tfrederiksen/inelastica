"""

:mod:`Inelastica.math`
=========================

.. autosummary::
   :toctree:

   SpectralMatrix

.. module:: Inelastica.math

"""

from math import *
from .gausskronrod import *
from .hilbert import *
from .misc import *
from .spectral import *
from .sphericalharmonics import *

__all__ = [s for s in dir() if not s.startswith('_')]
