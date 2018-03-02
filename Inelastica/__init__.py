"""
==============================
Inelastica (:mod:`Inelastica`)
==============================

.. module:: Inelastica

Inelastica is a `Python`_ package that provides modules/scripts for
electronic structure and elastic/inelastic transport calculations.

"""

__all__ = [s for s in dir() if not s.startswith('_')]

try:
    # Import version string and the major, minor, micro as well
    from . import info
    from .info import git_revision as __git_revision__
    from .info import version as __version__
    from .info import major as __major__
    from .info import minor as __minor__
    from .info import micro as __micro__

    __all__ = ['__{}__'.format(s) for s in ['version', 'major', 'minor', 'micro', 'git_revision']]
except Exception as e:
    print('Could not import info:')
    print(str(e))
