from __future__ import print_function

import numpy as N
import scipy.special as SS


def sphericalHarmonics(l, m, costh, sinfi, cosfi):
    r"""
    Spherical harmonics :math:`Y_\ell^m(\theta,\phi)` derived from ``scipy.special``.

    Note: This function has been checked up to d-orbitals with the Siesta overlap matrix.

    Parameters
    ----------
    l : int
    m : int
    costh : float
    sinfi : float
    cosfi : float

    Returns
    -------
    Ylm : complex
    """
    # New faster Spherical Harmonics. Checked up to d-orbitals with the Siesta overlap matrix.
    norm = N.sqrt((2*l+1)/(4*N.pi))*N.sqrt(float(N.math.factorial(l-m))/float(N.math.factorial(l+m)))
    if m == 0:
        ffi = norm
    else:
        expimfi = (cosfi+1.0j*sinfi)**m # Find sin(m fi) and cos(m fi) as im and re parts
        if m < 0:
            norm = -(-1)**(-m)*N.sqrt(2)*norm
            ffi = norm * expimfi.imag
        else:
            norm = N.sqrt(2)*norm
            ffi = norm * expimfi.real
    return SS.lpmv(m, l, costh)*ffi


def _OLD_sphericalHarmonics(sinth, costh, sinfi, cosfi):
    pi = 3.141592654

    # l=0 m=0
    Y00 = 1/(2.*N.sqrt(pi))
    # l=1 m=-1
    Y1m1 = -(N.sqrt(3/pi)*sinfi*sinth)/2.
    # l=1 m=0
    Y10 = (costh*N.sqrt(3/pi))/2.
    # l=1 m=1
    Y11 = -(cosfi*N.sqrt(3/pi)*sinth)/2.
    # l=2 m=-2
    Y2m2 = (cosfi*N.sqrt(15/pi)*sinfi)/4. - \
           (cosfi*costh**2*N.sqrt(15/pi)*sinfi)/4. + \
           (cosfi*N.sqrt(15/pi)*sinfi*sinth**2)/4.
    # l=2 m=-1
    Y2m1 = -(costh*N.sqrt(15/pi)*sinfi*sinth)/2.
    # l=2 m=0
    Y20 = N.sqrt(5/pi)/8. + (3*costh**2*N.sqrt(5/pi))/8. - (3*N.sqrt(5/pi)*sinth**2)/8.
    # l=2 m=1
    Y21 = -(cosfi*costh*N.sqrt(15/pi)*sinth)/2.
    # l=2 m=2
    Y22 = (cosfi**2*N.sqrt(15/pi))/8. - (cosfi**2*costh**2*N.sqrt(15/pi))/ \
          8. - (N.sqrt(15/pi)*sinfi**2)/8. + (costh**2*N.sqrt(15/pi)*sinfi**2)/ \
          8. + (cosfi**2*N.sqrt(15/pi)*sinth**2)/8. - (N.sqrt(15/pi)*sinfi**2*sinth**2)/    8.
    # l=3 m=-3
    Y3m3 = (-9*cosfi**2*N.sqrt(35/(2.*pi))*sinfi*sinth)/ \
           16. + (9*cosfi**2*costh**2*N.sqrt(35/(2.*pi))*sinfi*sinth)/ \
           16. + (3*N.sqrt(35/(2.*pi))*sinfi**3*sinth)/    16. - (3*costh**2*N.sqrt(35/(2.*pi))*sinfi**3*sinth)/  \
           16. - (3*cosfi**2*N.sqrt(35/(2.*pi))*sinfi*sinth**3)/    16. + (N.sqrt(35/(2.*pi))*sinfi**3*sinth**3)/16.
    # l=3 m=-2
    Y3m2 = (cosfi*costh*N.sqrt(105/pi)*sinfi)/8. - (cosfi*costh**3*N.sqrt(105/pi)*sinfi)/    8. + (3*cosfi*costh*N.sqrt(105/pi)*sinfi*sinth**2)/8.
    # l=3 m=-1
    Y3m1 = -(N.sqrt(21/(2.*pi))*sinfi*sinth)/    16. - (15*costh**2*N.sqrt(21/(2.*pi))*sinfi*sinth)/    16. + (5*N.sqrt(21/(2.*pi))*sinfi*sinth**3)/16.
    # l=3 m=0
    Y30 = (3*costh*N.sqrt(7/pi))/16. + (5*costh**3*N.sqrt(7/pi))/    16. - (15*costh*N.sqrt(7/pi)*sinth**2)/16.
    # l=3 m=1
    Y31 = -(cosfi*N.sqrt(21/(2.*pi))*sinth)/    16. - (15*cosfi*costh**2*N.sqrt(21/(2.*pi))*sinth)/    16. + (5*cosfi*N.sqrt(21/(2.*pi))*sinth**3)/16.
    # l=3 m=2
    Y32 = (cosfi**2*costh*N.sqrt(105/pi))/16. - (cosfi**2*costh**3*N.sqrt(105/pi))/ \
          16. - (costh*N.sqrt(105/pi)*sinfi**2)/    16. + (costh**3*N.sqrt(105/pi)*sinfi**2)/ \
          16. + (3*cosfi**2*costh*N.sqrt(105/pi)*sinth**2)/    16. - (3*costh*N.sqrt(105/pi)*sinfi**2*sinth**2)/16.
    # l=3 m=3
    Y33 = (-3*cosfi**3*N.sqrt(35/(2.*pi))*sinth)/    16. + (3*cosfi**3*costh**2*N.sqrt(35/(2.*pi))*sinth)/ \
          16. + (9*cosfi*N.sqrt(35/(2.*pi))*sinfi**2*sinth)/    16. - (9*cosfi*costh**2*N.sqrt(35/(2.*pi))*sinfi**2*sinth)/ \
          16. - (cosfi**3*N.sqrt(35/(2.*pi))*sinth**3)/    16. + (3*cosfi*N.sqrt(35/(2.*pi))*sinfi**2*sinth**3)/16.

    return [[Y00],
            [Y1m1, Y10, Y11],
            [Y2m2, Y2m1, Y20, Y21, Y22],
            [Y3m3, Y3m2, Y3m1, Y30, Y31, Y32, Y33]]
