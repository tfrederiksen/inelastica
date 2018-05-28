import numpy as N
import numpy.fft as FFT


def Hilbert(f, ker=None):
    r"""
    Hilbert transform :math:`\mathcal{H}[f](y)` of a function :math:`f(x)`
    with finite support, sampled on an equidistant grid :math:`\{x_1,x_2,x_{nh}\}`.

    Definition:

    :math:`\mathcal{H}[f](y)= \frac{1}{\pi} p.v.\int_{-\infty}^{\infty} dx \frac{f(x)}{x-y}`.

    Parameters
    ----------
    f : ndarray
        List of function values :math:`f(x)`.
    ker : ndarray (optional)

    Returns
    -------
    Hf : ndarray
        Hilbert transform of :math:`f(x)`.
    ker : ndarray
        Kernel function, can be reused for other transformations with the same number of grid points.
    """

    def kernel(f):
        'Hilbert transform kernel'
        nh = len(f)
        n = 2 * nh
        aux = N.empty(nh+1, N.float)
        aux[0] = 0
        tmp = N.arange(1, nh+1)
        aux[1:] = N.log(tmp) * tmp
        ker = N.empty(n, N.float)
        ker[0] = 0
        ker[1:nh] = aux[2:nh+1] - 2*aux[1:nh] + aux[:nh-1]
        ker[nh+1:] = -ker[1:nh][::-1]
        return -FFT.fft(ker)/N.pi

    def transform(f, ker):
        'Convolution with kernel'
        n = len(f)
        fpad = FFT.fft(N.array((f, N.zeros(n))).flat)
        r = FFT.ifft(fpad*ker)
        return r[0:n]

    if ker != None:
        # A kernel was specified at the function call
        return transform(f, ker), ker
    else:
        print 'Hilbert: Generating kernel'
        ker = kernel(f)
        return transform(f, ker), ker
