import numpy as N
import numpy.linalg as LA
import numpy.fft as FFT


# DEFINITION OF THE HILBERT TRANSFORM:
# H[f](y) = 1/\pi p.v.\int^{\infty}_{\infty} dx { f(x)/(x-y) }


def Hilbert(f, ker=None):
    'Hilbert transform'

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

    if ker!=None:
        # A kernel was specified at the function call
        return transform(f, ker), ker
    else:
        print 'Hilbert: Generating kernel'
        ker = kernel(f)
        return transform(f, ker), ker
