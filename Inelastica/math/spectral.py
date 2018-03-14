import numpy as N
import numpy.linalg as LA


def mm(* args):
    """
    Matrix multiplication with arbitrary number of arguments
    and the SpectralMatrix type.
    """
    args = list(args)

    # Look for SpectralMatrices
    where = N.where(N.array([isinstance(ii, SpectralMatrix) for ii in args]))[0]
    if len(where) > 0:
        res = __mmSpectralMatrix(args, where)
    else:
        res = __mm(args)

    return res


def __mmSpectralMatrix(args, where):
    # find smallest
    size = [args[where[ii]].L.shape[1] for ii in range(len(where))]
    smallest = N.where(N.array(size) == N.min(size))[0][0]

    # if more than one SpectralMatrix, expand others into normal matrices
    args = [[ii] for ii in args] # To get flatten to work
    for ii in where:
        if ii != where[smallest]:
            args[ii] = [args[ii][0].L, args[ii][0].R]
    args = sum(args, []) # Flatten list

    # Split at Spectral matrix
    where = N.where(N.array([isinstance(ii, SpectralMatrix) for ii in args]))[0]
    res = SpectralMatrix()
    res.L = __mm(args[:where[0]]+[args[where[0]].L])
    res.R = __mm([args[where[0]].R]+args[where[0]+1:])
    return res


def __mm(args):
    # Normal matrix mult, order important for speed if matrices has different sizes
    if len(args) == 1: return args[0]

    # Find smallest matrix
    Lsize = N.array([ii.shape[0] for ii in args[:-1]])
    if len(args[-1].shape) == 1:
        # Most rightmost is vector
        Rsize = N.array([ii.shape[1] for ii in args[1:-1]]+[1])
    else:
        Rsize = N.array([ii.shape[1] for ii in args[1:]])

    Lwhere = N.where(Lsize==N.min(Lsize))[0]
    Rwhere = N.where(Rsize==N.min(Rsize))[0]

    if N.min(Lsize) > N.min(Rsize):
        where = [ii+1 for ii in Rwhere]
        left = False
    else:
        where = Lwhere
        left = True

    if len(where) == len(args)-1 or len(args) == 2:
        # Order does not matter
        res = N.dot(args[0], args[1])
        for ii in range(len(args)-2):
            res = N.dot(res, args[ii+2])
    else:
        # Make sure the matrix mult is done in order of smaller
        # matrices first
        ii = where[0]
        if ii == 0:
            res = __mm([__mm([args[0], args[1]])]+args[2:])
        elif ii == len(args)-1:
            res = __mm(args[:ii-1]+[__mm([args[ii-1], args[ii]])])
        else:
            if left:
                res = __mm(args[:ii]+[__mm([args[ii]]+args[ii+1:])])
            else:
                res = __mm([__mm(args[:ii]+[args[ii]])]+args[ii+1:])
    return res


class SpectralMatrix(object):
    r"""
    Matrix class for spectral matrices.

    Idea: Split a spectral matrix :math:`\mathbf{A}^{(mxm)}` into two smaller matrices
    :math:`\mathbf{A} = \mathbf{L}.\mathbf{R}`
    where :math:`\mathbf{L}^{(mxn)} = [\lambda_1 \mathbf{v}_1, ..., \lambda_n \mathbf{v}_n]`
    and :math:`\mathbf{R}^{(nxm)} = [\mathbf{v}_1, ..., \mathbf{v}_n]^\dagger`
    are constructed from a call to ``LA.eigh(A)`` keeping only
    :math:`n` eigensolutions with eigenvalues above the specified cutoff.

    The magic takes place in Matrix multiply function.

    Todo: Adding two gives full matrix! Should be easy to fix.

    Attributes
    ----------
    L : ndarray
    R : ndarray

    Parameters
    ----------
    A : ndarray
    """
    # self.L/R : Left / right matrices

    def __init__(self, A=None, cutoff=1e-8):
        if isinstance(A, N.ndarray):
            # Initialize ... only Hermitian matrices
            ev, evec = LA.eigh(A)
            # Drop eigenvalues
            indx = N.where(N.abs(ev)>cutoff)[0]
            ev, evec = ev[indx], evec[:, indx]
            print "SpectralMatrix: Fraction of eigenvalues above cutoff (%.1e) is %i/%i"%(cutoff, len(indx), len(A))
            self.L = N.dot(evec, N.diag(ev))
            self.R = dagger(evec)
            if False:
                print N.allclose(A, N.dot(self.L, self.R))

    def full(self):
        r"""
        Returns the dense ndarray via matrix multiplication :math:`\mathbf{A}^{(mxm)} = \mathbf{L}.\mathbf{R}`.
        """
        return mm(self.L, self.R)

    def __add__(self, b, subtract=False):
        # Could be improved for addition of two spectral matrices
        if not isinstance(b, SpectralMatrix):
            return self.full()+b
        else:
            NN = self.L.shape[0]
            Na, Nb = self.L.shape[1], b.L.shape[1]
            if Na+Nb>0.3*NN:
                # Too many eigenvalues
                print "Too many eigenvalues"
                if not subtract:
                    return self.full()+b.full()
                else:
                    return self.full()-b.full()
            else:
                res=SpectralMatrix()
                res.L = N.zeros((NN, Na+Nb), N.complex) # Always assume complex !?
                res.R = N.zeros((Na+Nb, NN), N.complex) # Always assume complex !?
                res.L[:, 0:Na], res.R[0:Na, :] = self.L, self.R
                if not subtract:
                    res.L[:, Na:Na+Nb], res.R[Na:Na+Nb, :] = b.L, b.R
                else:
                    res.L[:, Na:Na+Nb], res.R[Na:Na+Nb, :] = -b.L, b.R
                return res

    def __radd__(self, b):
        return self+b

    def __sub__(self, b):
        if not isinstance(b, SpectralMatrix):
            return self.full()-b
        else:
            return self.__add__(b, subtract=True)

    def __rsub__(self, b):
        if not isinstance(b, SpectralMatrix):
            return b-self.full()
        else:
            return b.__add__(self, subtract=True)

    def __mul__(self, b):
        tmp = SpectralMatrix()
        tmp.L, tmp.R = self.L*b, self.R
        return tmp

    def __rmul__(self, b):
        return self*b

    def __dagger__(self):
        tmp = SpectralMatrix()
        tmp.L = dagger(self.R)
        tmp.R = dagger(self.L)
        return tmp


def trace(a):
    """
    Returns the trace of a normal or spectral matrix.
    """
    if isinstance(a, SpectralMatrix):
        return N.trace(mm(a.R, a.L)) # Switch sum around to make faster!
    else:
        return N.trace(a)


def dagger(x):
    """
    Returns the hermitian conjugation of a normal or spectral matrix.
    """
    if isinstance(x, SpectralMatrix):
        return x.__dagger__()
    else:
        return N.transpose(N.conjugate(x))
