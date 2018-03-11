import numpy as N
import numpy.linalg as LA


def outerAdd(* args):
    """
    Returns :math:`A_{ijk}=B_i+C_j+D_k`.
    """
    # A_ijk=B_i+C_j+D_k
    tmp = args[0].copy()
    for ii in range(1, len(args)):
        tmp = N.add.outer(tmp, args[ii])
    return tmp


def dist(x):
    """
    Euclidian norm :math:`||x||` of a vector :math:`x`.
    """
    return N.sqrt(N.dot(x, x))


def mysqrt(x):
    "Square root of matrix."
    ev, U = LA.eig(x)
    U = N.transpose(U)

    tmp = N.zeros((len(ev), len(ev)), N.complex)
    for ii,evi in enumerate(ev):
        tmp[ii, ii] = N.sqrt(evi)

    return mm(LA.inv(U), tmp, U)


def fermi(mu, E, kT):
    """
    Fermi function :math:`n_F(\mu) = 1/[exp((E-\mu)/k_BT)+1]`.
    """
    return 1/(N.exp(N.clip((E-mu)/kT, -70.0, 70.0))+1)


def box(mu1, mu2, grid, kT):
    """
    Window (box) function defined as
    :math:`n_F(\mu_2)-n_F(\mu_1)`
    """
    # f2-f1 (Box!)
    return fermi(mu2, grid, kT)-fermi(mu1, grid, kT)


def trapez(x, f, equidistant=False):
    """
    Integration of vector f on grid x using the 3rd degree polynomial.
    The grid x does not have to be equidistant.

    Parameters
    ----------
    x : ndarray
    f : ndarray
    equidistant : bool
        `False` = 3rd degree polynomial method.
        `True` = Linear trapez method.
    Returns
    -------
    res : complex number
        Result of the integration
    """
    if equidistant:
        # Trapez method!
        d = N.array((x[1]-x[0])*N.ones(len(x)), N.complex)
        d[0] = d[0]/2
        d[-1] = d[-1]/2
        return N.dot(d, f)
    else:
        # 3rd degree polynomial except for 1st and last bins
        res = (x[1]-x[0])*(f[0]+f[1])/2+(x[-1]-x[-2])*(f[-1]+f[-2])/2
        for ii in range(1, len(x)-2):
            x0, x1, x2, x3 = x[ii-1], x[ii], x[ii+1], x[ii+2]
            y0, y1, y2, y3 = f[ii-1], f[ii], f[ii+1], f[ii+2]
            res += ((x1-x2)*(-6*x0**2*(x0-x3)*x3**2*(y1+y2)+4*x0*x2*(x0-x3)*x3*(x0+x3)*(2*y1+y2)+\
                 3*x2**3*(x3**2*(y0-y1)+x0**2*(y1-y3))+\
                 x1**3*(-2*x2*(x3*(y0-y2)+x0*(y2-y3))+3*(x3**2*(y0-y2)+x0**2*(y2-y3))+\
                 x2**2*(-y0+y3))+x2**4*(x3*(-y0+y1)+x0*(-y1+y3))+\
                 x2**2*(-2*x3**3*(y0-y1)-3*x0**2*x3*(3*y1+y2)+3*x0*x3**2*(3*y1+y2)+\
                 2*x0**3*(-y1+y3))+x1**4*(x3*(-y0+y2)+x2*(y0-y3)+x0*(-y2+y3))+\
                 x1*(4*x0*(x0-x3)*x3*(x0+x3)*(y1+2*y2)-2*x2**3*(x3*(y0-y1)+x0*(y1-y3))+\
                 x2**4*(y0-y3)+x2*(-6*x0**2*x3*(y1+y2)+6*x0*x3**2*(y1+y2)+\
                 4*x3**3*(y0+y1+y2)-4*x0**3*(y1+y2+y3))-\
                 3*x2**2*(x3**2*(y0+2*y1+y2)-x0**2*(2*y1+y2+y3)))+\
                 x1**2*(-2*x3**3*(y0-y2)-3*x0**2*x3*(y1+3*y2)+3*x0*x3**2*(y1+3*y2)+\
                 x2**3*(-y0+y3)+2*x0**3*(-y2+y3)-\
                 3*x2*(x3**2*(y0+y1+2*y2)-x0**2*(y1+2*y2+y3))+\
                 3*x2**2*(x3*(2*y0+y1+y2)-x0*(y1+y2+2*y3)))))/\
                 (12.*(x0-x1)*(x0-x2)*(x0-x3)*(x1-x3)*(x2-x3))
        return res


def interpolate(nx, x, y):
    """
    Interpolate :math:`f(x)=y` to find :math:`f(nx)`. Extrapolation allowed.

    NB: Makes no checks for nx inside x region!!!

    Parameters
    ----------
    nx : ndarray
        Abscissa values to interpolate the function.
    x : ndarray
        Known abscissa values.
    y : ndarray
        Known function values.

    Returns
    -------
    ny : ndarray
        Interpolated function values on nx.
    """
    ny = N.array([0.0]*len(nx))
    Lpos = N.searchsorted(x, nx)
    for ix, pos in enumerate(Lpos):
        if pos < len(x):
            ny[ix] = y[pos-1]+(y[pos]-y[pos-1])/(x[pos]-x[pos-1])*(nx[ix]-x[pos-1])
        else:
            #TF: NB EXTRAPOLATION condition added!
            ny[ix] = y[-1]+(y[-1]-y[-2])/(x[-1]-x[-2])*(nx[ix]-x[-1])
    return ny
