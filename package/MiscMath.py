print "SVN $Id$"

import numpy as N
import numpy.linalg as LA
import numpy.fft as FFT

##############################################################

def mm(* args):
    # Matrix multiplication with arbitrary number of arguments                                                                                              
    tmp=N.dot(args[0],args[1])
    for ii in range(len(args)-2):
        tmp=N.dot(tmp,args[ii+2])
    return tmp

def dagger(x):
    # Hermitian conjugation
    return N.transpose(N.conjugate(x))

##############################################################

def outerAdd(* args):
    # A_ijk=B_i+C_j+D_k
    tmp=args[0].copy()
    for ii in range(1,len(args)):
        tmp=N.add.outer(tmp,args[ii])
    return tmp

def dist(x):
    return N.sqrt(N.dot(x,x))

def mysqrt(x):
    # Square root of matrix    
    ev,U = LA.eig(x)
    U = N.transpose(U)

    tmp=N.zeros((len(ev),len(ev)),N.complex)
    for ii in range(len(ev)):
        tmp[ii,ii]=N.sqrt(ev[ii])
        
    return mm(LA.inv(U),tmp,U)

##############################################################

def fermi(mu,E,kT):
    return 1/(N.exp(N.clip((E-mu)/kT,-70.0,70.0))+1)

def box(mu1,mu2,grid,kT):
    # f2-f1 (Box!)
    return fermi(mu2,grid,kT)-fermi(mu1,grid,kT)

##############################################################
def trapez(x,f,equidistant=False):
    '''
    Integration of vector f on grid x using the 3rd degree polynomial.
    The grid x does not have to be equidistant.
    '''
    if equidistant:
        # Trapez method!
        d = N.array((x[1]-x[0])*N.ones(len(x)),N.complex)
        d[0] = d[0]/2
        d[-1] = d[-1]/2
        return N.dot(d,f)
    else:
        # 3rd degree polynomial except for 1st and last bins
        sum=(x[1]-x[0])*(f[0]+f[1])/2+(x[-1]-x[-2])*(f[-1]+f[-2])/2
        for ii in range(1,len(x)-2):
            x0,x1,x2,x3=x[ii-1],x[ii],x[ii+1],x[ii+2]
            y0,y1,y2,y3=f[ii-1],f[ii],f[ii+1],f[ii+2]
            sum+=((x1-x2)*(-6*x0**2*(x0-x3)*x3**2*(y1+y2)+4*x0*x2*(x0-x3)*x3*(x0+x3)*(2*y1+y2)+\
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
        return sum
##############################################################

def interpolate(nx,x,y):
    """
    Interpolate f(x)=y to find f(nx)
    Makes no checks for nx inside x region!!!
    """
    ny=N.array([0.0]*len(nx))
    Lpos=N.searchsorted(x,nx)
    for ix,pos in enumerate(Lpos):
        if pos<len(x):
            ny[ix]=y[pos-1]+(y[pos]-y[pos-1])/(x[pos]-x[pos-1])*(nx[ix]-x[pos-1])
        else:
            #TF: NB EXTRAPOLATION condition added!
            ny[ix]=y[-1]+(y[-1]-y[-2])/(x[-1]-x[-2])*(nx[ix]-x[-1])
    return ny

##############################################################
# DEFINITION OF THE HILBERT TRANSFORM:
# H[f](y) = 1/\pi p.v.\int^{\infty}_{\infty} dx { f(x)/(x-y) }
def Hilbert(f,ker=None):
    'Hilbert transform'
    
    def kernel(f):
        'Hilbert transform kernel'
        n = 2*len(f)
        aux = N.zeros(n/2+1,N.float)
        for i in N.arange(1,n/2+1):
            aux[i] = i*N.log(i)
        ker = N.zeros(n,N.float)
        for i in N.arange(1,n/2):
            ker[i] = aux[i+1]-2*aux[i]+aux[i-1]
            ker[n-i] = -ker[i]
        return -FFT.fft(ker)/N.pi

    def transform(f,ker):
        'Convolution with kernel'
        n = len(f)
        fpad = FFT.fft(N.array((f,N.zeros(n))).flat)
        r = FFT.ifft(fpad*ker)
        return r[0:n]

    if ker!=None:
        # A kernel was specified at the function call
        return transform(f,ker), ker
    else:
        print 'Hilbert: Generating kernel'
        ker = kernel(f)
        return transform(f,ker), ker

##############################################################

def sphericalHarmonics(sinth,costh,sinfi,cosfi):
    pi=3.141592654

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
            [Y1m1,Y10,Y11],
            [Y2m2,Y2m1,Y20,Y21,Y22],
            [Y3m3,Y3m2,Y3m1,Y30,Y31,Y32,Y33]]

##############################################################
# GaussKronrod

def abwe1 ( n, m, tol, coef2, even, b, x ) :
#*****************************************************************************80
## ABWE1 calculates a Kronrod abscissa and weight.
#  Licensing:
#    This code is distributed under the GNU LGPL license.
#  Modified:
#    30 April 2013
#  Author:
#    Original FORTRAN77 version by Robert Piessens, Maria Branders.
#    Python version by John Burkardt.
#  Reference:
#    Robert Piessens, Maria Branders,
#    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
#    of Gauss and Lobatto,
#    Mathematics of Computation,
#    Volume 28, Number 125, January 1974, pages 135-139.
  from sys import exit
  if ( x == 0.0 ):    ka = 1
  else:    ka = 0
  for iter in range ( 1, 51 ):
    b1 = 0.0
    b2 = b[m+1-1]
    yy = 4.0 * x * x - 2.0
    d1 = 0.0
    if ( even ):
      ai = m + m + 1
      d2 = ai * b[m+1-1]
      dif = 2.0
    else:
      ai = m + 1
      d2 = 0.0
      dif = 1.0
    for k in range ( 1, m + 1 ):
      ai = ai - dif
      i = m - k + 1
      b0 = b1
      b1 = b2
      d0 = d1
      d1 = d2
      b2 = yy * b1 - b0 + b[i-1]
      if ( not even ):
        i = i + 1
      d2 = yy * d1 - d0 + ai * b[i-1]

    if ( even ):
      f = x * ( b2 - b1 )
      fd = d2 + d1
    else:
      f = 0.5 * ( b2 - b0 )
      fd = 4.0 * x * d2
#
#  Newton correction.
#
    delta = f / fd
    x = x - delta
    if ( ka == 1 ):
      break
    if ( abs ( delta ) <= tol ):
      ka = 1
#
#  Catch non-convergence.
#
  if ( ka != 1 ):
    print ''
    print 'ABWE1 - Fatal error!'
    print '  Iteration limit reached.'
    print '  Last DELTA was %e' % ( delta )
    sys.exit ( 'ABWE1 - Fatal error!' )
#
#  Computation of the weight.
#
  d0 = 1.0
  d1 = x
  ai = 0.0
  for k in range ( 2, n + 1 ):
    ai = ai + 1.0
    d2 = ( ( ai + ai + 1.0 ) * x * d1 - ai * d0 ) / ( ai + 1.0 )
    d0 = d1
    d1 = d2
  w = coef2 / ( fd * d2 )
  return x, w

def abwe2 ( n, m, tol, coef2, even, b, x ):
#*****************************************************************************80
## ABWE2 calculates a Gaussian abscissa and two weights.
#  Licensing:
#    This code is distributed under the GNU LGPL license.
#  Modified:
#    30 April 2013
#  Author:
#    Original FORTRAN77 version by Robert Piessens, Maria Branders.
#    PYTHON version by John Burkardt.
#  Reference:
#    Robert Piessens, Maria Branders,
#    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
#    of Gauss and Lobatto,
#    Mathematics of Computation,
#    Volume 28, Number 125, January 1974, pages 135-139.
  from sys import exit
  if ( x == 0.0 ):    ka = 1
  else:    ka = 0
  for iter in range ( 1, 51 ):
    p0 = 1.0
    p1 = x
    pd0 = 0.0
    pd1 = 1.0
    if ( n <= 1 ):
      if ( x != 0.0 ):
        p2 = ( 3.0 * x * x - 1.0 ) / 2.0
        pd2 = 3.0 * x
      else:
        p2 = 3.0 * x
        pd2 = 3.0

    ai = 0.0
    for k in range ( 2, n + 1 ):
      ai = ai + 1.0
      p2 = ( ( ai + ai + 1.0 ) * x * p1 - ai * p0 ) / ( ai + 1.0 )
      pd2 = ( ( ai + ai + 1.0 ) * ( p1 + x * pd1 ) - ai * pd0 ) / ( ai + 1.0 )
      p0 = p1
      p1 = p2
      pd0 = pd1
      pd1 = pd2
    delta = p2 / pd2
    x = x - delta
    if ( ka == 1 ):
      break
    if ( abs ( delta ) <= tol ):
      ka = 1
  if ( ka != 1 ):
    print ''
    print 'ABWE2 - Fatal error!'
    print '  Iteration limit reached.'
    print '  Last DELTA was %e' % ( delta )
    sys.exit ( 'ABWE2 - Fatal error!' )
  an = n
  w2 = 2.0 / ( an * pd2 * p0 )
  p1 = 0.0
  p2 = b[m+1-1]
  yy = 4.0 * x * x - 2.0
  for k in range ( 1, m + 1 ):
    i = m - k + 1
    p0 = p1
    p1 = p2
    p2 = yy * p1 - p0 + b[i-1]
  if ( even ):
    w1 = w2 + coef2 / ( pd2 * x * ( p2 - p1 ) )
  else:
    w1 = w2 + 2.0 * coef2 / ( pd2 * ( p2 - p0 ) )
  return x, w1, w2

def kronrod ( n, tol ):
#*****************************************************************************80
#
## KRONROD adds N+1 points to an N-point Gaussian rule.
#
#  Discussion:
#
#    This subroutine calculates the abscissas and weights of the 2N+1
#    point Gauss Kronrod quadrature formula which is obtained from the
#    N point Gauss quadrature formula by the optimal addition of N+1 points.
#
#    The optimally added points are called Kronrod abscissas.  The
#    abscissas and weights for both the Gauss and Gauss Kronrod rules
#    are calculated for integration over the interval [-1,+1].
#
#    Since the quadrature formula is symmetric with respect to the origin,
#    only the nonnegative abscissas are calculated.
#
#    Note that the code published in Mathematics of Computation
#    omitted the definition of the variable which is here called COEF2.
#
#  Storage:
#
#    Given N, let M = ( N + 1 ) / 2.
#
#    The Gauss-Kronrod rule will include 2*N+1 points.  However, by symmetry,
#    only N + 1 of them need to be listed.
#
#    The arrays X, W1 and W2 contain the nonnegative abscissas in decreasing
#    order, and the weights of each abscissa in the Gauss-Kronrod and
#    Gauss rules respectively.  This means that about half the entries
#    in W2 are zero.
#
#    For instance, if N = 3, the output is:
#
#    I      X               W1              W2
#
#    1    0.960491        0.104656         0.000000
#    2    0.774597        0.268488         0.555556
#    3    0.434244        0.401397         0.000000
#    4    0.000000        0.450917         0.888889
#
#    and if N = 4, (notice that 0 is now a Kronrod abscissa)
#    the output is
#
#    I      X               W1              W2
#
#    1    0.976560        0.062977        0.000000
#    2    0.861136        0.170054        0.347855
#    3    0.640286        0.266798        0.000000
#    4    0.339981        0.326949        0.652145
#    5    0.000000        0.346443        0.000000
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    30 April 2013
#
#  Author:
#
#    Original FORTRAN77 version by Robert Piessens, Maria Branders.
#    PYTHON version by John Burkardt.
#
#  Reference:
#
#    Robert Piessens, Maria Branders,
#    A Note on the Optimal Addition of Abscissas to Quadrature Formulas
#    of Gauss and Lobatto,
#    Mathematics of Computation,
#    Volume 28, Number 125, January 1974, pages 135-139.
#
#  Parameters:
#
#    Input, integer N, the order of the Gauss rule.
#
#    Input, real TOL, the requested absolute accuracy of the
#    abscissas.
#
#    Output, real X(N+1), the abscissas.
#
#    Output, real W1(N+1), the weights for the Gauss-Kronrod rule.
#
#    Output, real W2(N+1), the weights for the Gauss rule.
#
  from numpy import pi
  from numpy import zeros
  from math import floor
  from math import sin
  from math import sqrt

  tau_dim = ( n + 1 ) / 2
  tau_dim = int ( floor ( tau_dim ) )
  b = zeros ( ((n+1)/2)+1 )
  tau = zeros ( tau_dim )
  w1 = zeros ( n + 1 )
  w2 = zeros ( n + 1 )
  x = zeros ( n + 1 )
  m = int ( floor ( ( n + 1 ) / 2 ) )
  even = ( 2 * m == n )
  d = 2.0
  an = 0.0
  for k in range ( 1, n + 1 ):
    an = an + 1.0
    d = d * an / ( an + 0.5 )
  tau[1-1] = ( an + 2.0 ) / ( an + an + 3.0 )
  b[m-1] = tau[1-1] - 1.0
  ak = an
  for l in range ( 1, m ):
    ak = ak + 2.0
    tau[l+1-1] = ( ( ak - 1.0 ) * ak \
      - an * ( an + 1.0 ) ) * ( ak + 2.0 ) * tau[l-1] \
      / ( ak * ( ( ak + 3.0 ) * ( ak + 2.0 ) \
      - an * ( an + 1.0 ) ) )
    b[m-l-1] = tau[l+1-1]
    for ll in range ( 1, l + 1 ):
      b[m-l-1] = b[m-l-1] + tau[ll-1] * b[m-l+ll-1]
  b[m+1-1] = 1.0
  bb = sin ( 0.5 * pi / ( an + an + 1.0 ) )
  x1 = sqrt ( 1.0 - bb * bb )
  s = 2.0 * bb * x1
  c = sqrt ( 1.0 - s * s )
  coef = 1.0 - ( 1.0 - 1.0 / an ) / ( 8.0 * an * an )
  xx = coef * x1
  coef2 = 2.0 / ( 2 * n + 1 )
  for i in range ( 1, n + 1 ):
    coef2 = coef2 * 4.0 * i / ( n + i )
  for k in range ( 1, n + 1, 2 ):
    [ xx, w1[k-1] ] = abwe1 ( n, m, tol, coef2, even, b, xx )
    w2[k-1] = 0.0
    x[k-1] = xx
    y = x1
    x1 = y * c - bb * s
    bb = y * s + bb * c
    if ( k == n ):
      xx = 0.0
    else:
      xx = coef * x1
    [ xx, w1[k+1-1], w2[k+1-1] ] = abwe2 ( n, m, tol, coef2, even, b, xx )
    x[k+1-1] = xx
    y = x1
    x1 = y * c - bb * s
    bb = y * s + bb * c
    xx = coef * x1
  if ( even ):
    xx = 0.0
    [ xx, w1[n+1-1] ] = abwe1 ( n, m, tol, coef2, even, b, xx )
    w2[n+1-1] = 0.0
    x[n+1-1] = xx
  return x, w1, w2

def GaussKronrod(NN):
    # Rescale to [-0.5,0.5] and expand to get all points.
    x, w1, w2 = kronrod(NN,1e-10)
    indx = N.argsort(x)
    x, w1, w2 = x[indx], w1[indx], w2[indx]
    if N.abs(x[0])<1e-7:
        x, w1, w2 = N.concatenate((-x,x[1:])), N.concatenate((w1,w1[1:])), N.concatenate((w2,w2[1:]))
    else:
        x, w1, w2 = N.concatenate((-x,x)), N.concatenate((w1,w1)), N.concatenate((w2,w2))
    indx = N.argsort(x)
    x, w1, w2 = x[indx], w1[indx], w2[indx]
    return x/2, w1/2, w2/2

##############################################################
