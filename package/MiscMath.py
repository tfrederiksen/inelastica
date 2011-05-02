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
