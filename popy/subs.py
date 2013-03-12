'''
a collection of commonly used routines

Created on Jan 28, 2013

@author: Jinbo Wang <jinbow@gmail.com>
@organization: Scripps Institute of Oceanography
'''

def setlabel(ax,xlab='',ylab='',title=''):
    """
    set xlabel, ylabel, title for ax
    """
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    return

def savefig(i,frst='',figfn='',caption='',figpath='../figures/'):
    """
    save figure in png format
    write rst file
    """
    import sys,os,pylab
    sfn = sys.argv[0].split('.')[0]
    if frst == '':
        frst = sfn+'.rst'
    if figfn =='':
        figfn = figpath+sfn+'_%03i'%i+'.png'
    if not os.path.exists(frst):
        f = open(frst,'w')
        words = ['.. _'+sfn+':\n\n'+sfn+'\n=========================\n\n',
                 '.. note:: prepared by '+sys.argv[0]+'\n',
                 '.. literalinclude:: '+sys.argv[0] ]
        f.writelines(words)
        f.close()
    f = open(frst,'r')
    w0 = f.readlines()
    f.close()
    f = open(frst,'w')
    w1 = ['\n.. figure:: '+figfn+'\n',
          '    :width: 700px\n',
          '    :align: center\n\n',
          '    '+caption+ "\n\n",
          ':ref:`Top <'+sfn+'>`\n\n']
    f.writelines(w0[:-1])
    f.writelines(w1)
    f.writelines(w0[-1])
    pylab.savefig(figfn)
    pylab.clf()
    return i+1

def smooth2d(x,window_len=10,window='hanning'):
    from numpy import zeros,arange
    newny=(smooth(x[:,0],window_len,window).size)
    newnx=(smooth(x[0,:],window_len,window).size)
    t1=zeros((newny,newnx))
    t2=zeros((newny,newnx))
    for i in arange(newnx):
        t1[:,i]=smooth(x[:,i],window_len,window)
    for j in arange(newny):
        t2[j,:]=smooth(t1[j,:],window_len,window)
    return t2

def smooth(x,window_len=10,window='hanning'):
    import numpy
    """smooth the data using a window with requested size.
    copied from http://www.scipy.org/Cookbook/SignalSmooth -- Jinbo Wang

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

def yticknum(ax,n):
    import matplotlib.pylab as plt
    ax.yaxis.set_major_locator(plt.MaxNLocator(n))
    return
def xticknum(ax,n):
    import matplotlib.pylab as plt
    ax.xaxis.set_major_locator(plt.MaxNLocator(n))
    return
def noxtick(ax):
    import matplotlib.pylab as plt
    plt.setp(ax.get_xticklabels(),visible=False)
    return
def noytick(ax):
    import matplotlib.pylab as plt
    plt.setp(ax.get_yticklabels(),visible=False)
    return

def D2matrix(z, N2, f0):
    import numpy as np
    def twopave(x):
        return (x[0:-1]+x[1:])/2.
    """
    construct second-order centered difference matrix
    \frac{\partial}{\partial_z} \frac{f_0^2}{N^2}
    \frac{\partial}{\partial_z}
    """
    #f0 = 1e-4
    N = z.size
    N2c = f0 ** 2 / twopave(N2)  # N2c is the N2 at center point
    dzc = np.diff(z)
    dzf = np.r_[0, twopave(dzc)]
    A = np.zeros((N, N))
    for k in range(1, N - 1):
        A[k, k - 1] = - N2c[k - 1] / (dzc[k - 1] * dzf[k])
        A[k, k] = N2c[k - 1] / (dzc[k - 1] * dzf[k]) + \
                  N2c[k] / (dzc[k] * dzf[k])
        A[k, k + 1] = - N2c[k] / (dzc[k] * dzf[k])

    A[0, 0] = N2c[0] / dzc[0] ** 2
    A[0, 1] = -A[0, 0]
    A[-1, -1] = N2c[-1] / dzc[-1] ** 2
    A[-1, -2] = -A[-1, -1]
    return A


def qgdecomp(z,N2,f0):
    import numpy as np
    M = D2matrix(z, N2,f0)
    w, vi = np.linalg.eig(M)
    vi = vi[:, np.argsort(w)]
    vi = vi / vi[0,:].reshape(1,-1) #eigenfunction
    w = np.sort(w) #eigenvalue
    return w, vi

def anomaly(lon,lat,var):
    from numpy import meshgrid,zeros,arange
    if len(lon.shape) == 1:
        x1,y1=meshgrid(lon,lat)
    else:
        x1,y1 = lon, lat
    if len(var.shape) == 2:
        va,vm = fit2Dsurf(x1,y1,var)
        return va,vm
    elif len(var.shape) == 3:
        vatmp = zeros(var.shape)
        vmtmp = zeros(var.shape)
        for i in arange(var.shape[0]):
            va,vm = fit2Dsurf(x1,y1,var[i,:,:])[0]
            vatmp[i,:,:] = va
            vmtmp[i,:,:] = vm
        return vatmp, vmtmp

def fit2Dsurf(x,y,p):
    """
      given y0=f(t0), find the best fit
      p = a + bx + cy + dx**2 + ey**2 + fxy
      and return a,b,c,d,e,f
    """
    from scipy.optimize import leastsq
    import numpy as np
    x,y=abs(x),abs(y)
    def err(c,x0,y0,p):
        a,b,c,d,e,f=c
        return p - (a + b*x0 + c*y0 + d*x0**2 + e*y0**2 + f*x0*y0)
    def surface(c,x0,y0):
        a,b,c,d,e,f=c
        return a + b*x0 + c*y0 + d*x0**2 + e*y0**2 + f*x0*y0
    dpdy = (np.diff(p,axis=0)/np.diff(y,axis=0)).mean()
    dpdx = (np.diff(p,axis=1)/np.diff(x,axis=1)).mean()
    xf=x.flatten()
    yf=y.flatten()
    pf=p.flatten()
    c = [pf.mean(),dpdx,dpdy,1e-22,1e-22,1e-22]
    coef = leastsq(err,c,args=(xf,yf,pf))[0]
    vm = surface(coef,x,y) #mean surface
    va = p - vm #anomaly
    return va,vm

def fit2exp(x,y,method='exp'):
    """
      given y0=f(t0), find the best fit
      p = a + b*exp(c*z)
      and return a,b,c,d,e,f
    """
    from scipy.optimize import leastsq
    import numpy as np

    def fit2poly(cc,x,y):
        if method=='exp':
            err = y - ( cc[0]*np.exp(cc[1]*x))
        elif method=='tanh':
            err = y - ( cc[0]*(1. - np.tanh(cc[1]*x)**2))
        return err

    x=x.flatten()
    y=y.flatten()
    c = [1e-5,1./200.]
    coef = leastsq(fit2poly,c,args=(x,y))
    return coef

def fit2poly(c,x,y):
    a,b,c,d,e,f=c
    fit = (a + b*x + c*y + d*x**2 + e*y**2 + f*x*y)
    return fit

def lat2f(d):
    """ Calculate Coriolis parameter from latitude
    d: latitudes, 1- or 2-D
    """
    from math import pi, sin
    return 2.0*0.729e-4*sin(d*pi/180.0)

def lonlat2xy(lon,lat):
    """convert lat lon to y and x
    x, y = lonlat2xy(lon, lat)
    lon and lat are 1d variables.
    x, y are 2d meshgrids.
    """
    from pylab import meshgrid,cos,pi
    r = 6371.e3
    lon = lon-lon[0]
    if lon.ndim == 1:
        lon,lat = meshgrid(lon,lat)
    x = 2*pi*r*cos(lat*pi/180.)*lon/360.
    y = 2*pi*r*lat/360.
    return x,y


def xy2lonlat(x,y):
    from pylab import meshgrid,cos,pi
    r = 6371.e3
    if len(x.shape) == 1:
        x,y = meshgrid(x,y)
    lon = x*180./(pi*r*cos(y/r))
    lat = y*180./(pi*r)
    return lon,lat

        
def twopave(x):
    return (x[0:-1]+x[1:])/2.

def gradxy(lon,lat,ssh):
    from numpy import pi,c_,r_,cos,diff,zeros
    from pylab import meshgrid
    r = 6371.e3
    if lon.ndim == 1:
        lon,lat = meshgrid(lon,lat)
    x = 2*pi*r*cos(lat*pi/180.)*lon/360.
    y = 2*pi*r*lat/360.
    def uv(x,y,ssh):
        x = c_[x, x[:,-2]]
        y = r_[y, y[-2,:].reshape(1,-1)]
        
        sshx = c_[ssh, ssh[:,-1]]
        v = diff(sshx,axis=1)/diff(x,axis=1)
        v[:,-1] = v[:,-2]

        sshy = r_[ssh, ssh[-1,:].reshape(1,-1)]
        u = diff(sshy,axis=0)/diff(y,axis=0)
        u[-2,:]=u[-1,:]
        return u,v

    if len(ssh.shape)==2:
        return uv(x,y,ssh)
    elif len(ssh.shape)==3:
        u, v = zeros(ssh.shape), zeros(ssh.shape)
        for i in range(ssh.shape[0]):
            ut,vt = uv(x,y,ssh[i,:,:])
            v[i,:,:] = vt
            u[i,:,:] = ut
        return u, v

def ssh2uv(lon,lat,ssh):
    """ calculate geostrophic velocity from SSH data
    lon[Nx], lat[ny]
    ssh[lat,lon] """
    from numpy import pi,sin#,c_,r_,cos,diff
    from pylab import meshgrid
    g = 9.8; #r = 6371.e3
    omega = 0.729e-4
    lon,lat = meshgrid(lon,lat)
    f=2.0*omega*sin(lat*pi/180.0)
    psi = g/f * ssh
    u, v = gradxy(lon,lat,psi)
    u = -1*u
    return u,v



def findsubdomain(lon,lat,p):
    """find the index number for the subdomain
     p=[lon0,lon1,lat0,lat1] in lon[nx] and lat[ny] """
    from numpy import argmin
    lon0, lon1, lat0, lat1 = p
    i0,i1,j0,j1 = argmin(lon-lon0), argmin(lon-lon1), argmin(lat-lat0), argmin(lat-lat1)
    return [i0,i1+1,j0,j1+1]