'''
a collection of commonly used routines

Created on Jan 28, 2013

@author: Jinbo Wang <jinbow@gmail.com>
@organization: Scripps Institute of Oceanography
'''

def fourpointavg(p):
    """
    four points average.
    return p(ny-1,nx-1)
    """
    if p.ndim==2:
        return (p[:-1,:-1]+p[:-1,1:]+p[1:,:-1]+p[1:,1:])/4.
    elif p.ndim==3:
        return (p[:,:-1,:-1]+p[:,:-1,1:]+p[:,1:,:-1]+p[:,1:,1:])/4.
    
def loaddata(fns,varn = None, ismean=True):
    """load data from fn,
    automatically detect file format"""
    import cPickle, gzip
    import numpy as np
    #import numpy as np
    def loadone(fn):
        print "load data from "+fn
        fs = fn.split('.')
        if fs[-1] == 'gz':
            f=gzip.open(fn,'rb')
        else:
            f=open(fn,'rb')
        d = cPickle.load(f)
        f.close()
        return d
    if np.isscalar(fns):
        fns=[fns]
    
    def loopkeymean(d1,d2):
        if type(d1) == dict:
            for key in d1.keys():
                if np.ma.isMaskedArray(d1[key]):
                    d3=d1[key].data
                    d3[d1[key].mask]=0
                    
                    d4=d2[key].data
                    d4[d2[key].mask]=0
                    
                    d1[key]=d4+d3
                    del d3, d4
                else:
                    d1[key]+=d2[key]
                
            return d1
        else:
            return d1+d2
    def stacks(e1,e2):
        print e1.ndim, e2.ndim
        if e1.ndim - e2.ndim ==1:
            e1=np.r_[e1,e2[np.newaxis,:]]
        else:
            e1=np.r_[e1,e2]
        return e1
    
    def loopkeystack(d1,d2):
        if type(d1) == dict:
            for key in d1.keys():
                d1[key]=stacks(d1[key],d2[key])
            return d1
        else:
            return stacks(d1,d2)
    
    da = loadone(fns[0])
    
    for fnn in fns[1:]:
        
        da1=loadone(fnn)
        if ismean:
            da=loopkeymean(da,da1)  
        else:
            da = loopkeystack(da,da1)

    return da


def clean(path):
    import os 
    if os.path.exists(path):
        os.remove(path)
    return

def smooth2d(x,wx=10,wy=10,window='hanning'):
    from numpy import zeros,arange
    newny=(smooth(x[:,0],wy,window).size)
    newnx=(smooth(x[0,:],wx,window).size)
    t1=zeros((newny,newnx))
    t2=zeros((newny,newnx))
    for i in arange(newnx):
        t1[:,i]=smooth(x[:,i],wy,window)
    for j in arange(newny):
        t2[j,:]=smooth(t1[j,:],wx,window)
    return t2

def MaskedSmooth(data, xp=3, yp=2):
    from numpy import arange,c_
    dout = data.copy()

    data = c_[data[:, -xp:], data, data[:, :xp]]
    ny, nx = dout.shape

    for j in arange(yp, ny - yp):
        for i in arange(nx):
            if ~dout.mask[j, i]:
                dout[j, i] = data[j - yp:j + yp + 1, i:i + 2 * xp + 1].mean()
    return dout

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

    s=numpy.ma.array(numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]])
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]


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

def EllipticGaussian(x,y,coef):
    from numpy import exp
    """D exp(-(a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2))
    """
    a,b,c,d,x0,y0 = coef
    v= d*exp(-(a*(x-x0)**2+b*(x-x0)*(y-y0)+c*(y-y0)**2))
    #print v
    return v

def fit2Gaussian1d(x,v):
    from scipy.optimize import leastsq
    from numpy import ma
    
    def Gaussian1d(x,coef):
        from numpy import exp
        """a exp(-(x-x0)^2/sigma^2)
        """
        a,sigma,x0 = coef
        v= a*exp(-(x-x0)**2 / sigma**2)
        #print v
        return v
        
    def err(coef):
        if ma.is_masked(v):
            res = v.flatten() - Gaussian1d(x,coef).flatten()
            return res.data[res.mask==False]
        else:
            res = v.flatten() - Gaussian1d(x,coef).flatten()
            return res
        
    vm = v.sum()
    x0 = (x*v).sum()/vm
    sigma = ((x-x0)**2*v).sum()/vm
     
    coef = [v.max(), sigma**0.5, x0]
    coef_fit = leastsq(err,coef,epsfcn=1e-40)[0]
    
    v_fit = Gaussian1d(x,coef_fit)
    return coef_fit,v_fit

def fit2EllipticGaussian(x,y,v):
    from scipy.optimize import leastsq
    from numpy import ma
    def err(coef):
        if ma.is_masked(v):
            res = v.flatten() - EllipticGaussian(x,y,coef).flatten()
            return res.data[res.mask==False]
        else:
            res = v.flatten() - EllipticGaussian(x,y,coef).flatten()
            return res
        
    vm = v.sum()
    x0,y0 = (x*v).sum()/vm,(y*v).sum()/vm
    a = ((x-x0)**2*v).sum()/vm
    c = ((y-y0)**2*v).sum()/vm
     
    coef = [0.5/a, 0., 0.5/c, v.mean(), x0,y0]
    coef_fit = leastsq(err,coef,epsfcn=1e-40)[0]
    
    v_fit = EllipticGaussian(x,y,coef_fit)
    return coef_fit,v_fit

def lat2f(d):
    """ Calculate Coriolis parameter from latitude
    d: latitudes, 1- or 2-D
    """
    from numpy import pi, sin
    return 2.0*0.729e-4*sin(d*pi/180.0)

def lonlat2xy(lon,lat):
    """convert lat lon to y and x
    x, y = lonlat2xy(lon, lat)
    lon and lat are 1d variables.
    x, y are 2d meshgrids.
    """
    from pylab import meshgrid,cos,pi
    r = 6371.e3
    #lon = lon-lon[0]
    if lon.ndim == 1:
        lon,lat = meshgrid(lon,lat)
    x = 2*pi*r*cos(lat*pi/180.)*lon/360.
    y = 2*pi*r*lat/360.
    return x,y

def degree2km(lat):
    """ calculate degree to kilometers at a given latitude
    """
    from pylab import cos,pi
    r = 6371.e3
    return 2*pi*r*cos(lat*pi/180.)/360.

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

def fourpave(p):
    """
    four points average.
    return p(ny-1,nx-1)
    """
    if p.ndim==2:
        return (p[:-1,:-1]+p[:-1,1:]+p[1:,:-1]+p[1:,1:])/4.
    elif p.ndim==3:
        return (p[:,:-1,:-1]+p[:,:-1,1:]+p[:,1:,:-1]+p[:,1:,1:])/4.
    
def gradxy(lon,lat,ssh):
    """ return psi_y, psi_x
    """
    from numpy import pi,c_,r_,cos,diff,zeros
    from pylab import meshgrid
    r = 6371.e3
    if lon.ndim == 1:
        lon,lat = meshgrid(lon,lat)
    x = 2*pi*r*cos(lat*pi/180.)*lon/360.
    y = 2*pi*r*lat/360.
    def uv(x,y,ssh):
#         x = c_[x, x[:,-2]]
#         y = r_[y, y[-2,:].reshape(1,-1)]
#         
#        sshx = c_[ssh, ssh[:,-1]]
        sshx=ssh
        v = diff(sshx,axis=1)/diff(x,axis=1)
        vv = (v[:,0:-1]+v[:,1:])/2
        del v
        v=c_[vv[:,0:1],vv,vv[:,-1].reshape(-1,1)]
        del vv
        #v[:,-1] = v[:,-2]

        #sshy = r_[ssh, ssh[-1,:].reshape(1,-1)]
        sshy=ssh
        u = diff(sshy,axis=0)/diff(y,axis=0)
        u = (u[0:-1,:]+u[1:,:])/2
        u=r_[u[0:1,:],u,u[-1,:].reshape(1,-1)]
        #u[-2,:]=u[-1,:]
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
def ssh2psi(lon,lat,ssh):
    from numpy import pi,sin,ndim,ma
    from pylab import meshgrid
    g = 9.8; #r = 6371.e3
    omega = 0.729e-4
    if ndim(lon)==1:
        lon,lat = meshgrid(lon,lat)
    f=2.0*omega*sin(lat*pi/180.0)
    psi = g/f * ssh
    if ma.is_masked(ssh):
        psi=ma.array(psi,mask=ssh.mask)
    return psi

def ssh2uv(lon,lat,ssh):
    """ calculate geostrophic velocity from SSH data
    lon[Nx], lat[ny]
    ssh[lat,lon] """
    from numpy import pi,sin,ndim
    from pylab import meshgrid
    g = 9.8; #r = 6371.e3
    omega = 0.729e-4
    if ndim(lon)==1:
        lon,lat = meshgrid(lon,lat)
    f=2.0*omega*sin(lat*pi/180.0)
    psi = g/f * ssh
    u, v = gradxy(lon,lat,psi)
    u = -1*u
    return u,v

def findsubdomain(lon,lat,p):
    """find the index number for the subdomain
     p=[lon0,lon1,lat0,lat1] in lon[nx] and lat[ny] """
    from numpy import argmin,array
    lon = array(lon)
    lat = array(lat)
    if lon.min()<0 and array(p[:2]).min()>0:
        lon[lon<0]=lon[lon<0]+360
    if lon.min()>=0 and array(p[:2]).min()<0:
        if p[0]<0:
            p[0]=p[0]+360
        if p[1]<0:
            p[1]=p[1]+360
            
    lon0, lon1, lat0, lat1 = p
    i0,i1,j0,j1 = argmin(abs(lon-lon0)), argmin(abs(lon-lon1)), argmin(abs(lat-lat0)), argmin(abs(lat-lat1))
    return [i0,i1+1,j0,j1+1]

def subtractsubdomain(lon,lat,p,data,index=False):
    """find the index number for the subdomain
     p=[lon0,lon1,lat0,lat1] in lon[nx] and lat[ny] """
    import numpy as np
    if index:
        #ic, jc, di, dj = p
        ic,i1,jc,j1=p 
        di=i1-ic
        dj=j1-jc
        data = np.roll(np.roll(data, shift=-ic+di+1,axis=1),shift=-jc+dj+1,axis=0)
        lon = np.roll(lon, shift=-ic+di+1)[:2*di]
        lat = np.roll(lat, shift=-jc+dj+1)[:2*dj]
        d = data[:2*dj,:2*di]
    else:
        i0,i1,j0,j1 = findsubdomain(lon,lat,p)
        data = np.roll(np.roll(data, shift=-i0+1,axis=1),shift=-j0+1,axis=0)
        lon = np.roll(lon, shift=-i0+1)[:i1-i0]
        lat = np.roll(lat, shift=-j0+1)[:j1-j0]
        d = data[:j1-j0,:i1-i0]
        
    return lon, lat, d

def N2rho(rhos,n2,dz,rho0=1035., g = 9.81, direction='top2bottom'):
    """
    Calculate the mean density field from a surface density field and a stratification profile,
    :math:`\rho = - {\rho_0 over g} \int_z N^2 dz`
    """
    import numpy as np
    nz = n2.size; ny,nx = rhos.shape
    rho = np.zeros((nz+1,ny,nx))
    rho[0,:,:] = rhos
    dz = abs(dz)
    for i in range(1,nz+1):
        if direction=='top2bottom':
            rho[i,:,:] = rho[i-1,:,:] + dz[i-1]*n2[i-1]*rho0/g
        else:
            rho[i,:,:] = rho[i-1,:,:] - dz[i-1]*n2[i-1]*rho0/g
    return rho


def efromp(p,dx,dy):
    """
    calculate energy from streamfunction field
    """
    from numpy import diff
    if len(p.shape)==3:
        u=diff(p,axis=1)/dy
        v=diff(p,axis=2)/dx
    elif p.shape.size==2:
        u=diff(p,axis=0)/dy
        v=diff(p,axis=1)/dx
    e=efromuv(u,v)
    return e


def efromuv(u,v):
    "calculate the kinetic energy from u,v field."
    if len(u.shape)==3:
        vv=(v[:,1:,:]+v[:,:-1,:])/2.
        uu=(u[:,:,1:]+u[:,:,:-1])/2.
    elif u.shape.size==2:
        vv=(v[1:,:]+v[:-1,:])/2.
        uu=(u[:,1:]+u[:,:-1])/2.
    return (uu**2+vv**2)/2.


def bispec2d(yss,m1,m2,m3,window=1,norm=1):
    """
    calculate the 2d bicoherence
    specify m1,m2,m3 so that m1+m2=m3 (0...Nt/2)
    ys[nr , nm ,nt]
    b(k,l)
    return frequency 0:0.5 and b(0<k<0.5,0<l<0.5)
    """
    from scipy.linalg.basic import hankel
    from scipy.fftpack import fft2,fftfreq
    from numpy import zeros,ones,sin,linspace,pi,arange,r_
    ys=yss.copy()
    nr,ny,nt=ys.shape
    #windowing
    wind=sin(linspace(0,pi,ny)).reshape(ny,1)**2*sin(linspace(0,pi,nt)).reshape(1,nt)**2
    for i in arange(nr):
        print "detrending"
        y=ys[i,:,:]
        #2d detrend
        y=(y-y.mean(axis=1).reshape(ny,1)*ones((1,nt)))
        y=(y-ones((ny,1))*y.mean(axis=0).reshape(1,nt))
        if window ==1:
            ys[i,:,:]=y*wind
        else:
            ys[i,:,:]=y
    print "down detrending"
    #2D fft"
    print "calculating fft2d"
    spec2d=fft2(ys,axes=(-2,-1))/nt/ny
    print "done fft2d"
    mask=hankel(arange(nt),r_[nt-1,arange(nt-1)])
    bi=zeros((nr,nt,nt))
    pkl=zeros((nr,nt,nt))
    pm=zeros((nr,nt,nt))
    for i in arange(nr):
        print 'nr=',i
        sm1=spec2d[i,m1,:]
        sm2=spec2d[i,m2,:]
        sm3=spec2d[i,m3,:]
        pkl[i,:,:]=abs(sm1.reshape(nt,1)*(sm2.conj().reshape(1,nt)))**2
        pm[i,:,:]=(abs(sm3)**2)[mask]
        bi[i,:,:]=(sm1.reshape(nt,1) * sm2.reshape(1,nt)) * (sm3.conj()[mask])
    if norm==1:
        bi=abs(bi.mean(axis=0))**2/pkl.mean(axis=0)/pm.mean(axis=0)
    else:
        bi=abs(bi.mean(axis=0))**2
    print "max of b(k,l)=",bi.max()
    freq=fftfreq(nt)
    return freq,bi

def bispec(ys,window=1,norm=1):
    from scipy.linalg.basic import hankel
    from scipy.fftpack import fft,fftfreq
    from numpy import zeros,ones,sin,linspace,pi,arange,r_
    """
    calculate the bicoherence
    ys[number of realization, sample size]
    b(k,l)
    return frequency 0:0.5 and b(0<k<0.5,0<l<0.5)
    """
    nr,nt=ys.shape
    ys=(ys-ys.mean(axis=1).reshape(nr,1)*ones((1,nt)))#detrend
    if window ==1:
        wind=ones((nr,1))*sin(linspace(0,pi,nt)).reshape(1,nt)**2
        ys=ys*wind
    #nfft=2**floor(log2(nt))
    nfft=nt
    mask=hankel(arange(nfft),r_[nfft-1,arange(nfft-1)])
    bi=zeros((nr,nfft,nfft))
    pkl=zeros((nr,nfft,nfft))
    pm=zeros((nr,nfft,nfft))
    for i in arange(nr):
        y=ys[i,:]
        spec=fft(y)/nfft
        pkl[i,:,:]=abs(spec.reshape(nfft,1)*(spec.conj().reshape(1,nfft)))**2
        pm[i,:,:]=(abs(spec)**2)[mask]
        bi[i,:,:]=(spec.reshape(nfft,1) * spec.reshape(1,nfft)) * (spec.conj()[mask])
    if norm==1:
        bi=abs(bi.mean(axis=0))**2/pkl.mean(axis=0)/pm.mean(axis=0)
    else:
        bi=abs(bi.mean(axis=0))**2

    print "max of b(k,l)=",bi.max()
    freq=fftfreq(nfft)
    return freq,bi


def checkmake(fn):
    """
    Make a folder with name 'fn'.
    """
    import os
    if not os.path.exists(fn):
        os.mkdir(fn)
    return

def bindata(lons, lats, points, gridsize, distance=2.,p=[0,360,-74,-24]):
    """ bin random data points to grids
    loc_points: 2 x dpoints, lon, lat
    loc_grids: 2 x dgrids, x.ravel, y.ravel
    """
    from scipy import spatial
    from numpy import arange, meshgrid,sqrt, exp, array,c_
    locs = c_[lons.ravel(), lats.ravel()]
    tree = spatial.cKDTree(locs)
    x0,x1,y0,y1=p
    x,y = meshgrid(arange(x0,x1,gridsize), arange(y0, y1, gridsize))

    grids = zip(x.ravel(), y.ravel())

    index = tree.query_ball_point(grids, distance)

    Tmis=[]
    sample_size=[]
    lon1d,lat1d = lons.ravel(), lats.ravel()
    x1d, y1d = x.ravel(), y.ravel()
    for i in range(x.size):
        ip = index[i]
        if len(ip) == 0:
            Tmis.append(0)
            sample_size.append(0)
        else:
            dis = ((lon1d[ip]-x1d[i])**2+(lat1d[ip]-y1d[i])**2)
            weight = exp(-(dis/distance**2/4.))
            weight = weight/weight.sum()
            Tmis.append((weight*tcos.ravel()[ip]).sum())
            sample_size.append(len(ip))
            
    a = ma.masked_equal(array(Tmis).reshape(x.shape[0], x.shape[1]),0)
    sample_size = ma.masked_equal(array(sample_size).reshape(x.shape[0], x.shape
[1]),0)
    return x,y,a,sample_size

def find_overlap(xlong,xshort):
    """ find overlap between xlong and xshort,
    suppose xshort is a subset of xlong
    """
    ip = ~(xlong<xshort.min()) + (xlong>xshort.max())
    
    return ip, xlong[ip]

def kdtree_sample2d(xin, yin, z2d, xout, yout, distance=2.,method='linear'):
    """ bin random data points to grids
    loc_points: 2 x dpoints, lon, lat
    loc_grids: 2 x dgrids, x.ravel, y.ravel
    """
    from scipy import spatial, ma
    from numpy import meshgrid, exp, array,c_, where
    xip,xin  = find_overlap(xin, xout)
    yip,yin = find_overlap(yin, yout)
    z2ds = z2d[where(yip==True)[0],:][:,where(xip==True)[0]]
        
    ismask = ma.is_masked(z2ds)
    
    xin2d,yin2d = meshgrid(xin, yin)
    
    if ismask:
        xin1d, yin1d = xin2d[z2ds.mask==False].ravel(), yin2d[z2ds.mask==False].ravel()
        z1d = z2ds[z2ds.mask==False]
        locs = c_[xin1d, yin1d]
    else:
        xin1d, yin1d = xin2d.ravel(), yin2d.ravel()
        z1d = z2ds.ravel()
        locs = c_[xin1d, yin1d]
    
    tree = spatial.cKDTree(locs)
    
    grids = zip(xout, yout)

    index = tree.query_ball_point(grids, distance)

    Tmis=[]
    sample_size=[]
    
    for i in range(xout.size):
        ip = index[i]
        if len(ip) == 0:
            Tmis.append(999999)
            sample_size.append(0)
        else:
            dis = ((xin1d[ip]-xout[i])**2+(yin1d[ip]-yout[i])**2)
            if method=='linear':
                dis = ma.masked_greater(dis**0.5, distance)
                weight = distance - dis**0.5
            else:
                weight = exp(-(dis/distance**2))
            weight = weight/weight.sum()
            Tmis.append((weight*z1d[ip]).sum())
            sample_size.append(len(ip))
            
    zout = ma.masked_greater(array(Tmis),1e5)
    sample_size = ma.masked_equal(array(sample_size),0)
    return zout, sample_size

def line_sample2d(x,y,z,x1,y1):
    """sample z along a path [x,y]
    x,y in pixel coordinates"""
    from scipy.interpolate import RectBivariateSpline as rbs
    # Extract the values along the line, using cubic interpolation
    f = rbs(x,y,z.T)
    return f.ev(x1,y1)
    #return scipy.ndimage.map_coordinates(z, np.vstack((y,x)))

def D2(z, n2, f0):
    from numpy import diff, matrix, r_, diag, identity
    nz = len(n2)
        #n2 = r_[1, n2, 1] # Nz->Nz+2 N^2 points including buoyancy at the surface and bottom
    zf = r_[z[0]-(z[1]-z[0]), z, z[-1]+(z[-1]-z[-2])] #Nz->Nz+2
        
        #n2 = r_[n2[0], n2]
    def mp(a,b):
        return matrix(a)*matrix(b)
    
    zc = r_[zf[0]-0.5*(zf[1]-zf[0]), twopave(zf), zf[-1]+0.5*(zf[-1]-zf[-2])] #Nz+3 \psi points
    dzc = diff(zc) # d\psi /dz (Nz+2)
    dzf = diff(zf) # (Nz+1)
    
    #nz = n2.size
    nz = nz
    F = diag(r_[0, f0**2 / n2, 0])
    d1 = mp(diag(1. / dzc) , diff(-identity(nz+2),axis=1))
    d2 = mp(diag(1./ dzf) , diff(identity(nz+2),axis=0))
    
    eigD2 = mp(d2 , mp(F,d1))
    
    F = diag(r_[f0**2/n2[0], f0**2 / n2, f0**2/n2[-1]])
    d2 = mp(diag(r_[1., 1./dzf, 1.]), diff(identity(nz+3),axis=1))
    d1 = mp(diag(1. / dzc), diff(-identity(nz+3),axis=0))
    sqgD2 = mp(d2,mp(F,d1)); 
    sqgD2[-1,-2:]=-1./dzc[-1], 1./dzc[-1]
    sqgD2[0,:2]=-1./dzc[0], 1./dzc[0]
    
    return eigD2, sqgD2, zf, zc

def D2check(z,n2,f0):
    import numpy as np
    from numpy import r_,diff,diag
    nz=len(z)+2
    zf = r_[z[0]-(z[1]-z[0]), z, z[-1]+(z[-1]-z[-2])] #Nz->Nz+2
    zc = r_[zf[0]-0.5*(zf[1]-zf[0]), twopave(zf), zf[-1]+0.5*(zf[-1]-zf[-2])] #Nz+3 \psi points
    dzc = diff(zc) # d\psi /dz (Nz+2)
    dzf = diff(zf) # (Nz+1)
    F = diag(r_[f0**2/n2[0], f0**2 / n2, f0**2/n2[-1]])
    D1 = np.matrix(np.diag(1. / dzc)) * np.matrix(np.diff(-np.identity(nz+1), axis=0))
    D2 = np.matrix(np.diag(np.r_[0, 1. / dzf, 0])) * \
         np.matrix(np.diff(np.identity(nz+1), axis=1))
    M1 = np.matrix(D2) * np.matrix(F) * np.matrix(D1)
    M1[0,0:2] = 1./np.r_[dzc[0],-dzc[0]]
    M1[-1,-2:] = 1./np.r_[dzc[-1], -dzc[-1]]
    return M1
    
def revert(var,axis=0):
    import numpy as np
    """
    calculate values in the grid center
    use periodic boundary condition
    """
    if len(var.shape) ==3:
        if axis == 0:
            nz, ny, nx = var.shape
            tem = np.zeros((nz + 2, ny, nx))
            tem[1:-1, :, :] = var
            tem[0, :, :] = var[0, :, :]
            tem[-1, :, :] = var[-1, :, :]
            tem = np.ma.masked_values(tem, 0)
            return (tem[1:, :, :] + tem[:-1, :, :]) / 2.
        if axis == 1:
            nz, ny, nx = var.shape
            tem = np.zeros((nz, ny + 2, nx))
            tem[:, 1:-1, :] = var
            tem[:, 0, :] = var[:, 0, :]
            tem[:, -1, :] = var[:, -1, :]
            tem = np.ma.masked_values(tem, 0)
            return (tem[:, 1:, :] + tem[:, :-1, :]) / 2.
        if axis == 2:
            nz, ny, nx = var.shape
            tem = np.zeros((nz, ny, nx + 2))
            tem[:, :, 1:-1] = var
            tem[:, :, 0] = var[:, :, 0]
            tem[:, :, -1] = var[:, :, -1]
            tem = np.ma.masked_values(tem, 0)
            return (tem[:, :, 1:] + tem[:, :, :-1]) / 2.
    if len(var.shape) == 1:
        tem = np.zeros(var.shape[0]+2)
        tem[1:-1] = var
        tem[0] = var[0]
        tem[-1] = var[-1]
        return (tem[1:]+tem[:-1])/2.
    
def Time2Str(Time):
    if isinstance(Time, str):
        return Time
    elif isinstance(Time, list):
        return '%4i%02i%02i'%(Time[0],Time[1],Time[2])
    else:
        return Time.strftime('%Y%m%d')
    
def angle(u,v):
    from numpy import dot, norm, arccos
    c = dot(u,v)/norm(u)/norm(v) # -> cosine of the angle
    return arccos(c) # if you really want the angle

def NS(lat):
    """format latitude symbol
    NS(80): N80
    NS(-80): S80"""
    if lat>=0:
        return 'N%02i'%lat
    else:
        return 'S%02i'%(abs(lat))
    
    
    

def imirrormirror(v):
    """
    reverse process of mirromirror
    output sst(ny,nx)
    input sst(2ny-1,2nx-1)
    """
    from numpy import r_,c_,zeros
    def im2d(v):
        v = r_[v[-1,:].reshape(1,-1),v,v[0,:].reshape(1,-1)]
        v = c_[v[:,-1].reshape(-1,1),v,v[:,0].reshape(-1,1)]
        v = fourpointavg(v)
        ny,nx = v.shape
        return v[:ny/2,:nx/2]
    if len(v.shape) == 2:
        return im2d(v)
    elif len(v.shape) == 3:
        nz,ny,nx = v.shape
        vv = zeros((nz,(ny+1)/2,(nx+1)/2))
        for i in range(nz):
            vv[i,:,:] = im2d(v[i,:,:])
        return vv

def mirrormirror(sst):
    """
    input sst(ny,nx)
    output sst(2ny-1,2nx-1)
    """
    from numpy import r_,c_
    if len(sst.shape) == 1:
        sst = r_[sst.flatten(),sst[-1::-1].flatten()]
        return twopave(sst)
    else:
        sst =  r_[sst, sst[-1::-1,:]]
        sst =  c_[sst, sst[:,-1::-1]]
        return  fourpointavg(sst)

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str
    
