'''
Mapping tools built on Basemap
@author: Jinbo Wang <jinbow@gmail.com>
@Institution: Scripps Institution of Oceanography
Created on Mar 12, 2013
'''

def setmap(lon=[],lat=[],p=[0,360,-90,90],
           proj='cyl',fix_aspect=False,
           fillcontinents=True,continentcolor='k',
           dlat=10,dlon=30):
    """ set up a map and return a handle.
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
#     if crop !=[]:
#         p = crop
#     elif len(lon) !=0:
#         p[0],p[1],p[2],p[3]=lon.min(),lon.max(),lat.min(),lat.max()
#         
    def mint(p):
        return 5*int(p[1]-p[0])/6/5

#             
    m = Basemap(projection=proj,llcrnrlon=p[0],urcrnrlon=p[1],
                llcrnrlat=p[2],urcrnrlat=p[3],resolution='l',fix_aspect=fix_aspect)
    m.drawparallels(np.arange(-90.,90.,dlat),labels=[1,0,0,1],linewidth=0)
    m.drawmeridians(np.arange(0,361.,dlon),labels=[1,0,0,1],linewidth=0)
    m.drawcoastlines()
    if fillcontinents or continentcolor==None:
        m.fillcontinents(color=continentcolor)
    m.drawmapboundary(fill_color='w')
    return m

def mapwithtopo(p,ax=[],cutdepth=[],aspect=2.5,cmap=[],dlon=30,dlat=10,smooth=False):
    import popy,os
    from netCDF4 import Dataset
    from mpl_toolkits.basemap import shiftgrid
    from matplotlib.colors import LightSource
    import pylab as plt
    import numpy as np
    
    etopofn='/net/mazdata2/jinbo/mdata5-jinbo/obs/ETOPO/ETOPO1_Bed_g_gmt4.grd'
    etopo = Dataset(etopofn,'r').variables
    x,y,z=etopo['x'][1:],etopo['y'][1:],etopo['z'][1:,1:]
    dx,dy=5,5
    x=x.reshape(-1,dx).mean(axis=-1)
    y=y.reshape(-1,dx).mean(axis=-1)
    z=z.reshape(y.size,dx,x.size,dx).mean(axis=-1).mean(axis=1)
    if smooth:
        z=popy.utils.smooth2d(z,window_len=3)
    
    if cutdepth!=[]:
        z[z<cutdepth]=cutdepth
        
    if ax==[]:
        fig=plt.figure()
        ax=fig.add_subplot()
    if cmap==[]:
        cmap=plt.cm.Greys
        
    z,x = shiftgrid(p[0],z,x,start=True)
    lon,lat,z = popy.utils.subtractsubdomain(x,y,p,z)
    
    m = setmap(p=p,dlon=dlon,dlat=dlat)
    x, y = m(*np.meshgrid(lon, lat))
    ls = LightSource(azdeg=90, altdeg=45)
    rgb = ls.shade(z, cmap=cmap)
    m.imshow(rgb, aspect=aspect)
    
    plt.savefig('/tmp/tmp.png')
    os.popen('eog /tmp/tmp.png')
    return etopo

def polygon(m,coor,facecolor='gray',alpha=0.4):
    """draw a patch on a map m using coordinate in coor
    coor is a [n,2] array, [:,0] contains longitude and [:,1] latitude.
    coor also can be [lon0,lon1,lat0,lat1] for a rectangle box
    popy.map.polygon(m,coor)
    """
    from numpy import array, c_
    from matplotlib.patches import Polygon
    import matplotlib.pylab as plt
    coor = array(coor) # in case array is a list
    if coor.ndim ==1 and coor.size ==4:
        p=coor
        coor = array([[p[0],p[2]],
                      [p[1],p[2]],
                      [p[1],p[3]],
                      [p[0],p[3]] ])
    print coor
    x,y = m(coor[:,0],coor[:,1])
    poly=Polygon(c_[x,y],facecolor=facecolor,alpha=alpha,closed=True)
    plt.gca().add_patch(poly)
    return
