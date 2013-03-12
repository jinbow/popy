'''
several handy subroutines related to IO
Created on Mar 7, 2013
@author: Jinbo Wang <jinbow@gmail.com>
Scripps Institution of Oceanography
'''

def savenc(var,varname='d',fn='d.nc',lon=[],lat=[],dep=[],tim=[]):
    import numpy as np
    from scipy.io import netcdf
    import sys
    from datetime import date
    
    if varname !='d':
        fn = varname+'.nc'
    #detect dimensions"
    dims = [1,1,1,1]
    for i in range(var.ndim):
        n = var.ndim - i
        dims[n] = var.shape[-i-1]
    #create netcdf file"
    f = netcdf.netcdf_file(fn, 'w')
    f.history = "created by "+sys.argv[0] +" on "+date.today().strftime("%d/%m/%y")
    f.createDimension('time',dims[0])
    f.createDimension('depth',dims[1])
    f.createDimension('latitude',dims[2])
    f.createDimension('longitude',dims[3])
    longitude = f.createVariable('longitude','>f8',('longitude',))
    latitude = f.createVariable('latitude','>f8',('latitude',))
    depth = f.createVariable('depth','>f8',('depth',))
    time = f.createVariable('time','>f8',('time',))
    varn = f.createVariable(varname, '>f8', ('time','depth','latitude','longitude'))
    varn[:] = var[:]
    if len(lon) == 0 or len(lon)!= dims[3]:
        print "use evenly spaced integer for longitude"
        lon = np.arange(dims[3])
    if len(lat) == 0 or len(lat) != dims[2]:
        print "use evenly spaced integer for latitude"
        lat = np.arange(dims[2])
    if len(dep) == 0 or len(dep) != dims[1]:
        print "use evenly spaced integer for depth"
        dep = np.arange(dims[1])
    if len(tim) == 0 or len(tim) != dims[0]:
        print "use evenly spaced integer for time"
        tim = np.arange(dims[0])
    longitude[:] = lon[:]
    latitude[:] = lat[:]
    depth[:] = dep
    time[:] = tim
    f.flush()
    print "saved data to " + fn
    return

    