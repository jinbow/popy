'''
several handy subroutines related to IO
Created on Mar 7, 2013
@author: Jinbo Wang <jinbow@gmail.com>
Scripps Institution of Oceanography
'''


class array:
    def __init__(self,lon=[0],lat=[0],depth=[0],time=[0]):
        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.time = time
        self.data = {}
    def tonetcdf(self,varname='d',fn = ''):
        """ save numpy array to netcdf file """
        import sys, os, time
        from netCDF4 import Dataset
        
        if fn =='':
            fn = sys.argv[0].split('.')[0]+'_tmp.nc'
        if os.path.exists(fn):
            os.remove(fn)
        
        rootgrp = Dataset(fn, 'w', format='NETCDF4')
        rootgrp.history = 'Created by '+sys.argv[0] + time.ctime(time.time())
        
        dims = self.data.shape
        if len(dims) == 2:
            self.data = self.data.reshape(1,1,dims[0],dims[1])
        if len(dims) == 3:
            self.data = self.data.reshape(1,dims[-3],dims[-2],dims[-1])
            
        depth = rootgrp.createDimension('depth', None)
        time = rootgrp.createDimension('time', None)
        lat = rootgrp.createDimension('lat', dims[-2])
        lon = rootgrp.createDimension('lon', dims[-1])
        
        times = rootgrp.createVariable('time','f8',('time',))
        depths = rootgrp.createVariable('depth','f4',('depth',))
        latitudes = rootgrp.createVariable('latitude','f4',('lat',))
        longitudes = rootgrp.createVariable('longitude','f4',('lon',))
        
        latitudes.units = 'degrees north'
        longitudes.units = 'degrees east'
        times.units = 'hours since 0001-01-01 00:00:00.0'
        times.calendar = 'gregorian'

        if len(self.time)==1:
            self.time = range(dims[0])
        if len(self.depth) == 1:
            self.depth = range(dims[1])
        if len(self.lat) == 1:
                self.lat = range(dims[2])
        if len(self.lon) == 1:
                self.lon = range(dims[3])
                
        depths[:] = self.depth
        times[:] = self.time
        latitudes[:] = self.lat
        longitudes[:] = self.lon
            
        temp = rootgrp.createVariable(varname,'f4',('time','depth','lat','lon',))
        temp[:] = self.data
        rootgrp.close()

def saveh5(filename, datasetname, data,driver=None,):
    import h5py
    
    hfile = h5py.File(filename,driver=driver)
    print "save data to HDF5, data size=",datasetname,data.shape
    dst = hfile.require_dataset(datasetname, shape=data.shape, dtype=data.dtype)
    dst[:] = data[:]
    
    hfile.close()
    
    return
    
def save(fn, d):
    """save data to file with format determined by the filename extension"""
    
    ext = fn.split('.')[-1]
    if ext == '.mat':
        from scipy.io import savemat
        #f=open(fn,'wb')
        savemat(fn,d)
        #f.close()
        print 'save data into '+fn
    
    return
    
def load(fn, varn=''):
    """load data from file with format determined by the filename extension"""
    ext = fn.split('.')[-1]
    if ext == '.mat':
        from scipy.io import loadmat
        f=loadmat(fn,squeeze_me=True)
        print 'load data from '+fn
        if varn=='':
            return f 
        else:
            return f[varn]
    return

def remove(fn):
    import os
    if os.path.exists(fn):
        os.remove(fn)
        return

def savevtk(v, fn = '', lons=[], lats=[], levels=[]):
    """save tracer field into vtk format,
    """
    import pyvtk  
    def norm(c):
        return (c-c[0])/(c[-1]-c[0])
    z = levels/abs(levels).max() * (lons.max()-lons.min()) * 0.7
    
    point_data = pyvtk.PointData(pyvtk.Scalars(v.T.flatten()))
    vtk_object = pyvtk.VtkData(pyvtk.RectilinearGrid(x=lons,y=lats,z=z), point_data)
    vtk_object.tofile(fn)
    return
    
def loadnc(fn, var='All'):
    from netCDF4 import Dataset
    
    f = Dataset(fn, 'r').variables
    d={}
    if var=='All':
        varn=f.keys()
        
        for var in varn:
            d[var] = f[var][:].squeeze()
        return d
    
    else:
        d=f[var][:].squeeze()

        return d


def savenetcdf4(v, fn = '', ndims=4, lons=[], lats=[], levels=[], records=[]):
        """ save numpy array to netcdf file 
        v={vname:vvalue, vname:vvalue}
        """
        import sys, os
        import time as mod_time
        from netCDF4 import Dataset
        
        if fn =='':
            fn = sys.argv[0].split('.')[0]+'_tmp.nc'
        if os.path.exists(fn):
            os.remove(fn)
        
        rootgrp = Dataset(fn, 'w', format='NETCDF4')
        rootgrp.history = 'Created by '+sys.argv[0] + mod_time.ctime(mod_time.time())

        if len(levels) !=0:
            depth = rootgrp.createDimension('depth', len(levels))
            depths = rootgrp.createVariable('depth','f8',('depth',))
            depths[:] = levels
        if len(records)!=0:
            time = rootgrp.createDimension('time', None)
            times = rootgrp.createVariable('time','f8',('time',))
            times.units = 'hours since 0001-01-01 00:00:00.0'
            times.calendar = 'gregorian'
            times[:] = records
        rootgrp.createDimension('lat', len(lats))
        latitude = rootgrp.createVariable('lat','f8',('lat',))
        latitude.units = 'degrees north'
        latitude[:] = lats[:]
        
        rootgrp.createDimension('lon', len(lons))
        longitude = rootgrp.createVariable('lon','f8',('lon',))
        longitude.units = 'degrees east'
        longitude[:] = lons[:]
        
        
        for varname in v.keys():
            if ndims == 4:
                rootgrp.createVariable(varname,'f8',('time','depth','lat','lon',), fill_value=99999)[:] = v[varname][:]
            elif ndims == 3:
                rootgrp.createVariable(varname,'f8',('depth','lat','lon',), fill_value=99999)[:] = v[varname][:]
            else:
                rootgrp.createVariable(varname,'f8',('lat','lon',), fill_value=99999)[:] = v[varname][:]
        rootgrp.close()
        
def savenc(var,varname='d',fn='',lon=[],lat=[],dep=[],time=[]):
    import numpy as np
    from scipy.io import netcdf
    import sys
    from datetime import date
    
    if fn =='':
        fn = varname+'.nc'
    #detect dimensions"
    dims = [1,1,1,1]
    for i in range(var.ndim):
        n = var.ndim - i
        dims[-i-1] = var.shape[-i-1]
    #create netcdf file"
    f = netcdf.netcdf_file(fn, 'w')
    f.history = "created by "+sys.argv[0] +" on "+date.today().strftime("%d/%m/%y")
    f.Filled_Value = -9999
    f.createDimension('time',dims[0])
    f.createDimension('depth',dims[1])
    f.createDimension('latitude',dims[2])
    f.createDimension('longitude',dims[3])
    longitude = f.createVariable('longitude','>f8',('longitude',))
    latitude = f.createVariable('latitude','>f8',('latitude',))
    depth = f.createVariable('depth','>f8',('depth',))
    tim = f.createVariable('time','>f8',('time',))
    varn = f.createVariable(varname, '>f8', ('time','depth','latitude','longitude'))
    varn[:] = var[:]
    if len(lon) == 0 or len(lon)!= dims[3]:
        print "use evenly spaced integer for longitude"
        lon = np.arange(dims[3])
    if len(lat) == 0 or len(lat) != dims[2]:
        print "use evenly spaced integer for latitude"
        lat = np.arange(dims[2])
    if len(dep) == 0 or len(dep) != dims[1]:
        print 'len(dep)=',len(dep),'dims[1]=',dims[1], "use evenly spaced integer for depth"
        dep = np.arange(dims[1])
    if len(time) == 0 or len(time) != dims[0]:
        print "use evenly spaced integer for time"
        time = np.arange(dims[0])
    longitude[:] = lon[:]
    latitude[:] = lat[:]
    depth[:] = dep
    tim[:] = time
    f.flush()
    print "saved data to " + fn
    return

def savedata(v, varname, step, fmt = 'mitgcm'):
    """ write numpy array to binary data following the rule of MITgcm output.
    the first three dimensions of array (v) should always be (..., z,y,x)
    The fourth dimension is referred to as record.
    """ 
    import numpy as np
    if v.ndim == 2:
        v = np.expand_dims(np.expand_dims(v,axis=0),axis=0)
    if v.ndim == 3:
        v = np.expand_dims(v,axis=0)
    dims = v.shape
    
    if fmt == 'mitgcm':
        
        v.astype('>f4').tofile(varname+'%10i.data'%step)
        f = open(varname+'%10i.meta'%step,'w')
        f.writelines('nDims = [   %3i ];\n'%len(dims))
        f.writelines('dimList = [\n')
        f.writelines('%6i, 1, %6i, \n'%(dims[-1],dims[-1]))
        f.writelines('%6i, 1, %6i, \n'%(dims[-2],dims[-2]))
        if len(dims) == 3:
            f.writelines('%6i, 1, %6i, \n];'%(dims[-3],dims[-3]))
        f.writelines('dataprec = [ ''float32'' ];\n')
        f.writelines('nrecords = [      %4i ];\n'%dims[0])
   
        f.writelines('timeStepNumber = [    0 ];\n')
        
         
    return
