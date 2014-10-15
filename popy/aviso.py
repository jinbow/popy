'''
Created on Sep 3, 2013

@author: Jinbo Wang
@contact: <jinbow@gmail.com>
@organization: Scripps Institution of Oceanography

'''

class aviso:
    import pickle
    from numpy import loadtxt
    from datetime import datetime
    
    def __init__(self):
        import numpy as np
        self.root = '/net/mazdata2/jinbo/mdata5-jinbo/obs/AVISO/ftp.aviso.oceanobs.com/global/dt/ref/madt/merged/'
        f1 = open(self.root + 'filelist_ssh.data')
        f2 = open(self.root + 'filelist_uv.data')
        self.FileList = {'h':[line.rstrip() for line in f1],
                         'uv':[line.rstrip() for line in f2]}
        self.lat = np.arange(721) * 0.25 - 90.
        self.lon = np.arange(1440) * 0.25
        self.date = self.LoadDate()
        f1.close()
        f2.close()
        return
    

    def LoadDate(self):
        sshdays = []
        for fn in self.FileList['h']:
            sshdays.append(DateFromFile(fn))
        return sshdays
    
    def FileName(self, YMD, vname='h'):
        """timestampe in format :YYYYMMDD"""
        import popy
        flist = self.FileList[vname]
        if YMD=='mean':
            return self.root + vname + '_mean.pkl'
        YMD = popy.utils.Time2Str(YMD)
        for line in flist:
            if YMD in line[:-11]:
                return self.root + line
            
    def loadgrid(self):
        import os, pickle
        from netCDF4 import Dataset
        fnpkl = self.root + 'grid.pkl'
        if os.path.exists(fnpkl):
            d = pickle.load(open(fnpkl, 'r'))
            return d
        else:
            fn = self.root + self.FileList['ssh'][10].strip()
            f = Dataset(fn, 'r').variables
            lat, lon = f['NbLatitudes'][:], f['NbLongitudes'][:]
            date = self.timestamps()
            g = {'lat':lat, 'lon':lon, 'date':date}
            pickle.dump(g, open(fnpkl, 'wb'))
            return g
        
    def LoadData(self, YMD='19921223', subdomain=[], index=True):
        from numpy import ma,array
        import popy
        from netCDF4 import Dataset
        import pickle,os,glob
        
        def load(fn):
            ff = Dataset(fn, 'r')
            f = ff.variables
            ssh = f['Grid_0001'][:].T
            ssh = ma.masked_greater(ssh, 500)
            ssh = (ssh[:-1, :] + ssh[1:, :]) / 2.
            ff.close()
            return ssh
        
        fn = self.FileName(YMD)
        if YMD=='allinone':
            fno = self.root+'h_allinone1993-2012.pkl'
            if os.path.exists(fno):
                da = ma.masked_greater(pickle.load(open(fno,'r')),1e10)
                return da
            else:
                fnlist=sorted(glob.glob(self.root+'h/dt_ref_global_merged_madt_h_qd_*.nc'))[12:]
                ens = []
                for fn in fnlist:
                    print fn
                    ens.append(load(fn))
                ens = array(ens)
                pickle.dump(ens, open(fno,'wb'))
                return ens
            
        elif YMD=='mean':
            fn = self.root+'ssh_mean1993-2012.pkl'
            
            if os.path.exists(fn):
                
                ssh=pickle.load(open(fn,'r'))
                print fn
            else:
                mean=0
                for fns in self.FileList['h']:
                    print 'read '+fns
                    mean+=load(self.root+fns)
                ssh = mean/len(self.FileList['h'])
                pickle.dump(ssh, open(fn,'wb'))
        elif YMD=='eke':
            fno = self.root+'h_allinone1993-2012.pkl'
            da = ma.masked_greater(pickle.load(open(fno,'r')),1e10)
            da=da[782:938,:,:]
            
            fn = self.root+'ssh_eke1993-2012.pkl'
            if os.path.exists(fn):
                ssh=pickle.load(open(fn,'r'))
            else:
                mean=0
                for fns in self.FileList['h']:
                    print 'read '+fns
                    mean+=load(self.root+fns)
                ssh = mean/len(self.FileList['h'])
                pickle.dump(ssh, open(fn,'wb'))
        else:
            ssh = load(fn)
            self.lat = (self.lat[1:] + self.lat[:-1]) / 2.  # shift to oisst2 data grid
        
        lon, lat = self.lon, self.lat
        
        if subdomain != []:
            lon, lat, ssh = popy.utils.subtractsubdomain(self.lon, self.lat, subdomain, ssh, index=index)
        
        return lon, lat, ssh
    
def DateFromFile(fn):
    """subtract time information from a filename"""
    import datetime
    i = 33
    return datetime.date(int(fn[i:i + 4]), int(fn[i + 4:i + 6]), int(fn[i + 6:i + 8]))
        
if __name__ == '__main__':
    from pylab import contourf, savefig
    a = aviso()    
    lon, lat, ssh = a.LoadData(subdomain=[100, 100, 10, 10])
    print ssh.shape, ssh.max(), ssh.min()
    contourf(ssh.T)
    savefig('tmp.png')
    
