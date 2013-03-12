'''
A collection of the basic parameters for self sciences

Created on Jan 24, 2013
@author: Jinbo Wang 
Email: jinbow@gmail.com
'''

class planet:
    ''' A collection of the constants for planets
    '''
    def __init__(self):
        from math import pi
        self.radius_mean = 0 # meter 
        self.radius_equatorial = 0 # meter 
        self.radius_polar = 0 # meter 
        self.surface = 0 # meter^2
        self.surface_land = 0 # m^2 
        self.surface_water = 0 #^2
        self.volume = 1.08321e12 # m^3
        self.mass = 5.9736e24 # kg 
        self.density_mean = 5.515  #kg/m^3 
        self.rotation_period = 23.*3600 + 56.*60 + 4.1 # second
        self.omega = 2*pi/(23.*3600 + 56.*60 + 4.1) # rad / second
        self.g = 0 # m/s^2
        return
    def f(self,latitude = 45.):
        import math
        return 2.*self.omega*math.sin(latitude)
        
class earth():
    def __init__(self):
        import math        
        """ Earth """
        self.radius_mean = 6371.0e3
        self.radius_equatorial = 6378.1e3
        self.radius_polar = 6356.8e3
        self.surface = 510072000.e6
        self.surface_land = 148940000e6
        self.surface_water = 361132000e6
        self.volume = 1.08321e21
        self.mass = 5.9736e24
        self.density_mean = 5.515e3
        self.rotation_period = 23.*3600 + 56.*60 + 4.1
        self.omega = 2*math.pi/(23.*3600 + 56.*60 + 4.1)
        self.g = 9.780327
        return
    