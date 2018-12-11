import numpy as np
from astropy import units as u
from astropy import constants as const

def orbit_vel(distance, mass, radius):
    
    """
    
    Returns maximum radial velocity (i.e. keplerian velocity) of object around another massive object.
    
    -----------------------------------------------------------------------
    
    Parameters: 
    
    distance: Distance to the objects. Units in astropy units of choice or parsecs if left unitless.
    
    mass: Mass of the central object. Units in astropy units of choice or solar masses if left unitless.
    
    radius: Radius of the orbit. Units in astropy units of choice or arcseconds if left unitless.
    
    -----------------------------------------------------------------------
    
    """
    if type(distance) != type(1*u.m):
        distance *= u.parsec
        
    if type(mass) != type(1*u.kg): 
        mass *= u.solMass
     
    if type(radius) != type(1*u.arcsec):
        radius *= u.arcsec
        
    if radius.decompose().unit == u.rad:
        radius = distance*np.tan(radius)

    V=np.sqrt(const.G*mass/radius)
    
    return V.to(u.km/u.s)
