import numpy as np

def orbit_vel(distance,mass,radius):
    """
    
    Returns maximum radial velocity (i.e. keplerian velocity) of object around a star.
    
    -----------------------------------------------------------------------
    
    Parameters: 
    
    Distance of the star from the Earth in parsec.
    
    Mass of the star in solar masses.
    
    Radius of the orbit in arcseconds.
    
    -----------------------------------------------------------------------
    
    """
    
    R = distance*radius*1.5e11 #calculate raidus in m
    M = mass*2e30#calculate mass in kg
    G = 6.674e-11   #gravitational constant
    V=np.sqrt(G*M/R)#calculate velocity
    
    return V/1000 #V in km/s

#TODO
#allow distance or angle input units, use astropy constants.
