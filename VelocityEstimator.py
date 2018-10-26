import numpy as np

def est(distance,radius,mass):
    """
    
    Returns maximum radial velocity (i.e. keplerian velocity) of gas in a disk.
    
    -----------------------------------------------------------------------
    
    Parameters: 
    
    Distance of the star from the Earth in parsec.
    
    Radius of the disk in arcseconds.
    
    Mass of the star in solar masses.
    
    -----------------------------------------------------------------------
    
    """
    R = dist*radius*1.5e11 #calculate raidus in m
    M = mass*2e30#calculate mass in kg
    G = 6.674e-11   #gravitational constant
    V=np.sqrt(G*M/R)#calculate velocity
    
    return V/1000 #V in km/s
