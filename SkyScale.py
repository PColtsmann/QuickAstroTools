from astropy import units as u
import numpy as np

def sky_scale(distance=None,angle=None,size=None):
    
    """

    Allows easy calculation of astronomical scales using a variety of units.
    Input two of: distance, angle and size.
    Arguments can be entered as astropy Quantities to specify their units, otherwise parsecs, arcsecs and AU will be assumed respectively.
    If you wish to specify units for the output, enter just the astropy unit in the argument.
    
    -----------------------------------------------------------------------
    
    Parameters: 
    
    distance: Distance of the object from the Earth. Units in astropy units of choice or parsecs if left unitless.
    
    angle: Angular size of the object. Units in astropy units of choice or arcsecs if left unitless.
    
    size: Physical size of the object. Units in astropy units of choice or AU if left unitless.
    
    -----------------------------------------------------------------------
    
    """
    # type error will be raised if more than one argument is sleft as None    
    if size == None:
        
        size = u.AU
        
        if type(distance) != type(1*u.parsec):
            distance *= u.parsec
            
        if type(angle) != type(1*u.arcsec):
            angle *= u.arcsec
            
        size = (distance.to(u.m)*np.tan(angle.to(u.degree))).to(size)
        
        return size    
    
    elif type(size*u.m) == type(u.m*u.m): # avoids Unit != IrreducibleUnit != PrefixUnit by comparing for composite units

        if type(distance) != type(1*u.parsec):
            distance *= u.parsec
            
        if type(angle) != type(1*u.arcsec):
            angle *= u.arcsec
                        
        size = (distance.to(u.m)*np.tan(angle.to(u.degree))).to(size)# will fail if inappropriate units were entered.
        
        return size
 
    
    elif angle == None:
        
        angle = u.arcsec
        
        if type(distance) != type(1*u.parsec):
            distance *= u.parsec
            
        if type(size) != type(1*u.AU):
            size *= u.AU
            
        angle = (np.arctan(size.to(u.m)/distance.to(u.m))).to(angle)
        
        return angle
    
    elif type(angle*u.m) == type(u.m*u.m):
        
        if type(distance) != type(1*u.parsec):
            distance *= u.parsec

        if type(size) != type(1*u.AU):
            size *= u.AU
        
        angle = (np.arctan(size.to(u.m)/distance.to(u.m))).to(angle)
        
        return angle
    
    
    elif distance == None:
        
        distance = u.parsec
        
        if type(angle) != type(1*u.arcsec):
            angle *= u.arcsec
            
        if type(size) != type(1*u.AU):
            size *= u.AU
 
        distance = (size.to(u.m)/np.tan(angle.to(u.degree))).to(distance)
        
        return distance
    
    elif type(distance*u.m) == type(u.m*u.m):
        
        if type(angle) != type(1*u.arcsec):
            angle *= u.arcsec
            
        if type(size) != type(1*u.AU):
            size *= u.AU
 
        distance = (size.to(u.m)/np.tan(angle.to(u.degree))).to(distance)
        
        return distance
    
    else:
        raise ValueError('Please specify only two of distance, angle and size as astropy Quantities or floats/integers.')

        
    
