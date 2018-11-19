from astropy import units as u
import numpy as np

def skyscale(distance=None,distanceunit=u.parsec,angle=None,angleunit=u.arcsec,size=None,sizeunit=u.au):
    """
    
    
    Allows easy calculation of astronomical scales using a variety of units.
    Input two of: distance, angle and size, and optionally include your preferred units for both input and output.
    
    
    -----------------------------------------------------------------------
    
    
    Parameters: 
    
    distance: distance of the object from the Earth.
    
    angle: angular size of the object.
    
    size: physical size of the object.
       
    distanceunit: the astropy units you wish to specify for distance, e.g. u.parsec for parsecs.
                  Defaults to parsecs if left unspecified.

    angleunit: the astropy units you wish to specify for angle, e.g. u.arcsec for arcseconds.
               Defaults to arcseconds if left unspecified.
               
    sizeunit: the astropy units you wish to specify for size, e.g. u.au for astronomical units.
               Defaults to astronomical units if left unspecified.
    
    
    -----------------------------------------------------------------------
    
    """
    if distance and angle and size:
        raise ValueError('Please specify only two of distance, angle and size.')    
        
    if distance and angle:
        
        distance = distance*distanceunit
        angle = angle*angleunit
        return (distance.to(u.m)*np.tan(angle.to(u.degree))).to(sizeunit)
    
    if distance and size:
        
        distance = distance*distanceunit
        size = size*sizeunit
        return (np.arctan(size.to(u.m)/distance.to(u.m))).to(angleunit)
    
    if angle and size:
        
        angle = angle*angleunit
        size = size*sizeunit
        return (size.to(u.m)/np.tan(angle.to(u.degree))).to(distanceunit)
    
    else:
        raise ValueError('Please specify only two of distance, angle and size.')

        
    