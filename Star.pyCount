import numpy as np
import SkyScale as ss
from astropy import units as u
from astropy.coordinates import SkyCoord

def GalacticCoords(Ra,Dec):
    
    """
    
    Quickly converts icrs Ra, Dec coordinates to galactic l, b coordinates.
    Assumes units of degrees if no astropy units specified.
    Returns floats l and b.
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    Ra:           Right Ascension of sky location to be converted. Units in astropy units of choice or degrees if left unitless.
    
    Dec:          Declination of sky location to be converted. Units in astropy units of choice or degrees if left unitless.
    
    --------------------------------------------------------------------------------------------------------
    
    """
    
    if type(Ra) != type(1*u.degree): #Check to see if input type == astropy.units.quantity.Quantity
        Ra *= u.degree
        
    if type(Dec) != type(1*u.degree):
        Dec *= u.degree
        
    Coords = SkyCoord(ra=Ra,dec = Dec, frame = 'icrs')
    
    return Coords.galactic.l.value,Coords.galactic.b.value

def star_count(Ra, Dec, maxdist, interval=150, Area=None, Radius=None, Side=None ):
    
    """
    
    Estimates the number of each type of star in a volume of specified direction and area/radius/side and depth.
    Uses Jurić M, et al. 2008. ApJ. 673. 864 for galactic stellar density distribution.
    Uses http://www.pas.rochester.edu/~emamajek/memo_star_dens.html for stellar type distribution.
    Most density parameters used have errors ~20%, stellar densities are lower limits and halo contribution not yet accounted for.
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    Ra: icrs Right Ascension of volume centre. Units in astropy units of choice or degrees if left unitless.
    
    Dec: icrs Declination of volume centre. Units in astropy units of choice or degrees if left unitless.
    
    maxdist: The distance up to which you wish to know the star count. Units in astropy units of choice or parsecs if left unitless.
        
    interval: Steps in numeric integration, 150 is numerically accurate to 1%.
    
    Area: Area of sky contained within volume. Units in astropy units of choice or square degrees if left unitless.
    
    Radius: The radius of the circle on the sky in which you wish to count stars up to 'maxdist'. Units in astropy size/length units of choice or degrees if left unitless.
    
    Side: The side of the square on the sky in which you wish to count stars up to 'maxdist'. Units in astropy size/length units of choice or degrees if left unitless.
    
    --------------------------------------------------------------------------------------------------------
    
    """

    
    #calculate the area of the shape on sky
    if not Area:
        if Radius != None and Side != None:
            raise ValueError('Please specify only one of Area, Radius or Side.')
            
        if Radius != None:
            if type(1.0*Radius) == float:
                Area = np.pi*Radius**2
            elif Radius.decompose().unit == u.rad:
                Area = np.pi*(Radius.to(u.degree).value)**2
            elif Radius.decompose().unit == u.m:
                Area = np.pi*(ss.sky_scale(size=Radius, distance=maxdist,angle=u.degree).value)**2

        elif Side != None:
            if type(1.0*Side) == float:
                Area = Side**2
            elif Side.decompose().unit == u.rad:
                Area = (Side.to(u.degree).value)**2
            elif Side.decompose().unit == u.m:
                Area = (ss.sky_scale(size=Side, distance=maxdist,angle=u.degree).value)**2
    else:
        if Radius != None or Side != None :
            raise ValueError('Please specify only one of Area, Radius or Side.')
        elif type(Area) == type(1*u.rad):
            Area = Area.to(u.degree*u.degree).value
        
    TypeList = ['Os','Bs','As','Fs','FDs','Gs','GDs','KDs','MDs','WDs','ESs','RGs','ALL']#D stands for Dwarf, RG stands for Red Giant
    DensList = [4.4e-8,3.2e-5,4.9e-4,0.0025,0.0024,0.0048,0.0033,0.0135,0.0917,0.0048,8.8e-4,2.7e-4,0.0984]
    Dictionary = {'Count':0}
    Dict = {}
    
    #table 10 Juric et al 2008
    #         #error
    R0 = 8000
    L = 2600  #+-20
    Lt = 3600 #+-20%
    Z0 = 25   #+-20%
    H = 300   #+-20%
    Ht = 900  #+-20%
    f = 0.12  #+-10%
    
    l, b = GalacticCoords(Ra,Dec) # convert to galactic coordinates
    
    #prep conversions into cylindrical coordinates for use in density equation
    sphtheta = b + 90
    sphphi = l
    cylrhofact = np.sin(sphtheta*np.pi/180)
    cylzfact =-1*np.cos(sphtheta*np.pi/180)
    cylphi = sphphi
    Rfact = -1*np.cos(cylphi*np.pi/180)
    
    if type(maxdist) == type(1*u.parsec):
        maxdist = maxdist.to(u.parsec).value
    
    smalldist=maxdist/interval
    areafraction = (Area/41252.96)
    totvolume = 0
    StarCount = np.zeros(len(TypeList))
    
    #numerical integration of star densities
    for i in range(1,interval+1):

        sphr = i*smalldist
        cylrho = sphr*cylrhofact
        cylz = sphr*cylzfact  
        R = cylrho*Rfact  + 8000#galactocentric
        Z = cylz
        
        dens =[]
        #calculate normalised density at interval location
        normdens = (np.exp((R0-R)/L)*np.exp(-(Z+Z0)/H) +f*np.exp((R0-R)/Lt)*np.exp(-(Z+Z0)/Ht))

        #calculate individual densities, the interval volume and star counts.
        for i in range(0, len(TypeList)):  
            dens.append(DensList[i]*normdens)
        area = areafraction*4*np.pi*sphr**2
        volume = area*smalldist
        totvolume += volume
        StarCount += np.array(dens)*volume
    
    #create output dictionary
    for i in range(0,len(TypeList)):
        Dictionary['Density'] = DensList[i]
        Dictionary['Count'] = StarCount[i]
        Dict[TypeList[i]]= Dictionary.copy()
    
    print('Total volume integrated is {} parsecs cubed'.format(totvolume))    
    return Dict
