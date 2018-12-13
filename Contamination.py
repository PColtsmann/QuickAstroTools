import numpy as np
import mpmath as mp
import scipy.integrate as integrate
from astropy import units as u

def contprob(flux, wavelength, separation, survey = None):
    
    """
    
    Provides an estimate for probabilities of background galaxies contaminating your unbiased pointed surveys.
    Uses Carniani et al. 2015. A&A. 584. A78 and Stach et al. 2018. ApJ. 860. 161.
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    flux:                Limiting detection flux in astropy units of choice or mJy if left unitless.
    
    wavelength:          Wavelength of observation in mm, only accepts 1.1 or 1.3 or 0.87.
    
    separation:          Minimum separation between centre of image and contaminating galaxy in astropy unit of choice or arcsec if left blank.
    
    survey:              Number of stars in your survey
    
    --------------------------------------------------------------------------------------------------------
    
    """
    
    if wavelength != 1.1 and  wavelength != 1.3 and wavelength != 0.87:
        raise ValueError('Wavelength must be 1.3, 1.1 or 0.87')
        
    if type(flux*u.millijansky) == type(u.millijansky*u.millijansky):# avoids Unit != IrreducibleUnit != PrefixUnit by comparing for composite units
        S = flux.to(u.millijansky).value
    else:
        S = flux
    
    if type(separation*u.degree) == type(u.degree*u.degree): 
        area = np.pi*(separation.to(u.degree).value)**2 
    else:
        area = np.pi*(separation*u.arcsec.to(u.degree).value)**2 
        
    if wavelength == 1.1 or wavelength == 1.3: #using Schechter function from Carniani et al. 2015. A&A. 584. A78
        if wavelength == 1.3:
            #Ideally should be used above 0.06 mJy.
            phi = 1.8e3 #per square degree
            S0 = 1.7 #mJy
            a = -2.08
        elif wavelength == 1.1:
            #Ideally should be used above 0.1 mJy.
            phi = 2.7e3
            S0 = 2.6
            a = -1.81        
        
        siglim = S/S0
        
        N = phi*(mp.gammainc(a+1,a=siglim,b='inf')) # integrating the Schechter function gives the incomplete gamma function scaled by phi
    
    if wavelength == 0.87: #using double power law from Stach et al. 2018. ApJ. 860. 161.
        #Ideally should be used between 2 and 8 mJy.
        N0 = 1200 #per square degree
        S0 = 5.1 #mJy
        siglim = S/S0
        a = 5.9
        b = 0.4
        N = integrate.quad(lambda sig: (N0)*((sig**a + sig**b)**-1), siglim, +np.inf)

    lam = N[0]*area #average number of galaxies in your area of interest
    prob0 = (np.e**(-1*lam)) #probability of no galaxies, poisson distribution
    prob = 1 - prob0 #probability of at least one galaxy contaminating a source
    
    if survey:
        probsurvey = prob * survey #expected contaminated sources, binomial distribution expectation
        print(' ')
        print('Chance that there is a contaminating galaxy in THIS image is {} %'.format(100*prob))
        print(' ')
        print('Chance that there is a contaminating galaxy in AN image is {} %'.format(100*(1-(1-prob)**survey)))
        print('')
        print('Expected number of contaminated images out of {} is {}'.format(survey,probsurvey))
        print('')
        print('Expected number of galaxies in THIS image is {}'.format(lam))
        print(' ')
        return [float(prob), float(probsurvey)]
    else:
        print(' ')
        print('Chance that there is a contaminating galaxy in THIS image is {} %'.format(100*prob))
        print('')
        print('Expected number of galaxies in THIS image is {}'.format(lam))
        print(' ')
        return float(prob)
