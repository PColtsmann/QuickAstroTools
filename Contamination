import numpy as np
import mpmath as mp

def schechter(flux,wavelength,separation,survey):#from Carniani et al. 2015. A&A. 584. A78
    """
    
    Provides an estimate for probabilities of background galaxies contaminating your unbiased ALMA surveys.
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    Flux:                Limiting detection flux in mJy.
    
    Wavelength:          Wavelength of observation in mm, only accepts 1.1 or 1.3.
    
    Separation:          Minimum separation between centre of image and contaminating galaxy in arcsec.
    
    Survey:              Number of stars in your survey
    
    --------------------------------------------------------------------------------------------------------
    
    """
    S = flux # mJy
    area = np.pi*(separation)**2 #area in sq arcsec
    
    if wavelength == 1.3:
        phi = 1.8e3
        Sstar = 1.7
        a = -2.08
    elif wavelength == 1.1:
        phi = 2.7e3
        Sstar = 2.6
        a = -1.81        
    else:
        print('Only "1.1" and "1.3" are accepted')
        return()
    
    siglim = S/Sstar
    
    N = phi*(mp.gammainc(a+1,a=siglim,b='inf'))#square degrees
    n = N/12960000
    
    lam = n*area #average number of galaxies in your area of interest)
    prob0 = (np.e**(-1*lam)) #probability of no galaxies, poisson distribution
    prob = 1 - prob0 #probability of at least one galaxy contaminating a source
    probsurvey = prob * survey #expected contaminated sources, binomial distribution expectation
    print(' ')
    print('Chance that there is a contaminating galaxy in THIS image is '+str(100*prob)+'%')
    print(' ')
    print('Chance that there is a contaminating galaxy in AN image is '+str(100*(1-(1-prob)**survey))+'%')
    print('')
    print('Expected number of contaminated images out of '+str(survey) +' is '+str(probsurvey))
    print('')
    return [float(prob), float(probsurvey)]
