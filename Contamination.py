import numpy as np
import mpmath as mp

def contprob(flux,wavelength,separation,survey):
    """
    
    Provides an estimate for probabilities of background galaxies contaminating your unbiased pointed surveys.
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    Flux:                Limiting detection flux in mJy.
    
    Wavelength:          Wavelength of observation in mm, only accepts 1.1 or 1.3 or 0.87.
    
    Separation:          Minimum separation between centre of image and contaminating galaxy in arcsec.
    
    Survey:              Number of stars in your survey
    
    --------------------------------------------------------------------------------------------------------
    
    """
    if wavelength != 1.1 and  wavelength != 1.3 and wavelength != 0.87:
        raise ValueError('Wavelength must be 1.3, 1.1 or 0.87')
    
    S = flux # mJy
    area = np.pi*(separation)**2 #area in sq arcsec
    
    if wavelength == 1.1 or wavelength == 1.3: #from Carniani et al. 2015. A&A. 584. A78
        if wavelength == 1.3:
            #Ideally should be used above 0.06 mJy.
            phi = 1.8e3 #per square degree
            Sstar = 1.7 #mJy
            a = -2.08
        elif wavelength == 1.1:
            #Ideally should be used above 0.1 mJy.
            phi = 2.7e3
            Sstar = 2.6
            a = -1.81        
        
        siglim = S/Sstar
        
        N = phi*(mp.gammainc(a+1,a=siglim,b='inf')) # integrating the Schechter function gives the incomplete gamma function scaled by phi
        n = N/12960000 # put into square arcseconds from square degrees
    
    if wavelength == 0.87: #from Simpson et al. 2015. ApJ. 807. 128.
        #Ideally should be used between 2 and 8 mJy.
        N0 = 390 #per square degree
        S0 = 8.4 #mJy
        siglim = S/S0
        a = 1.9
        b = 10.5
        
        N = (N0/S0)*((siglim**a + siglim**b)**-1)
        n = N/12960000 # put into square arcseconds from square degrees
    
    
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
    


"""
TODO
Include probabilities for multiple contaminating galaxies i.e. chance of multiple unrelated or multiplicity e.g. https://arxiv.org/pdf/1304.4266.pdf
Include expected number of galaxies in an image.
"""
