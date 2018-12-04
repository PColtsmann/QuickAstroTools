import numpy as np
import mpmath as mp

def contprob(flux,wavelength,separation,survey):#from Carniani et al. 2015. A&A. 584. A78
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
    
    if wavelength == 0.87: #from Stach et al. 2018. ApJ. 860. 161.
        #Ideally should be used between 2 and 8 mJy.
        N0 = 1200 #per square degree
        S0 = 5.1 #mJy
        siglim = S/S0
        a = 5.9
        b = 0.4
        N = integrate.quad(lambda sig: (N0)*((sig**a + sig**b)**-1), siglim, +np.inf)
        n = N[0]/12960000 # put into square arcseconds from square degrees
        
    
    
    lam = n*area #average number of galaxies in your area of interest)
    prob0 = (np.e**(-1*lam)) #probability of no galaxies, poisson distribution
    prob = 1 - prob0 #probability of at least one galaxy contaunit option - deg or arcsecminating a source
    probsurvey = prob * survey #expected contaminated sources, binomial distribution expectation
    print(' ')
    print('Chance that there is a contaminating galaxy in THIS image is '+str(100*prob)+'%')
    print(' ')
    print('Chance that there is a contaminating galaxy in AN image is '+str(100*(1-(1-prob)**survey))+'%')
    print('')
    print('Expected number of contaminated images out of '+str(survey) +' is '+str(probsurvey))
    print('')
    print('Expected number of galaxies in THIS image is ' +str(lam))
    print(' ')
    return [float(prob), float(probsurvey)]

"""
TODO
Include probabilities for multiple contaminating galaxies i.e. chance of multiple unrelated or multiplicity e.g. https://arxiv.org/pdf/1304.4266.pdf
Rewrite for {} use
restrict sig figs
Galaxy in a ring option?
opt arg for area
unit option - deg or arcsec
make schechter and dpl int/no int separate functions
include more papers specifications - based on galaxy range?
include errors
make survey an optional argument
take care of multiplicity
add warnings for being above/below recommended range 
"""
