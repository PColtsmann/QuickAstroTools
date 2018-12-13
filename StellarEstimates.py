import numpy as np
from astropy import units as u

def create_dictionary(filename):
    
    """
    
    Constructs a dictionary of spectral types for use in the following functions.
    
    Sourced from http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
    MamajekTable.txt in this repository can be used.
    
    Asterix/* have been removed from logAges
    ellipses/... have been replaced by inf
    '19.52:' in T9.5V M_J has been replaced with '19.52'
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    filename: The name/location of the .txt file containing the stellar parameters.
              The file should be that linked above with the above adjustments made or MamajekTable.txt in this repository.

    --------------------------------------------------------------------------------------------------------
    
    """
    
    file = np.loadtxt(filename, dtype=str)
    
    Dictionary = {}
    for i in range(0,len(file)):
        Dict = {}
        Dict['SpT'] = file[i][0]
        Dict['Teff'] = float(file[i][1])
        Dict['logT'] = float(file[i][2])
        Dict['BCv'] = float(file[i][3])
        Dict['Mv'] = float(file[i][4])
        Dict['logL'] = float(file[i][5])
        Dict['B-V'] = float(file[i][6])
        Dict['Bt-Vt'] = float(file[i][7])
        Dict['V-G'] = float(file[i][8])
        Dict['U-B'] = float(file[i][9])
        Dict['V-Rc'] = float(file[i][10])
        Dict['V-Ic'] = float(file[i][11])
        Dict['V-Ks'] = float(file[i][12])
        Dict['J-H'] = float(file[i][13])
        Dict['H-Ks'] = float(file[i][14])
        Dict['Ks-W1'] = float(file[i][15])
        Dict['W1-W2'] = float(file[i][16])
        Dict['W1-W3'] = float(file[i][17])
        Dict['W1-W4'] = float(file[i][18])
        Dict['Msun'] = float(file[i][19])
        Dict['logAge'] = float(file[i][20])
        Dict['b-y'] = float(file[i][21])
        Dict['M_J'] = float(file[i][22])
        Dict['M_Ks'] = float(file[i][23])
        Dict['Mbol'] = float(file[i][24])
        Dict['i-z'] = float(file[i][25])
        Dict['z-Y'] = float(file[i][26])
        Dict['R_Rsun'] = float(file[i][27])
        Dictionary[file[i][0]] = Dict
    return Dictionary

def search(value, column, dictionary,thresh = 0.2):
    
    """
    
    Searches for stellar types that have similar properties to that specified.
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    value: The approximate value of the property you are trying to match to a type.
    
    column: The type of property 'value' corresponds to, a string defining the column in the dictionary.
    
    dictionary: The dictionary with a column 'column' to be searched for entries with values close to 'value'. See create_dictionary.
    
    thresh: The factor to which you wish your matches to be accurate by. e.g. thresh = 0.2 will return matches within 20% of 'value'.
            Larger factors are needed if 'value' is close to 0.

    --------------------------------------------------------------------------------------------------------
    
    """
    
    matches = []
    residuals = []
    closest = []
    
    for i in dictionary:
        testvalue = dictionary[i][str(column)]
        residual = testvalue - value
        absresidual = abs(residual)
        residuals.append(absresidual)
        if absresidual <= value*thresh:
            matches.append([i,testvalue,residual])
    
    for i in range(0,len(matches)):
        if abs(matches[i][2]) == min(residuals):
            closest.append(matches[i][0])
    
    return matches,closest
            
def stellar_distance(dictionary, SpT=None, distance=None, magnitude=None, thresh=0.1):
    
    """
    
    Given two of spectral type, distance and apparent magnitude will calculate the other.
    
    ---------------------------------------------------------------------------------------------------------
    
    Parameters:
        
    dictionary: The dictionary with which your parameters will be compared. See create_dictionary.
    
    SpT: The spectral type of the star you wish to find the distance/magnitude of. If left blank, type/s with most similar properties will be calculated to within factor 'thresh'.
         Must be included in dictionary, see your dictionary or http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt for all available.
         e.g. 'A0V'
    
    distance: The distance to the star you wish to find the magnitude/type of, as an astropy Quantity or in parsecs if left unitless. If left blank or only astropy units are specified, this will be calculated in parsecs or specified units.
    
    magnitude: The magnitude of the star you wish to find the distance/type of. If left blank, this will be calculated.
    
    thresh: The factor to which you wish your spectral type matches to be accurate by. e.g. thresh = 0.2 will return matches within 20% of your specified parameters.
            Larger factors are needed if this value is close to 0.

    --------------------------------------------------------------------------------------------------------
    
    """
    print('')
    
    if SpT != None and SpT not in dictionary:
        raise ValueError("That spectral type is not supported, remember to include luminosity class, e.g. 'A0V', else please refer to http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt for all available types")
    
    elif distance != None and SpT and not magnitude:
        if type(distance) == type(1*u.parsec):
            distance = distance.to(u.parsec).value
        print('Calculated apparent magnitude is {}'.format(dictionary[SpT]['Mv'] + 5*(np.log10(distance)-1)))
        return 
    
    elif SpT and (magnitude or magnitude == 0) and (type(distance) != float or type(distance) != int):
        if distance == None:
            print('Calculated distance is {}'.format(10*np.sqrt(100**((magnitude-dictionary[SpT]['Mv'])/5))*u.parsec))
            return
        elif type(distance) == type(u.parsec):
            print('Calculated distance is {}'.format((10*np.sqrt(100**((magnitude-dictionary[SpT]['Mv'])/5))*u.parsec).to(distance)))
            return
        else:
            print('Please enter two and only two of SpT, distance and magnitude')
    
    elif distance != None and (magnitude or magnitude ==0) and not SpT:
        if type(distance) == type(1*u.parsec):
            distance = distance.to(u.parsec).value        
        Mv = magnitude - 5*(np.log10(distance)-1)
        print('Spectral type estimates, their absolute maginitudes and their residuals are {}'.format(search(Mv,'Mv',dictionary,thresh)[0]))
        print('')
        print('Closest spectral type estimate is {}'.format(search(Mv,'Mv',dictionary,thresh)[1]))
        return 
    
    else:
        print('Please enter two and only two of SpT, distance and magnitude')
