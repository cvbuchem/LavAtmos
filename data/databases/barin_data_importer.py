import pandas as pd
from scipy.interpolate import interp1d

def barin_data_importer():
    '''
    Imports Barin (1997) data for specified species.
    
    Returns
    -------
    logK : dict
        LogK data for species.
    
    '''

    species = ['K2O(l)']
    logK = {}

    for spec in species:
        data = pd.read_csv(f'data/databases/barin/data/barin_data_{spec}.csv')
        logK[spec] = interp1d(data['T'],data[spec])

    return logK