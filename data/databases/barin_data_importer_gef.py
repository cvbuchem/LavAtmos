import pandas as pd
from scipy.interpolate import interp1d

class barin_data_class:
    
    def __init__(self,gef,DeltaH,T):
        
        self.gef = gef
        self.DeltaH = DeltaH
        self.T = T

def barin_data_importer():
    '''
    Imports Barin (1997) data for specified species.
    
    Returns
    -------
    logK : dict
        LogK data for species.
    
    '''

    species = ['K2O(l)']
    barin_data = {}

    for spec in species:
        data = pd.read_csv(f'data/databases/barin/data/barin_data_{spec}_gef.csv')
                
        gef = interp1d(data['T'],data['gef'])
        DeltaH = interp1d(data['T'],data['dHf'])
        T = data['T']

        barin_data[spec] = barin_data_class(gef,DeltaH,T)
        
    return barin_data