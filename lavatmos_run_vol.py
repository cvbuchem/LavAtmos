# Standard python packages
import os 
import numpy as np 
import pandas as pd
import sys
import warnings
warnings.filterwarnings("ignore")

# Local packages and paths
sys.path.append(os.getcwd())
from input.paths import paths_importer
paths = paths_importer()
comp_path = paths.lava_comps
wkdir = paths.lavatmos_dir+'ThermoEngine/LavAtmos'
os.chdir(wkdir)
sys.path.append(wkdir)

import lavatmos
import lavatmos_vol

# Import input
T_boa = float(sys.argv[1])
P_volatile = float(sys.argv[2])
comp_melt_name = sys.argv[3]
comp_vol_name = sys.argv[4]
output_dir = sys.argv[5]

# Importing melt composition
comp_fname = comp_path+comp_melt_name+'.csv'
print(f'Magma composition read from: {comp_fname}')
comp_df = pd.read_csv(comp_fname,names=['spec','abund'])
comp = {}
for i in comp_df.index:
    comp[comp_df['spec'].loc[i]] = comp_df['abund'].loc[i]

# Running for melt without volatile atmosphere
if comp_vol_name == 'None':
    print('Running LavAtmos without volatiles')
    # Initiate and run instance of LavAtmos
    system = lavatmos.melt_vapor_system()
    lavatmos_output = system.vaporise(T_boa, comp)

    # Save results
    output_name = 'degassed_partial_pressure.csv'
    lavatmos_output.to_csv(output_dir+output_name)

# Runnning for melt with volatile atmosphere
else: 
    print('Running LavAtmos with volatiles')
    comp_vol = pd.read_csv(comp_vol_name,index_col=0).to_dict()['mole_fraction']

    # Initiate and run instance of LavAtmos
    system = lavatmos_vol.melt_vapor_system()
    lavatmos_output = system.vaporise(T_boa, P_volatile, comp, comp_vol)

    # Save results
    output_name = 'degassed_partial_pressure.csv'
    lavatmos_output.to_csv(output_dir+output_name)




