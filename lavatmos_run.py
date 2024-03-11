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
# import lavatmos_vol3 as lavatmos
import lavatmos

# Import input
T_surf = float(sys.argv[1])
P_surf = float(sys.argv[2])
comp_name = sys.argv[3]
output_dir = sys.argv[4]

# Import composition
comp_fname = comp_path+comp_name+'.csv'
print(f'Magma composition read from: {comp_fname}')
comp_df = pd.read_csv(comp_fname,names=['spec','abund'])
comp = {}
for i in comp_df.index:
    comp[comp_df['spec'].loc[i]] = comp_df['abund'].loc[i]

# Initiate and run instance of LavAtmos
system = lavatmos.melt_vapor_system()
lavatmos_output = system.vaporise(T_surf, comp, P_melt=P_surf)

# Save results
output_name = 'degassed_partial_pressure.csv'
lavatmos_output.to_csv(output_dir+output_name)