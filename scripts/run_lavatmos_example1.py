# Standard python packages
import os
import numpy as np
import pandas as pd
import sys

# LavAtmos
mod_dir = '/home/jovyan/ThermoEngine/LavAtmos'
sys.path.append(mod_dir)
import lavatmos

# Import compositions
os.chdir(mod_dir)
vf13_comps_df = pd.read_csv('/home/jovyan/ThermoEngine/LavAtmos/data/input/vf2013_comps.csv',index_col=0)
vf13_comps = {}
for name in vf13_comps_df.columns:
    vf13_comps[name] = vf13_comps_df[name].to_dict()
print(vf13_comps['BSE'])

# Initiate instance of LavAtmos
system = lavatmos.melt_vapor_system()

# Temperature values
# T = np.arange(1500,4050,50)
T = np.array([2000,3000])

# Run the vaporisation calculations
lavatmos_bse = system.vaporise(T, vf13_comps['BSE'])

# Save results
output_dir = 'output/'
name = 'script_example1_output.csv'
print(f'Saving results to: {output_dir+name}')
lavatmos_bse.to_csv(output_dir+name)