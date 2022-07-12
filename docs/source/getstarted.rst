Get started
============

.. _getstarted:

Starting a Jupyter Lab environment
----------------------------------

In order to get started with LavAtmos we first need to run the docker
image that contains Thermoengine (see :doc:`installation`). Open your
terminal and change your working directory to wherever you cloned
Thermoengine::

   cd path/to/Thermoengine

Note: make sure to change ``path/to/`` to the actual path on your system.

Next, run the following command::

   ./run_docker_locally.sh

If needed (check with your system admin if allowed!) run the above 
command using ``sudo``. If working correctly, you should get a URL in your 
terminal that you can then copy and paste into a web browser. This will bring
you to a Jupyter lab environment in which you will now be able to run LavAtmos.

Importing LavAtmos
------------------

In order to run the instructions given below directly from a Jupyter
notebook file, check out the `notebook directory <https://github.com/cvbuchem/LavAtmos/tree/master/notebooks>`_
on the Github page.

First of all, ensure that you are working within the LavAtmos directory::

   import os
   os.chdir('path/to/Thermoengine/LavAtmos')

Note: Make sure to change ``path/to/`` to the actual path in your system. 

You should then be able to import LavAtmos using::
   
   import lavatmos

Defining a melt composition
---------------------------

The melt compositions can be passed to LavAtmos in the form of dictionaries (just as for MELTS).
An often used composition is Bulk Silicate Earth (BSE)::

   comp_BSE = {'SiO2': 45.4,\
               'MgO': 36.76,\
               'Al2O3': 4.48,\
               'TiO2': 0.21,\
               'FeO': 8.1,\
               'CaO': 3.65,\
               'Na2O': 0.349,\
               'K2O': 0.031}

The LavAtmos repository also contains some preset compositions which can
be imported directly using::

   vf13_comps_df = pd.read_csv('/home/jovyan/ThermoEngine/LavAtmos/data/input/vf2013_comps.csv',index_col=0)
   vf13_comps = {}
   for name in vf13_comps_df.columns:
       vf13_comps[name] = vf13_comps_df[name].to_dict()
   print(vf13_comps['BSE'])

Note: MELTS, which is used to calculate the melt thermodynamics, may not always be able to converge on a composition. 
It generally works well for naturally occuring compositions.

Running the equilibrium calculations
------------------------------------

Initialize a LavAtmos object::

   system = lavatmos.melt_vapor_system()

Define the temperatures for which you want to perform the calculations::

   T = np.arange(1500,4050,50)

And perform the gas-melt equilibrium calculations using the vaporise function::

   results = system.vaporise(T, comp_BSE)
   print(results)
