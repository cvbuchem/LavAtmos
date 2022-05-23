# LavAtmos

This code performs gas-melt equilibrium calculations for a given temperature and melt composition. The thermodynamics of the melt are modelled by the MELTS (Ghiorso and Sack 1995) code as presented in the Thermoengine package.

## Installation

Since this code uses the MELTS code as provided by the [Thermoengine](https://enki-portal.gitlab.io/ThermoEngine/) package, this needs to be installed before LavAtmos may be used. 

We suggest to follow the instructions on the [Thermoengine GitLab page](https://gitlab.com/ENKI-portal/ThermoEngine) on how to run a container image locally. 

Once you have installed Thermoengine, clone this repository into the Thermoengine directory.

## Usage

Ensure that you are working within the LavAtmos directory. Import LavAtmos:

``import LavAtmos``

You should then define a composition. An often used example is bulk silicate earth (BSE):

``comp_BSE = {'SiO2': 45.4, 'MgO': 36.76, 'Al2O3': 4.48, 'TiO2': 0.21, 'FeO': 8.1, 'CaO': 3.65, 'Na2O': 0.349, 'K2O': 0.031} ``

Initialise a LavAtmos object:

``system = lavatmos.melt_vapor_system()``

Define a temperature range:

``T = np.arange(1500,4050,50)``

And perform the gas-melt equilibrium calculations using the vaporise function:

``results = system.vaporise(T, comp_BSE)``

Checkout the example notebooks in the notebook directory for more comprehensive examples on how to use LavAtmos.

## Troubleshooting

If you are experiencing any issues when installing or using LavAtmos, please don't hesitate to open a new issue on Github. If you do, please leave as many details about your issue/error as you can.


## License and citation

License: [MIT](https://choosealicense.com/licenses/mit/)

If you use LavAtmos results, please cite:

- Van Buchem (in prep) - will add link to paper as soon as it is published

Since LavAtmos makes use of MELTS, MELTS should also be cited as specified on [the MELTS webpage](https://melts.ofm-research.org/).