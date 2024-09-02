# Standard modules
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import optimize
import collections
import sys

# Thermoengine modules
from thermoengine import equilibrate 
from thermoengine import model

# Local modules
from data.databases.janaf_data_importer_gef import janaf_data_importer
from data.databases.barin_data_importer_gef import barin_data_importer

class melt_vapor_system:
    
    '''
    
    Class in which the identities and properties of the melt vapor system 
    are saved. Also loads the important functions from thermoengine. Refrain 
    from initialising more than necessary due to time it takes to load 
    thermo data.
    
    '''
    
    def __init__(self):
        
        
        # Constants
        self._R = 8.314
        self.Tr = 298.15
        
        # Importing stoichiometries for reactions
        fname_cdef_values = 'cdef_values_vapor_reactions.csv'
        dname_cdef_values = 'data/'
        self.cdef = pd.read_csv(dname_cdef_values+fname_cdef_values)\
                    .set_index('vapor')
        
        # Importing thermo data
        self.thermo_data = janaf_data_importer() # janaf data
        self.thermo_data.update(barin_data_importer()) # barin data
        
        # Initialise states
        self.melt = MeltState()
        # self.vapor = VaporState()
        
        # Allowed input oxides
        self.oxides = ['SiO2','MgO','Al2O3','FeO','CaO','Na2O','K2O',\
                       'TiO2','Fe2O3']
        
    def vaporise(self, T, melt_comp, P_melt = 0.01,\
                          fO2_initial_guess = 1e-5,\
                          verbose = True):

        '''
        Main function with which the vaporisation reaction is run.
        
        
        Parameters
        ----------
        T : array, float, or int
            Equilibrium temperatures in Kelvin for which the vaporisation
            reactions must be run. Advised range: 1500 - 4000 K

        melt_comp : dict
            Composition of the melt in terms of percentage oxides.
        
        P_melt : float
            Pressure at which the melt equilibrium chemistry calculations are
            done. Advised range: 1e-3 - 100 bar (Results usually pressure 
            independent withing advised range)

        fO2_initial_guess : float
            Initial guess for the equilibrium value of the oxygen fugacity. 
            Used by the equation solver. Should not have to be changed. 
            
        Returns
        -------
        pressures : DataFrame
            The partial pressures of all the vapor species in the vapor above
            the melt for the given equilibrium temperatures.
        
        '''
        
        # Ensures that T is iterable even if just one value is given
        if type(T) != np.ndarray:
            T = np.array([T])
        
        # Calculate thermodynamic values
        self.melt_comp = self.melt.set_melt_comp(melt_comp,verbose)        
        self.melt.calculate_melt_chemical_potentials(T,P_melt,self.thermo_data)
        self.melt.calculate_melt_activities(T,P_melt)            
        logKr = self.logKr_calc(T,self.melt.mu0_liquid)
        
        # Calculate fO2 for all given temperatures
        fO2 = np.zeros(len(T))
        
        # TODO: clean this up
        if verbose:
            print('\nCalculating vapor pressures for given temperature values')
            # Progressbar settings
            with tqdm(total=len(T), file=sys.stdout) as pbar: 
                for i,t in enumerate(T):
                    pbar.set_description(f'Calculated')
                    
                    # Select thermodynamic values
                    a = self.melt.a.loc[t]
                    logKr_t = logKr.loc[t]

                    # Calculate fO2
                    results_opt = optimize.fsolve(self.mass_balance_equation,\
                                                  fO2_initial_guess,\
                                                  args=(a,logKr_t))
                    fO2[i] = results_opt[0]
                    
                    pbar.update(1) # Update progress bar
            
            print('Done!')
        
        else:
            for i,t in enumerate(T):                
                # Select thermodynamic values
                a = self.melt.a.loc[t]
                logKr_t = logKr.loc[t]

                # Calculate fO2
                results_opt = optimize.fsolve(self.mass_balance_equation,\
                                              fO2_initial_guess,\
                                              args=(a,logKr_t))
                fO2[i] = results_opt[0]
                

        # Save partial pressures 
        partial_pressures = self.calculate_partial_pressures(fO2,self.melt.a,\
                                                       logKr,final_output=True)
        partial_pressures.insert(0,'O2(g)',fO2) 
        
        return partial_pressures
        
        
    def logKr_calc(self, T, mu0_liquid):
        '''
        Calculates the logK values for all of the vapor reactions at a 
        given temperature.
        
        Assumes:

        logKr = logKf_products - v*logKf_reactants
        
        And since logKf of O2 is zero for all relevant T:

        logKr_vapor = logKf_vapor - ci*logKf_melt_oxide

        Parameters
        ----------
        T : array, float, or int
            Equilibrium temperature in kelvin.

        mu0_liquid : DataFrame
            Pandas dataframe containing Gibbs standard free energy values
            for all endmembers in the melt.

            
        Returns
        -------
        logK : dict
            LogK values for all vapor reactions.

        '''
        
        logKr = {}
        
        # Defines temperatures above which thermo values need to be 
        # extrapolated from JANAF values 
        T_interp = T[T<=6000]
        T_extrap = T[T>6000]

        # Interpolating for values within JANAF range
        gibbs_O2 = - T_interp * self.thermo_data['O2(g)'].gef(T_interp)\
                   + self.thermo_data['O2(g)'].DeltaH(self.Tr)

        # Extrapolating from last few values for T outside JANAF range 
        # (not advised)
        if T_extrap.size > 0:
            T_temp = self.thermo_data['O2(g)'].T 
            gibbs_temp = - T_temp * self.thermo_data['O2(g)'].gef(T_temp)\
                         + self.thermo_data['O2(g)'].DeltaH(self.Tr)
            gibbs_O2 = np.append(gibbs_O2,np.poly1d(np.polyfit(T_temp[-4:],\
                                                 gibbs_temp[-4:],1))(T_extrap))
        
        for i in range(len(self.cdef.index)):
            
            #Calculate gibbs energy of a vapor for a given T
            vapor = self.cdef.index[i]
            
            # Importing names and values from cdef file
            endmember = self.cdef.iloc[i]['endmember']
            liq_oxide1 = self.cdef.iloc[i]['liq_oxide1']
            
            if isinstance(liq_oxide1, str):
                liq_oxide1 = liq_oxide1.replace('(l)','')
            else:
                liq_oxide1 = 'None'
                
            liq_oxide2 = self.cdef.iloc[i]['liq_oxide2']
            if isinstance(liq_oxide2, str):
                liq_oxide2 = liq_oxide2.replace('(l)','')
            else: 
                liq_oxide2 = 'None'
            
            ci = self.cdef.iloc[i]['ci']
            di = self.cdef.iloc[i]['di']
            ei = self.cdef.iloc[i]['ei'] 
            fi = self.cdef.iloc[i]['fi']

            # Interpolating for values within JANAF range
            gibbs_vapor = - T_interp * self.thermo_data[vapor].gef(T_interp)\
                          + self.thermo_data[vapor].DeltaH(self.Tr)
            
            # Extrapolating from last few values for T outside JANAF range 
            # (not advised)
            if T_extrap.size > 0:
                T_temp = self.thermo_data[vapor].T 
                gibbs_temp = - T_temp * self.thermo_data[vapor].gef(T_temp)\
                             + self.thermo_data[vapor].DeltaH(self.Tr)
                gibbs_vapor = np.append(gibbs_vapor,\
                                        np.poly1d(np.polyfit(T_temp[-4:],
                                                  gibbs_temp[-4:],1))\
                                        (T_extrap))


            # Calculate K value of a vapor for a given T
            if vapor != 'O(g)':
                # Retrieving gibbs energies
                gibbs_endmember = np.array(mu0_liquid[endmember])

                if liq_oxide1 == 'None':
                    gibbs_liq_oxide1 = 0
                else:
                    gibbs_liq_oxide1 = np.array(mu0_liquid[liq_oxide1])

                if liq_oxide2 == 'None':
                    gibbs_liq_oxide2 = 0
                else:
                    gibbs_liq_oxide2 = np.array(mu0_liquid[liq_oxide2])
                
                logKr[i] = (-np.log10(np.e))/(self._R*T)\
                            * (gibbs_vapor\
                               - ci*gibbs_endmember\
                               - di*gibbs_O2\
                               - ei*gibbs_liq_oxide1\
                               - fi*gibbs_liq_oxide2)
            else: 
                logKr[i] = (-np.log10(np.e))/(self._R*T) * (gibbs_vapor\
                                                            - di * gibbs_O2)
        
        # Convert to dataframe
        logKr = pd.DataFrame(logKr,index = T)
        
        return logKr
    
    
    def mass_balance_equation(self,fO2,a,logKr):

        '''
        Equation that needs to be solved for fO2.
        
        Parameters
        ----------
        input : float
            fO2 in bar.

        a : dict
            Activities of oxides in the melt.
        
        logKr : dict
            LogK values for all vapor reactions.
            
        Returns
        -------
        eq1 : float
            The result of the mass balance equation.    
        '''
        
        # Calculate partial pressures
        partial_pressures = self.calculate_partial_pressures(fO2,a,logKr)
        
        # Calculate mass balance eq.
        mass_balance_eq = 0
        mass_balance_eq += fO2
        
        for i,gas in enumerate(self.cdef.index):
            mass_balance_eq += self.cdef['di'].iloc[i] * partial_pressures[i]
        
        return mass_balance_eq
    
    
    def calculate_partial_pressures(self,fO2,a,logKr,final_output=False): 
        '''
        Calculates the partial pressure all the vapor species.
        
        Parameters
        ----------
        fO2 : float
            Oxygen fugacity.

        a : dict
            Activities of oxides in the melt.
        
        logKr : dict
            LogK values for all vapor reactions.
            
        Returns
        -------
        pP : dict
            Partial pressures of all the vapor species
        
        '''

        partial_pressures_all = pd.DataFrame()

        
        # O does not require oxide, hence separate
        partial_pressures_all[0] = 10**logKr[0] * fO2**self.cdef['di'].iloc[0]

        for i in range(1,len(self.cdef.index)):
            
            # Importing names and values from cdef file
            endmember = self.cdef.iloc[i]['endmember']
            liq_oxide1 = self.cdef.iloc[i]['liq_oxide1']
            
            if isinstance(liq_oxide1, str):
                liq_oxide1 = liq_oxide1.replace('(l)','')
            else:
                liq_oxide1 = 'None'
                
            liq_oxide2 = self.cdef.iloc[i]['liq_oxide2']
            if isinstance(liq_oxide2, str):
                liq_oxide2 = liq_oxide2.replace('(l)','')
            else: 
                liq_oxide2 = 'None'
            
            ci = self.cdef.iloc[i]['ci']
            di = self.cdef.iloc[i]['di']
            ei = self.cdef.iloc[i]['ei'] 
            fi = self.cdef.iloc[i]['fi']
            
            # Retrieving activities            
            if endmember in a:
                a_endmember = a[endmember]
            else: 
                a_endmember = 0
            
            if liq_oxide1 == 'None':
                a_liq_oxide1 = 1
            else:
                a_liq_oxide1 = a[liq_oxide1]

            if liq_oxide2 == 'None':
                a_liq_oxide2 = 1
            else:
                a_liq_oxide2 = a[liq_oxide2]
            
            # Calculating partial pressure
            partial_pressures_all[i] = 10**logKr[i]\
                                  * a_endmember**ci\
                                  * fO2**di\
                                  * a_liq_oxide1**ei\
                                  * a_liq_oxide2**fi

        # Ensures that if this is the final output, all pressures of same 
        # species are summed. Necessary because there may be more than one 
        # reaction producing the same species (for Fe for example).
        # TODO: DEPRECATED - REMOVE
        if final_output:
            partial_pressures = pd.DataFrame()
            for i,vapor in enumerate(self.cdef.index):
                if vapor in partial_pressures.columns:
                    partial_pressures[vapor] += partial_pressures_all[i]
                else:
                    partial_pressures[vapor] = partial_pressures_all[i]

            return partial_pressures

        else:
            return partial_pressures_all
        
        
class MeltState:
    '''
    Class for holding the vapor state for the thermodynamic parameters 
    determined by melts. 
    
    Parameters
    ----------
    
    Attributes
    ----------
    
    '''

    def __init__(self,melts_version='1.0.2'): 

        # Constants
        self._R = 8.314 # Gas constant
        self.Tr = 298.15

        # Initialising parameters to be calculated
        self.a = None
        
        # Initialising thermoengine classes
        self.melts = equilibrate.MELTSmodel(melts_version)
        modelDB = model.Database(liq_mod='v1.0')
        self.liq_phs = modelDB.get_phase('Liq')
        self.endmember_names = self.liq_phs.endmember_names
        self.nonendmember_oxide_names = ['MgO','FeO','CaO','Na2O',\
                                         'K2O','Cr2O3']
        
        # Importing data for liquid reactions
        fname_cdef_liq_values = 'cdef_values_liquid_reactions.csv'
        dname_cdef_liq_values = 'data/'        
        self.cdef_liq = pd.read_csv(dname_cdef_liq_values\
                                + fname_cdef_liq_values).set_index('liq_oxide')
        
        # Oxides that are included in calculations
        self.used_oxides = ['SiO2','MgO','Al2O3','FeO','CaO','Na2O','K2O',\
                            'TiO2','Fe2O3']

        
    def set_melt_comp(self,input_comp,verbose):
        '''
        Checks if input composition is valid before setting melt 
        composition in thermoengine class (MELTS). 
        
        Parameters
        ----------
        input_comp : dict
            Input composition.

        verbose : bool
            True/False for printing text or not. 

        Returns
        -------
        melt_comp : dict
            Melt composition as accepted by MELTS.
                
        '''

        # Set dictionary labels to oxides allowed by melts
        melt_comp = {} # dict([(ox,0) for ox in self.liq_phs.OXIDES])

        for species in input_comp:

            if species not in self.liq_phs.OXIDES: # Check validity of input
                print(f'ERROR: {species} not allowed. Composition not set.')
                return  
            else: 
                melt_comp[species] = input_comp[species]
        
            if species not in self.used_oxides:
                print(f'WARNING: input value for {species} passed. Species not'
                      + f' (yet) included in vaporisation calculations.')
            
        # Set melt comp in thermoengine class
        self.melts.set_bulk_composition(melt_comp)
        if verbose:
            print(f'Melt composition set to: {melt_comp}')

        # Save included oxides
        self.included_oxides = {}
        for ox in melt_comp:
            if melt_comp[ox] == 0:
                self.included_oxides[ox] = False
            else:
                self.included_oxides[ox] = True

        return melt_comp
    
    
    def calculate_melt_chemical_potentials(self, T, P_melt, thermo_data):
        '''
        Calculates the chemical potentials of the oxide species in the melt. 
        For the oxides that are not included in MELTS as endmemebers, data 
        from the thermo databases is used. 
        
        Parameters
        ----------
        T : float
            Equilibrium temperature in kelvin
            
        P_melt : float
            Pressure at which the melt equilibrium chemistry calculations are
            done.
            
        thermo_data : dict
            Dictionary containing Gibbs energy function values for each species 
            relevant to the calculations.
                
        '''
        # Calculates endmember chemical potentials using MELTS function
        mu0_endmember = np.array([self.liq_phs.gibbs_energy(T,P_melt,mol=imol)\
                                  for imol in np.eye(15)])
        self.mu0_liquid = pd.DataFrame(mu0_endmember,\
                                index=self.liq_phs.endmember_names,columns=T).T

        for oxide in self.nonendmember_oxide_names:

            # Try to interpolate, if outside range of database uses linear fit 
            # of last four data points to extrapolate to higher temperatures. 
            max_temp = thermo_data[oxide+'(l)'].T.iloc[-1]
            T_interp = T[T<=max_temp]
            T_extrap = T[T>max_temp]

            gibbs = - T_interp * thermo_data[oxide+'(l)'].gef(T_interp) \
                    + thermo_data[oxide+'(l)'].DeltaH(self.Tr)
            if T_extrap.size > 0: 
                T_temp = thermo_data[oxide+'(l)'].T 
                gibbs_temp = - T_temp * thermo_data[oxide+'(l)'].gef(T_temp)\
                             + thermo_data[oxide+'(l)'].DeltaH(self.Tr)
                gibbs = np.append(gibbs,np.poly1d(np.polyfit(T_temp[-4:],
                                        gibbs_temp[-4:],1))(T_extrap))

            self.mu0_liquid[oxide] = gibbs
            
    
    def calculate_melt_activities(self, T, P_melt):
        '''
        Retrieves the activity of the endmembers using a Thermoengine (MELTS)
        function.
        
        Parameters
        ----------
        T : float
            Equilibrium temperature in kelvin
            
        P_melt : float
            Pressure at which the melt equilibrium chemistry calculations are
            done.                
        '''
        a = {}
        # Calculates endmember activities using MELTS function
        for t in T:
            # Calculate excess chemical potential
            output = self.melts.equilibrate_tp(t,P_melt,initialize=True)
            status,t,p,xmlout = output[0]
            print('OUTPUT:',status,t,p)
            excs_dict = self.melts.get_thermo_properties_of_phase_components(xmlout,\
                                'Liquid',mode='excess')
            print(excs_dict)
            included_endmembers = list(excs_dict.keys())
            excs = list(excs_dict.values())            
            
            # Convert to activity
            a[t] = np.exp(excs/(self._R*t))
        
        self.a = pd.DataFrame(a,index=included_endmembers).T