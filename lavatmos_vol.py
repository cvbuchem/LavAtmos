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

    def __init__(self,melts_version='1.0.2'):
        
        
        # Constants
        self._R = 8.314
        self.Tr = 298.15

        # Volatile species
        self.all_vol_species = ['H2O','H2','CO2','CO','SO2','SO','NO','N2','P','PO']
        oxidized_vol_species = ['H2O','CO2','SO2','NO','PO']
        non_oxidized_vol_species = ['H2','CO','SO','N2','P']
        
        # Importing stoichiometries for reactions
        fname_cdefg_values = 'cdefg_values_vapor_reactions.csv'
        fname_cd_values_vol = 'cd_values_volatile_reactions.csv'
        dname_cdef_values = 'data/'
        self.cdefg = pd.read_csv(f'data/{fname_cdefg_values}').set_index('vapor')
        self.cd = pd.read_csv(f'data/{fname_cd_values_vol}').set_index('oxidized')

        # Importing thermo data
        self.thermo_data = janaf_data_importer() # janaf data
        self.thermo_data.update(barin_data_importer()) # barin data
        
        # Initialise states
        self.melt = MeltState(melts_version)

    def vaporise(self, T, P_volatile, melt_comp, volatile_comp,\
                          fO2_initial_guess = 1e-10,\
                          verbose = True):

        '''
        Main function with which the vaporisation reaction is run.
        
        
        Parameters
        ----------
        T : array, float, or int
            Equilibrium temperatures in Kelvin for which the vaporisation
            reactions must be run.
        
        P_melt : float
            Pressure at which the melt equilibrium chemistry calculations are
            done.

        fO2_initial_guess : float
            Initial guess for the equilibrium value of the oxygen fugacity. 
            Used by the equation solver.
            
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
        self.melt.calculate_melt_chemical_potentials(T,P_volatile,self.thermo_data)
        self.melt.calculate_melt_activities(T,P_volatile)            
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
                                                  args=(a,logKr_t,P_volatile,volatile_comp))
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
                                              args=(a,logKr_t,P_volatile,volatile_comp))
                fO2[i] = results_opt[0]

        # Save partial pressures 
        partial_pressures, partial_pressures_volatiles = self.calculate_partial_pressures(fO2,self.melt.a,\
                                                       logKr,P_volatile, volatile_comp,final_output=True)
        partial_pressures = pd.concat([partial_pressures,partial_pressures_volatiles], axis=1)
        partial_pressures.insert(0,'O2(g)',fO2) 

        return partial_pressures#, self.melt.a, logKr
        
        
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

        TODO: Rewrite in terms of matrix multiplications.
              If needed, add option to use MELTS chempot values.
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


        # Calculating gibbs energies for volatiles species (TODO: ADD EXTRAP)
        gibbs_volatiles = {}
        for spec in self.all_vol_species:
            gibbs_volatiles[spec+'(g)'] = - T_interp\
                                    * self.thermo_data[f'{spec}(g)'].gef(T_interp)\
                                    + self.thermo_data[f'{spec}(g)'].DeltaH(self.Tr)

        for i,oxidised_spec in enumerate(self.cd.index):

            # Importing names and values from cdef file
            reduced_spec = self.cd.iloc[i]['reduced']
            ci = self.cd.iloc[i]['ci']
            di = self.cd.iloc[i]['di']
 
            # Retrieving gibbs energies
            gibbs_reduced_spec = np.array(gibbs_volatiles[reduced_spec])

            # Calculate K value of a vapor for a given T            
            logKr[oxidised_spec] = (-np.log10(np.e))/(self._R*T)\
                        * -(gibbs_volatiles[oxidised_spec]\
                           - ci*gibbs_reduced_spec\
                           - di*gibbs_O2)

        #Calculate gibbs energy of a vapor for a given T    
        for i, vapor in enumerate(self.cdefg.index):    
            
            # Importing names and values from cdef file
            endmember = self.cdefg.iloc[i]['endmember']
            liq_oxide1 = self.cdefg.iloc[i]['liq_oxide1'].replace('(l)','')
            liq_oxide2 = self.cdefg.iloc[i]['liq_oxide2'].replace('(l)','')
            ci = self.cdefg.iloc[i]['ci']
            di = self.cdefg.iloc[i]['di']
            ei = self.cdefg.iloc[i]['ei'] 
            fi = self.cdefg.iloc[i]['fi']
            gi = self.cdefg.iloc[i]['gi']

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
                
                logKr[vapor] = (-np.log10(np.e))/(self._R*T)\
                            * (gibbs_vapor\
                               - ci*gibbs_endmember\
                               - di*gibbs_O2\
                               - ei*gibbs_volatiles['H2(g)']\
                               - fi*gibbs_liq_oxide1\
                               - gi*gibbs_liq_oxide2)
            else: 
                logKr[vapor] = (-np.log10(np.e))/(self._R*T) * (gibbs_vapor\
                                                            - di * gibbs_O2)
        
        # Convert to dataframe
        logKr = pd.DataFrame(logKr,index = T)
        
        return logKr
    
    
    def mass_balance_equation(self,fO2,a,logKr,P_volatile,volatile_comp):

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
        partial_pressures, partial_pressures_volatiles = self.calculate_partial_pressures(fO2,a,logKr,P_volatile,volatile_comp)
        
        # Calculate mass balance eq.
        mass_balance_eq = 0
        mass_balance_eq += fO2

        # Sum O bearing volatile species
        for spec in partial_pressures_volatiles:
            if 'O2' in spec:
                mass_balance_eq += partial_pressures_volatiles[spec]
            elif 'O' in spec:
                mass_balance_eq += 0.5*partial_pressures_volatiles[spec]
        
        # Sum O bearing silicate species
        for i,gas in enumerate(self.cdefg.index):
            mass_balance_eq += self.cdefg['di'].iloc[i] * partial_pressures[i]
        
        
        return mass_balance_eq
    
    
    def calculate_partial_pressures(self,fO2,a,logKr, P_volatile, volatile_comp,final_output=False): 
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

        # Calculating total pressure per el
        tot_volatile_comp = sum(volatile_comp.values()) # ensures norm to 1
        P0_vol = {}
        for el in volatile_comp:
            P0_vol[el] = P_volatile * volatile_comp[el]/tot_volatile_comp

        # Calculate the volatile partial pressures
        partial_pressures_volatiles = pd.DataFrame()

        # H species
        partial_pressures_volatiles['H2O(g)'] = fO2**0.5*P0_vol['H']/ (2*(10**logKr['H2O(g)']+fO2**0.5))
        partial_pressures_volatiles['H2(g)'] = 0.5*(P0_vol['H'] - 2*partial_pressures_volatiles['H2O(g)'])
        
        # C species
        partial_pressures_volatiles['CO2(g)'] = fO2**0.5*P0_vol['C']/ (10**logKr['CO2(g)']+fO2**0.5)
        partial_pressures_volatiles['CO(g)'] = P0_vol['C'] - partial_pressures_volatiles['CO2(g)']
        
        # S species
        partial_pressures_volatiles['SO2(g)'] = fO2**0.5*P0_vol['S']/ (10**logKr['SO2(g)']+fO2**0.5)
        partial_pressures_volatiles['SO(g)'] = P0_vol['S'] - partial_pressures_volatiles['SO2(g)']

        # N species
        partial_pressures_volatiles['NO(g)'] = 0.5 * ((fO2/(2*10**logKr['NO(g)']**2))**0.5\
                                                * ((fO2/(2*10**logKr['NO(g)']**2) + 4*P0_vol['N']))**0.5\
                                                - fO2/(2*10**logKr['NO(g)']**2))
        partial_pressures_volatiles['N2(g)'] = 0.5*(P0_vol['N'] - partial_pressures_volatiles['NO(g)'])

        # P species
        partial_pressures_volatiles['PO(g)'] = fO2**0.5*P0_vol['P']/ (10**logKr['PO(g)']+fO2**0.5)
        partial_pressures_volatiles['P(g)'] = P0_vol['P'] - partial_pressures_volatiles['PO(g)']

        partial_pressures_vapor = pd.DataFrame()
        
        # O
        partial_pressures_vapor[0] = 10**logKr['O(g)'] * fO2**self.cdefg['di'].iloc[0]
            
        for i,spec in enumerate(self.cdefg.index[1:]):
            i += 1
            # Importing names and values from cdef file
            endmember = self.cdefg.iloc[i]['endmember']
            liq_oxide1 = self.cdefg.iloc[i]['liq_oxide1'].replace('(l)','')
            liq_oxide2 = self.cdefg.iloc[i]['liq_oxide2'].replace('(l)','')
            ci = self.cdefg.iloc[i]['ci']
            di = self.cdefg.iloc[i]['di']
            ei = self.cdefg.iloc[i]['ei'] 
            fi = self.cdefg.iloc[i]['fi']
            gi = self.cdefg.iloc[i]['gi']

            
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
            partial_pressures_vapor[i] = 10**logKr[spec]\
                                  * a_endmember**ci\
                                  * fO2**di\
                                  * partial_pressures_volatiles['H2(g)']**ei\
                                  * a_liq_oxide1**fi\
                                  * a_liq_oxide2**gi
            

        # Ensures that if this is the final output, all pressures of same 
        # species are summed. Necessary because there may be more than one 
        # reaction producing the same species (for Fe for example).
        # TODO: DEPRECATED - REMOVE

        if final_output:
            partial_pressures = pd.DataFrame()
            for i,vapor in enumerate(self.cdefg.index):
                if vapor in partial_pressures_vapor.columns:
                    partial_pressures[vapor] += partial_pressures_vapor[i]
                else:
                    partial_pressures[vapor] = partial_pressures_vapor[i]

            return partial_pressures,partial_pressures_volatiles

        else:
            return partial_pressures_vapor, partial_pressures_volatiles
        
        
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

        # Removing all non-liquid phases
        included_phases = {}
        for phase in self.melts.get_phase_inclusion_status():
            if phase != 'Liquid':
                included_phases[phase] = False
            else:
                included_phases[phase] = True
        self.melts.set_phase_inclusion_status(included_phases)

        # Importing data for liquid reactions
        fname_cdef_liq_values = 'cdef_values_liquid_reactions.csv'
        dname_cdef_liq_values = 'data/'        
        self.cdefg_liq = pd.read_csv(dname_cdef_liq_values\
                                + fname_cdef_liq_values).set_index('liq_oxide')
        
        # Oxides that are included in calculations
        self.used_oxides = ['SiO2','MgO','Al2O3','FeO','CaO','Na2O','K2O',\
                            'TiO2','Fe2O3','H2O']

        
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
            excs_dict = self.melts.get_thermo_properties_of_phase_components(xmlout,\
                                'Liquid',mode='excess')
            included_endmembers = list(excs_dict.keys())
            excs = list(excs_dict.values())            
            
            # Convert to activity
            a[t] = np.exp(excs/(self._R*t))
        
        self.a = pd.DataFrame(a,index=included_endmembers).T