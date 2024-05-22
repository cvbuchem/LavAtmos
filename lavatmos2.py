# Standard modules
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import optimize
from scipy.interpolate import interp1d
from copy import copy
import collections
import sys
import subprocess
import os

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

    def __init__(self, paths=None):
        
        
        # Constants
        self._R = 8.314
        self.Tr = 298.15

        # Points used as a smart start for fO2
        t_dep_points = {}
        t_dep_points['T'] = [2000,2500,3000,3500,4000]
        t_dep_points['fO2'] = np.log10([1e-16,1e-11,1e-5,1e-2,1e0])
        self.fO2_interp_func = interp1d(t_dep_points['T'], t_dep_points['fO2'], fill_value='extrapolate')

        
        # Importing stoichiometries for reactions
        fname_cdef_values = 'cdef_values_vapor_reactions_2.csv'
        dname_cdef_values = 'data/'
        self.cdef = pd.read_csv(dname_cdef_values+fname_cdef_values).set_index('vapor')

        fname_ml_values = 'mass_law_values.csv'
        dname_ml_values = 'data/'
        self.ml_values = pd.read_csv(dname_ml_values+fname_ml_values).set_index('species')
        
        # Importing thermo data
        self.thermo_data = janaf_data_importer() # janaf data
        self.thermo_data.update(barin_data_importer()) # barin data

        # FastChem 
        ######################################################################
        ## This needs to be changed depending on where FastChem is located! ##
        ########## Comment out whichever lines suit your usecase. ############
        ######################################################################

        # For if LavAtmos2 is being run as part of the BigPipe
        self.fastchem_dir = paths.fastchem3_dir
        self.abundances_location = paths.element_abundances3
        self.fastchem_column_names = ['Pbar','Tk','n_<tot>','n_g','mu']

        # For if LavAtmos2 is standalone
        # self.fastchem_dir = 'FastChem/' 
        # self.abundances_location = self.fastchem_dir+'input/element_abundances/' 

        #######################################################################

        # Initialise states
        self.melt = MeltState()
        
        # Allowed input oxides
        self.oxides = ['SiO2','MgO','Al2O3','FeO','CaO','Na2O','K2O','TiO2','Fe2O3']

        # Vaporised elements
        self.vaporised_elements = ['Al', 'Ca', 'Fe', 'K', 'Mg', 'Na',\
                                   'O', 'Si', 'Ti']
                            
        self.all_elements = ['C','H','He','N','O','P','S','Si','Ti',\
                             'V','Cl','K','Na','Mg','F','Ca','Fe','Al']

        self.mass_law_contribution = {'C' : 0,
                                      'H' : 0,
                                      'He': 0,
                                      'N' : 0,
                                      'O' : 0.5,
                                      'P' : 0,
                                      'S' : 0,
                                      'Si':-1,
                                      'Ti':-1,
                                      'V' : 0,
                                      'Cl': 0,
                                      'K' :-0.25,
                                      'Na':-0.25,
                                      'Mg':-0.5,
                                      'F' : 0,
                                      'Ca':-0.5,
                                      'Fe':-0.75,
                                      'Al':-0.75}
                                        


        # Solar elemental abundance
        self.solar_abundance = {
                        'e-' : 0.00,
                        'Al' : 6.45,
                        'Ar' : 0.00,
                        'C' : 8.43,
                        'Ca' : 6.34,
                        'Cl' : 5.50,
                        'Co' : 4.99,
                        'Cr' : 5.64,
                        'Cu' : 4.19,
                        'F' : 4.56,
                        'Fe' : 7.50,
                        'Ge' : 3.65,
                        'H' : 12.0,        
                        'He' : 10.93,
                        'K' : 5.03,
                        'Mg' : 7.60,
                        'Mn' : 5.43,
                        'N' : 7.83,
                        'Na' : 6.24,
                        'Ne' : 7.93,
                        'Ni' : 6.22,
                        'O' : 8.69,
                        'P' : 5.41,
                        'S' : 7.12,
                        'Si' : 7.51,
                        'Ti' : 4.95,
                        'V' : 3.93,
                        'Zn' : 4.56
                        }

    def vaporise(self, T, P_volatile, melt_comp, volatile_comp, P_melt = 0.01,\
                          fO2_initial_guess = 1e-10,\
                          verbose = True):

        self.P_volatile = P_volatile

        # Ensures that T is iterable even if just one value is given
        if type(T) != np.ndarray:
            T = np.array([T])
        
        # Calculate thermodynamic values
        self.melt_comp = self.melt.set_melt_comp(melt_comp,verbose)        
        self.melt.calculate_melt_chemical_potentials(T,P_melt,self.thermo_data)
        self.melt.calculate_melt_activities(T,P_melt)            
        self.logKr = self.logKr_calc(T,self.melt.mu0_liquid)
        
        # Calculate fO2 for all given temperatures
        fO2 = np.zeros(len(T))

        # Dev test commands below
        # fO2_old = [1e-2]
        # fastchem_partial_pressures = self.calculate_partial_pressures_fastchem(fO2_old,T[0])
        
        # print('\nCalculating vapor pressures for given temperature values')
        # print(T)
        
        # Progressbar settings
        with tqdm(total=len(T), file=sys.stdout) as pbar: 
            for i,t in enumerate(T):
                
                pbar.set_description(f'Calculated')
                
                # if t < 2400:
                #     fO2_initial_guess = 1e-12
                # else:
                #     fO2_initial_guess = 1e-8
                
                # if P_volatile < 1e1:
                #     fO2_initial_guess = 1e-8
                # else:
                #     fO2_initial_guess = 1e-4

                self.mass_balance_eq = 1e20
                best_mb_output = 1e20        
                count = 0

                # if t >= 3000 and P_volatile < 1e0:
                #     fO2_tries = [1e-6,1e-5,1e-4,1e-2,1e-1,1e0,1e1,1e2,1e3]                
                # else:
                #     fO2_tries = [1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-2,1e-1,1e0,1e1,1e2,1e3]                
    
                # fO2_tries_less = [1e-12,1e-14,1e-16,1e-18,1e-20,1e-22,1e-24,1e-26,1e-28,1e-30]

                # if P_volatile > 1e2:
                #     factor = 2
                # else:
                #     factor = 90
                # print('Factor:',factor)

                factor = 90
                fO2_tries = 10**self.fO2_interp_func(t)*np.array([1e-3,1e-2,1e-1,1e0,1e1,1e2])


                maxfev = 1000

                while np.abs(self.mass_balance_eq) > 1e-9 and count < len(fO2_tries):

                    fO2_initial_guess = fO2_tries[count]
                    print(f'\nTry #{count}')
                    print('fO2 initial guess:', fO2_initial_guess)
                    
                    # Calculate fO2
                    results_opt = optimize.fsolve(self.mass_balance_equation_fastchem,\
                                                  fO2_initial_guess,args=([t],volatile_comp),xtol=1e-9,
                                                  factor=factor, maxfev=maxfev)
                    print('Results opt:',results_opt)
                    print('Mass balance equation:',self.mass_balance_eq)

                    if np.abs(self.mass_balance_eq) < np.abs(best_mb_output)+np.abs(best_mb_output)*0.1:
                        print('Found new best solution!')
                        fO2[i] = results_opt[0]
                        best_mb_output = self.mass_balance_eq
                        O_abun_best = copy(self.O_abun)
                        count += 1

                        # if self.mass_balance_eq > 10 and count == 1:
                        #     fO2_tries = fO2_tries_less

                        # elif self.mass_balance_eq > 0:
                        #     print('Going back to prev try')
                        #     count -= 2
                        #     factor = factor/2

                        # elif count < len(fO2_tries) and fO2_tries[count] > 1e-12:
                        #     # Avoids trying values that are too low
                        #     while fO2[i]/100 > fO2_tries[count] and count < len(fO2_tries):
                        #         # print('count +1')
                        #         count += 1


                    elif np.abs(self.mass_balance_eq) > 1e-6:
                        count += 1                    
                    else:
                        count = len(fO2_tries)+1


                    # print(results_opt[0])
                
                # print(results_opt[1])
                # print(results_opt[2])
                # print(results_opt[3])
                
                pbar.update(1) # Update progress bar   
    
                # if self.is_iterable(P_boa):
                #     P_boa = P_boa[0]
                
                # fO2 = fO2[i]
                
                vapor_partial_pressures = self.vapor_partial_pressure_calc(fO2[i], [t])
                # print(vapor_partial_pressures.sum(axis=1).iloc[0])
                P_outgassed = vapor_partial_pressures.sum(axis=1).iloc[0]+fO2[i]
                P_boa = P_outgassed + self.P_volatile
                partial_pressures = self.calculate_partial_pressures_fastchem_loop([O_abun_best],[t],fO2[i],[P_boa],vapor_partial_pressures,volatile_comp)
        
        # print(partial_pressures)    
        # print('Done!')

        return partial_pressures

    def mb_output(self, fO2_array, T, P_volatile, melt_comp, volatile_comp, P_melt = 0.01,\
                          fO2_initial_guess = 1e-10,\
                          verbose = True):

        self.P_volatile = P_volatile

        # Ensures that T is iterable even if just one value is given
        if type(T) != np.ndarray:
            T = np.array([T])
        
        # Calculate thermodynamic values
        self.melt_comp = self.melt.set_melt_comp(melt_comp,verbose)        
        self.melt.calculate_melt_chemical_potentials(T,P_melt,self.thermo_data)
        self.melt.calculate_melt_activities(T,P_melt)            
        self.logKr = self.logKr_calc(T,self.melt.mu0_liquid)
        

        mb = np.zeros(len(fO2_array))     
        pps = {}

        for i,fO2 in enumerate(fO2_array):
            print('\nfO2',fO2)
            # mb[i], pps[fO2] = self.mass_balance_equation_fastchem_altered(np.array([fO2]),T,volatile_comp)
            mb[i] = self.mass_balance_equation_fastchem(np.array([fO2]),T,volatile_comp)
            print('mb',mb[i])

        return mb, pps

    def fastchem_debug(self, fO2_array, T, P_volatile, melt_comp, volatile_comp, P_melt = 0.01,\
                          fO2_initial_guess = 1e-10,\
                          verbose = True):

        self.P_volatile = P_volatile

        # Ensures that T is iterable even if just one value is given
        if type(T) != np.ndarray:
            T = np.array([T])
        
        # Calculate thermodynamic values
        self.melt_comp = self.melt.set_melt_comp(melt_comp,verbose)        
        self.melt.calculate_melt_chemical_potentials(T,P_melt,self.thermo_data)
        self.melt.calculate_melt_activities(T,P_melt)            
        self.logKr = self.logKr_calc(T,self.melt.mu0_liquid)
        

        mb = np.zeros(len(fO2_array))     
        # pps = {}

        for i,fO2 in enumerate(fO2_array):

            mb[i] = self.mass_balance_equation_fastchem(np.array([fO2]),T,volatile_comp)

        return mb

    def mass_balance_equation_fastchem(self,fO2,T, volatile_comp):

        '''
        Equation that needs to be solved for fO2.
        
        Parameters
        ----------
        input : float
            fO2 in bar.
        
        Returns
        -------
        eq1 : float
            The result of the mass balance equation.    
        '''
        
        if fO2 < 0:
            # print('fO2 less than 0! punished',fO2)
            return 1e20

        # Calculate partial pressures
        # partial_pressures = self.calculate_partial_pressures_fastchem_loop(fO2,T)

        # Partial pressures
        # Calculate vaporisation partial pressures
        vapor_partial_pressures = self.vapor_partial_pressure_calc(fO2, T)
        # for i,name in enumerate(self.cdef.index):
        #     print(name)
        #     print('vapor_pp',vapor_partial_pressures[i])

        # for i,vapor in enumerate(self.cdef.index):
        #     print(f'{vapor:>10}: {vapor_partial_pressures[i].iloc[0]:.4e}')
        
        print('fO2',fO2)
        P_outgassed = vapor_partial_pressures.sum(axis=1).iloc[0]+fO2
        P_boa = P_outgassed + self.P_volatile
        print('P_outgassed',P_outgassed)
        print('P BOA',P_boa)
        # print('O_abun init',O_abun_init)
        # O_abun_init = 1e3
        # print('############# Innter loop ###############')
        xtol = 1e-20
        # O_abun_init = fO2*1e7

        self.inner_loop_output = 1e20
        best_inner_loop_output = 1e20
        count = 0

        # O_abun_tries = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12]
        O_abun_tries = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12]


        while np.abs(self.inner_loop_output) > 1e-2 and count < len(O_abun_tries):

            O_abun_init = fO2*O_abun_tries[count]

            # print(f'\nTry #{count}')
            # print('O initial guess:', O_abun_init)
            
            results_opt = optimize.fsolve(self.inner_loop,\
                                          O_abun_init, args=(T,fO2,P_boa,vapor_partial_pressures, volatile_comp),\
                                          factor=90,xtol=xtol)
                
            print('Results opt:',results_opt)
            print('Inner loop output:',self.inner_loop_output)

            if np.abs(self.inner_loop_output) < np.abs(best_inner_loop_output)+np.abs(best_inner_loop_output)*0.1:
                print('Found new best solution!')
                self.O_abun = results_opt[0]
                best_inner_loop_output = self.inner_loop_output
                count += 1

            elif np.abs(self.inner_loop_output) > 1e-10:
                count += 1                    
            else:
                count = len(O_abun_tries)+1




        # print('O_abun final', self.O_abun)
        partial_pressures = self.calculate_partial_pressures_fastchem_loop([self.O_abun],T,fO2,P_boa,vapor_partial_pressures,volatile_comp)

        # Calculate mass balance eq.
        self.mass_balance_eq = 0
        # mass_balance_eq += fO2_new

        def el_stoichiometry(gas,el):
            index = gas.find(el)
            retgas[index+1]


        # print('\n  MASS BALANCE')#, mass_balance_eq)
        for gas in partial_pressures.columns:
            # print('\nGas:',gas)
            contrib = 0
            
            for el in self.mass_law_contribution:

                if el == gas:

                    # print('el',el,'is',gas)
                    # print('monoatomic')
                    # print('1 *',self.mass_law_contribution[el])

                    contrib += self.mass_law_contribution[el]

                elif el in gas:

                    index = gas.find(el)

                    # print('el',el,'in',gas)

                    # Check if there's a coefficient present
                    if index + len(el) < len(gas) and gas[index + len(el)].isdigit():
                        stoi = float(gas[index + len(el)])
                    
                    # Ensures inclusion of sinle elements at end of species
                    elif index + len(el) == len(gas):
                        # print('END CHECK')
                        # print(index + len(el), '=', len(gas))
                        stoi = 1 

                    # Avoid counting Si as S
                    elif gas[index + len(el)].isupper():
                        stoi = 1  # Default coefficient is 1
            
                    else: 
                        stoi = 0

                    # print(stoi,'*',self.mass_law_contribution[el])

                    contrib += self.mass_law_contribution[el]*stoi

            # print('total contribution:',contrib)
            # print(partial_pressures[gas].iloc[0])
            self.mass_balance_eq += contrib*partial_pressures[gas].iloc[0]


            # if gas in self.ml_values.index:
            #     # print('\n',gas)
            #     # print(self.ml_values['ml_value'].loc[gas])
            #     # print(partial_pressures[gas].iloc[0])
            #     self.mass_balance_eq += self.ml_values['ml_value'].loc[gas]*partial_pressures[gas].iloc[0]
            #     # print('mass balance', mass_balance_eq)
        
        print('Tried fO2:', fO2[0])
        print('MASS BALANCE EQ:',self.mass_balance_eq)
        print('\n')

        return self.mass_balance_eq

    def mass_balance_equation_fastchem_altered(self,fO2,T, volatile_comp):

        '''
        Equation that needs to be solved for fO2.
        
        Parameters
        ----------
        input : float
            fO2 in bar.
        
        Returns
        -------
        eq1 : float
            The result of the mass balance equation.    
        '''
        
        if fO2 < 0:
            # print('fO2 less than 0! punished',fO2)
            return 1e20

        # Calculate partial pressures
        # partial_pressures = self.calculate_partial_pressures_fastchem_loop(fO2,T)

        # Partial pressures
        # Calculate vaporisation partial pressures
        vapor_partial_pressures = self.vapor_partial_pressure_calc(fO2, T)
        # print('vapor_pp',vapor_partial_pressures)

        # for i,vapor in enumerate(self.cdef.index):
        #     print(f'{vapor:>10}: {vapor_partial_pressures[i].iloc[0]:.4e}')
        
        # print('fO2',fO2)
        P_outgassed = vapor_partial_pressures.sum(axis=1).iloc[0]+fO2
        P_boa = P_outgassed + self.P_volatile
        # print('P_outgassed',P_outgassed)
        # print('P BOA',P_boa)
        O_abun_init = fO2*2
        # O_abun_init = 1e3
        # print('############# Innter loop ###############')
        results_opt = optimize.fsolve(self.inner_loop,\
                                      O_abun_init, args=(T,fO2,P_boa,vapor_partial_pressures, volatile_comp),\
                                      factor=1)
        self.O_abun = results_opt[0]

        partial_pressures = self.calculate_partial_pressures_fastchem_loop([self.O_abun],T,fO2,P_boa,vapor_partial_pressures,volatile_comp)

        # Calculate mass balance eq.
        self.mass_balance_eq = 0
        # mass_balance_eq += fO2_new

        def el_stoichiometry(gas,el):
            index = gas.find(el)
            retgas[index+1]


        # print('\n  MASS BALANCE', mass_balance_eq)
        for gas in partial_pressures.columns:
            print('\nGas:',gas)
            contrib = 0
            for el in self.mass_law_contribution:
                if el == gas:
                    print('monoatomic')
                    print('1 *',self.mass_law_contribution[el])
                    contrib += self.mass_law_contribution[el]
                elif el in gas:
                    index = gas.find(el)
                    if gas[index+len(el)].isdigit():
                        print('el',el,'in',gas)
                        stoi = float(gas[index+len(el)])
                        print(stoi,'*',self.mass_law_contribution[el])
                        contrib += self.mass_law_contribution[el]*stoi

            print('total contribution:',contrib)
            print(partial_pressures[gas].iloc[0])
            self.mass_balance_eq += contrib*partial_pressures[gas].iloc[0]


            # if gas in self.ml_values.index:
            #     # print('\n',gas)
            #     # print(self.ml_values['ml_value'].loc[gas])
            #     # print(partial_pressures[gas].iloc[0])
            #     self.mass_balance_eq += self.ml_values['ml_value'].loc[gas]*partial_pressures[gas].iloc[0]
            #     # print('mass balance', mass_balance_eq)
        
        print('\nTried fO2:', fO2)
        print('MASS BALANCE EQ:',self.mass_balance_eq)

        return self.mass_balance_eq, partial_pressures

    def inner_loop(self,O_abun,T,fO2,P_boa,vapor_partial_pressures, volatile_comp):

        if O_abun < 0:
            # print('O2_abun less than 0! punished',O2_abun)
            return 1e20

        fastchem_partial_pressures = self.calculate_partial_pressures_fastchem_loop(O_abun,T,fO2,P_boa,vapor_partial_pressures, volatile_comp)

        fO2_output = fastchem_partial_pressures['O2'].iloc[0]
        dfO2 = np.log10(fO2)-np.log10(fO2_output)
        dfO2 = (fO2-fO2_output)/fO2
        # print('\n')
        # print('P_boa', P_boa)
        # print('fO2_target',fO2)
        # print('fO2_output',fO2_output)
        # print('O_abun',O_abun)
        # print('dfO2',dfO2)
        
        self.inner_loop_output = dfO2

        return dfO2

    def calculate_partial_pressures_fastchem_loop(self,O_abun,T,fO2,P_boa,vapor_partial_pressures, volatile_comp):
    
        # Calculate gas partial pressure using FastChem
        # print('\nP_boa', P_boa)
        O_abun = O_abun[0]
        # print('O abun in calulcate fastchem pressure function:', O_abun)
        self.calculate_fastchem_abundance(O_abun, vapor_partial_pressures, P_boa, volatile_comp)
        self.edit_fastchem_configs(T,P_boa)
        self.run_fastchem()

        fastchem_partial_pressures = self.read_fastchem_partial_pressures()
        
        # print(f'fO2: {fastchem_partial_pressures["O2"].iloc[0]:.6e}')

        total_H_pressure = 0
    
    
        # for gas in fastchem_partial_pressures:
        #     if 'H4' in gas:
        #         total_H_pressure += fastchem_partial_pressures[gas]*4
        #     elif 'H3' in gas:
        #         total_H_pressure += fastchem_partial_pressures[gas]*3
        #     elif 'H2' in gas:
        #         total_H_pressure += fastchem_partial_pressures[gas]*2
        #     elif 'H' in gas and not 'He' in gas:
        #         # print(gas)
        #         # print(output_pressure_var['H'][gas])
        #         total_H_pressure += fastchem_partial_pressures[gas]*1

        # print(f'Total H pressure: {total_H_pressure.iloc[0]:.6e}')

        total_pressure = fastchem_partial_pressures.sum(axis=1)

        # print('total_pressure',total_pressure)

        return fastchem_partial_pressures


    def calculate_fastchem_abundance(self, O_abun, vapor_partial_pressures, P_boa,\
                                     volatile_comp, final_output=False):

        # print('DEBUUUUUGGG!!!!!!!!!!!')
        # print(vapor_partial_pressures)
        # print(P_boa)

        if self.is_iterable(P_boa):
            P_boa = P_boa[0]

        mole_fractions = vapor_partial_pressures/P_boa
        frac = {}
        
        for el in self.vaporised_elements:
            frac[el] = 0

        for i,spec in enumerate(self.cdef.index):
            # print(spec)
            for el in self.vaporised_elements:

                if (el+'3') in spec:
                    frac[el] += 3*mole_fractions[i].iloc[0]
                    # print(f'{spec}, 3*{el}, ')

                elif (el+'2') in spec:
                    frac[el] += 2*mole_fractions[i].iloc[0]
                    # print(f'{spec}, 2*{el}, ')

                elif (el) in spec:
                    frac[el] += mole_fractions[i].iloc[0]
                    # print(f'{spec}, 1*{el}, ')
            # print(frac)

        # print('\nBEFORE fO2 add!!')
        # print(frac)
        # print(2*fO2,fO2)
        
        frac['O'] += O_abun

        # print('AFTER fO2 add!!')
        # print(frac)
        # print(self.P_volatile)

        for vol in volatile_comp:
            # print(vol, volatile_comp[vol], self.P_volatile)
            frac[vol] = self.P_volatile*volatile_comp[vol]/P_boa
            # print(frac[vol])

        '''
        Normalize elemental abundances
        '''
        fracN = {}
        # print(frac.values())
        sum_fracs = sum(frac.values())
        # print(sum_fracs)
        for el in frac:
            fracN[el] = frac[el]/sum_fracs

        verbose = False
        if verbose:
            print('Normalized fractional abundances:')
            for el in fracN:
                print(el)
                print(fracN[el])


        # Open abundance file template
        # TODO: Consider not hardcoding this (moving to paths file)
        template_name = 'element_abundances_template2.dat'
        output_name = 'element_abundances_output.dat'
        solar_abund_name = 'element_abundances_solar.dat'

        template = open(self.abundances_location+template_name, 'r')
        elab_file = template.read()
        template.close()

        # Set abundances to outgassed abundances
        for el in self.all_elements:
            if el in self.vaporised_elements or el in volatile_comp and fracN[el] != 0:
                abund = 12 + np.log10(fracN[el]*1e20)

                if self.is_iterable(abund):
                    abund = abund[0]
                
                elab_file = elab_file.replace(f'<<{el}>>', str(abund))
            else:
                # abund = self.solar_abundance[el]
                abund = 0
                elab_file = elab_file.replace(f'<<{el}>>', str(abund))
        
        # Save abundance file
        self.abundance_fname = self.abundances_location+output_name
        g = open(self.abundance_fname, 'w')
        g.write(elab_file)
        g.close()

        return P_boa

    def vapor_partial_pressure_calc(self, fO2, T):

        partial_pressures_vapor = pd.DataFrame()
        # print(partial_pressures_vapor)
        a = self.melt.a.loc[T]
        logKr = self.logKr.loc[T]
        # print(logKr)

        # O does not require oxide, hence separate
        partial_pressures_vapor[0] = 10**logKr[0] * fO2**self.cdef['di'].iloc[0]

        for i in range(1,len(self.cdef.index)):
            
            # Importing names and values from cdef file
            endmember = self.cdef.iloc[i]['endmember']
            liq_oxide1 = self.cdef.iloc[i]['liq_oxide1'].replace('(l)','')
            liq_oxide2 = self.cdef.iloc[i]['liq_oxide2'].replace('(l)','')
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
            
            # print('DEBUG')
            # print(logKr[i])
            # print(a_endmember**ci)
            # print(a_liq_oxide1**ei)
            # print(a_liq_oxide2**fi)
            # print(fO2)

            # Calculating partial pressure
            partial_pressures_vapor[i] = 10**logKr[i]\
                                  * a_endmember**ci\
                                  * fO2**di\
                                  * a_liq_oxide1**ei\
                                  * a_liq_oxide2**fi

            # print('pp', partial_pressures_vapor[i])

        return partial_pressures_vapor

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
        
        for i in range(len(self.cdef.index)):
            
            #Calculate gibbs energy of a vapor for a given T
            vapor = self.cdef.index[i]
            
            # Importing names and values from cdef file
            endmember = self.cdef.iloc[i]['endmember']
            liq_oxide1 = self.cdef.iloc[i]['liq_oxide1'].replace('(l)','')
            liq_oxide2 = self.cdef.iloc[i]['liq_oxide2'].replace('(l)','')
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

    def edit_fastchem_configs(self,T,P):
        # print(T,P)
        # print(f'Running FastChem for single point at T: {T[0]} [K] and P: {P[0]:.3e} [bar]')
         # Open parameter template file
        tp_point_path = 'input/tp_point.dat' 
        tp_file = open(self.fastchem_dir+tp_point_path, 'w')
        tp_file.write(f'P\tT\n{P[0]:.6e}\t{T[0]:.6e}')
        tp_file.close()

        '''
        Editing config file
        '''
        # print('\nEditing FastChem config')
        config_path = self.fastchem_dir+'input/config.input'
        param_path = self.fastchem_dir+'input/parameters.dat'

        # Config file
        template_path_config = self.fastchem_dir+'input/config_template.input'

        # Open parameter template file
        template = open(template_path_config)
        configurations = template.read()
        template.close()

        fastchem_config = {
            'param_path' : param_path,
            'tp_grid_path' : tp_point_path, 
            'output_abundance_fname' : 'output/boa_chem.dat',
            'element_abundance_file' : self.abundance_fname.replace('FastChem/','')
            }
        for config in fastchem_config:
            configurations = configurations.replace(f'<<{config}>>', fastchem_config[config])

        # Save config file
        config_file = open(config_path, 'w')
        config_file.write(configurations)
        config_file.close()

        '''
        Editing param file
        '''
        # print('Editing FastChem param')
        # Parameter file

        '''
        template_path_param = self.fastchem_data+'parameters_template.dat'

        # Open parameter template file
        template = open(template_path_param)
        parameters = template.read()
        template.close()

        fastchem_param = {\
            'element_abundance_file' : self.abundance_fname,
            'species_data_file' : self.fastchem_data+'logK/mixed_ion.dat'
        }

        for param in fastchem_param:
            parameters = parameters.replace(f'<<{param}>>',fastchem_param[param])

        # Save abundance file
        param_file = open(param_path, 'w')
        param_file.write(parameters)
        param_file.close()
        '''


    def run_fastchem(self):
        # Check call instead of call can catch the error 
        try: 
            subprocess.check_call([f'./fastchem input/config.input'],\
                                  shell=True,\
                                  cwd=f'{self.fastchem_dir}/',stdout=subprocess.DEVNULL) 
        except: 
            print(f'\nFastChem cannot run properly.')
            print(f'Try compile it by running make under {self.fastchem_dir}\n'); raise 

    def read_fastchem_partial_pressures(self):

        fname = self.fastchem_dir+'output/boa_chem.dat'
        fastchem_partial_pressures = pd.read_csv(fname, sep='\s+')
        return fastchem_partial_pressures.drop(columns=self.fastchem_column_names)\
               *fastchem_partial_pressures[self.fastchem_column_names[0]].iloc[0]

    def is_iterable(self,variable):
        try:
            iter(variable)
            return True
        except TypeError:
            return False


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
            excs_dict = self.melts.get_thermo_properties_of_phase_components(xmlout,\
                                'Liquid',mode='excess')
            included_endmembers = list(excs_dict.keys())
            excs = list(excs_dict.values())            
            
            # Convert to activity
            a[t] = np.exp(excs/(self._R*t))
        
        self.a = pd.DataFrame(a,index=included_endmembers).T