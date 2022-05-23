from data.databases.janaf.janaf import Janafdb


def janaf_data_importer():
    '''
    Imports Janaf data for specified species using a modified version of thermochem.
    
    Returns
    -------
    logK : dict
        LogK data for species.
    
    '''
    db = Janafdb()

    janaf_data = {}

    # O
    janaf_data['O2(g)'] = db.getphasedata(formula='O2')
    janaf_data['O(g)'] = db.getphasedata(formula='O',phase='g')

    # Si
    janaf_data['SiO2(l)'] = db.getphasedata(formula='O2Si',phase='l') # used for test

    janaf_data['SiO2(g)'] = db.getphasedata(formula='O2Si',phase='g')
    janaf_data['Si(g)'] = db.getphasedata(formula='Si',phase='g')
    janaf_data['Si2(g)'] = db.getphasedata(formula='Si2',phase='g')
    janaf_data['Si3(g)'] = db.getphasedata(formula='Si3',phase='g')
    janaf_data['SiO(g)'] = db.getphasedata(formula='OSi',phase='g')

    # Mg
    janaf_data['MgO(l)'] = db.getphasedata(formula='MgO',phase='l')
    janaf_data['Mg2SiO4(l)'] = db.getphasedata(formula='Mg2O4Si',phase='l') # used for test

    janaf_data['MgO(g)'] = db.getphasedata(formula='MgO',phase='g')
    janaf_data['Mg(g)'] = db.getphasedata(formula='Mg',phase='g')
    janaf_data['Mg2(g)'] = db.getphasedata(formula='Mg2',phase='g')
    
    # Al
    janaf_data['Al2O3(l)'] = db.getphasedata(formula='Al2O3',phase='l')

    janaf_data['Al(g)'] = db.getphasedata(formula='Al',phase='g')
    janaf_data['Al2(g)'] = db.getphasedata(formula='Al2',phase='g')
    janaf_data['AlO(g)'] = db.getphasedata(formula='AlO',phase='g')
    janaf_data['Al2O(g)'] = db.getphasedata(formula='Al2O',phase='g')
    janaf_data['AlO2(g)'] = db.getphasedata(formula='AlO2',phase='g')
    janaf_data['Al2O2(g)'] = db.getphasedata(formula='Al2O2',phase='g')

    # Fe
    janaf_data['FeO(l)'] = db.getphasedata(formula='FeO',phase='l')
    
    janaf_data['FeO(g)'] = db.getphasedata(formula='FeO',phase='g')
    janaf_data['Fe(g)'] = db.getphasedata(formula='Fe',phase='g')
    
    # Ca 
    janaf_data['CaO(l)'] = db.getphasedata(formula='CaO',phase='l')
    
    janaf_data['CaO(g)'] = db.getphasedata(formula='CaO',phase='g')
    janaf_data['Ca(g)'] = db.getphasedata(formula='Ca',phase='g')
    janaf_data['Ca2(g)'] = db.getphasedata(formula='Ca2',phase='g')

    # K
    janaf_data['K(g)'] = db.getphasedata(formula='K',phase='g')
    janaf_data['K2(g)'] = db.getphasedata(formula='K2',phase='g')
    janaf_data['KO(g)'] = db.getphasedata(formula='KO',phase='g')

    # Na
    janaf_data['Na2O(l)'] = db.getphasedata(formula='Na2O',phase='l') 
    janaf_data['Na2O3Si(l)'] = db.getphasedata(formula='Na2O3Si',phase='l') # used for test
    
    janaf_data['Na(g)'] = db.getphasedata(formula='Na',phase='g')
    janaf_data['Na2(g)'] = db.getphasedata(formula='Na2',phase='g')
    janaf_data['NaO(g)'] = db.getphasedata(formula='NaO',phase='g')

    # Ti
    janaf_data['TiO2(l)'] = db.getphasedata(formula='O2Ti',phase='l')
    
    janaf_data['TiO2(g)'] = db.getphasedata(formula='O2Ti',phase='g')
    janaf_data['Ti(g)'] = db.getphasedata(formula='Ti',phase='g')
    janaf_data['TiO(g)'] = db.getphasedata(formula='OTi',phase='g')

    # Cr
    janaf_data['Cr2O3(l)'] = db.getphasedata(formula='Cr2O3',phase='l')

    janaf_data['Cr(g)'] = db.getphasedata(formula='Cr',phase='g')
    janaf_data['CrO(g)'] = db.getphasedata(formula='CrO',phase='g')
    janaf_data['CrO2(g)'] = db.getphasedata(formula='CrO2',phase='g')
    janaf_data['CrO3(g)'] = db.getphasedata(formula='CrO3',phase='g')
    
    return janaf_data