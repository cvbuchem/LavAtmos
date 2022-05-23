from data.databases.janaf.janaf import Janafdb


def janaf_data_importer():
    '''
    Imports Janaf data for specified species using a modified version of thermochem.
    
    Returns
    -------
    janafData : dict
        Phase data for all relevant species.
    
    '''
    db = Janafdb()

    janafData = {}

    janafData['O2(g)'] = db.getphasedata(formula='O2')
    janafData['O(g)'] = db.getphasedata(formula='O',phase='g')

    janafData['SiO2(l)'] = db.getphasedata(formula='O2Si',phase='l')
    janafData['SiO2(g)'] = db.getphasedata(formula='O2Si',phase='g')
    janafData['Si(g)'] = db.getphasedata(formula='Si',phase='g')
    janafData['Si(l)'] = db.getphasedata(formula='Si',phase='l')
    janafData['Si2(g)'] = db.getphasedata(formula='Si2',phase='g')
    janafData['Si3(g)'] = db.getphasedata(formula='Si3',phase='g')
    janafData['SiO(g)'] = db.getphasedata(formula='OSi',phase='g')

    janafData['MgO(l)'] = db.getphasedata(formula='MgO',phase='l')
    janafData['MgO(g)'] = db.getphasedata(formula='MgO',phase='g')
    janafData['Mg(g)'] = db.getphasedata(formula='Mg',phase='g')
    janafData['Mg(l)'] = db.getphasedata(formula='Mg',phase='l')
    janafData['Mg2(g)'] = db.getphasedata(formula='Mg2',phase='g')

    janafData['Al2O3(l)'] = db.getphasedata(formula='Al2O3',phase='l')
    janafData['Al(g)'] = db.getphasedata(formula='Al',phase='g')
    janafData['Al(l)'] = db.getphasedata(formula='Al',phase='l')
    janafData['Al2(g)'] = db.getphasedata(formula='Al2',phase='g')
    janafData['AlO(g)'] = db.getphasedata(formula='AlO',phase='g')
    janafData['Al2O(g)'] = db.getphasedata(formula='Al2O',phase='g')
    janafData['Al2O2(g)'] = db.getphasedata(formula='Al2O2',phase='g')

    janafData['FeO(l)'] = db.getphasedata(formula='FeO',phase='l')
    # janafData['Fe2O3(l)'] = db.getphasedata(formula='Fe2O3',phase='l') # Only available in crystal phase
    janafData['FeO(g)'] = db.getphasedata(formula='FeO',phase='g')
    janafData['Fe(g)'] = db.getphasedata(formula='Fe',phase='g')
    janafData['Fe(l)'] = db.getphasedata(formula='Fe',phase='l')

    janafData['CaO(l)'] = db.getphasedata(formula='CaO',phase='l')
    janafData['CaO(g)'] = db.getphasedata(formula='CaO',phase='g')
    janafData['Ca(g)'] = db.getphasedata(formula='Ca',phase='g')
    janafData['Ca(l)'] = db.getphasedata(formula='Ca',phase='l')
    janafData['Ca2(g)'] = db.getphasedata(formula='Ca2',phase='g')

    # janafData['K2O(l)'] = db.getphasedata(formula='K2O',phase='l') # Only available in crystal phase
    # janafData['K2O(g)'] = db.getphasedata(formula='K2O',phase='g') # NO K2O GAS?
    janafData['K(g)'] = db.getphasedata(formula='K',phase='g')
    janafData['K(l)'] = db.getphasedata(formula='K',phase='l')
    janafData['K2(g)'] = db.getphasedata(formula='K2',phase='g')
    janafData['KO(g)'] = db.getphasedata(formula='KO',phase='g')

    janafData['Na2O(l)'] = db.getphasedata(formula='Na2O',phase='l') 
    # janafData['Na2O(g)'] = db.getphasedata(formula='K2O',phase='g') # NO Na2O GAS?
    janafData['Na(g)'] = db.getphasedata(formula='Na',phase='g')
    janafData['Na(l)'] = db.getphasedata(formula='Na',phase='l')
    janafData['Na2(g)'] = db.getphasedata(formula='Na2',phase='g')
    janafData['NaO(g)'] = db.getphasedata(formula='NaO',phase='g')

    janafData['TiO2(l)'] = db.getphasedata(formula='O2Ti',phase='l')
    janafData['TiO2(g)'] = db.getphasedata(formula='O2Ti',phase='g')
    janafData['Ti(g)'] = db.getphasedata(formula='Ti',phase='g')
    janafData['Ti(l)'] = db.getphasedata(formula='Ti',phase='l')
    janafData['TiO(g)'] = db.getphasedata(formula='OTi',phase='g')

    janafData['Cr2O3(l)'] = db.getphasedata(formula='Cr2O3',phase='l')
    janafData['Cr(g)'] = db.getphasedata(formula='Cr',phase='g')
    janafData['Cr(l)'] = db.getphasedata(formula='Cr',phase='l')
    janafData['CrO(g)'] = db.getphasedata(formula='CrO',phase='g')
    janafData['CrO2(g)'] = db.getphasedata(formula='CrO2',phase='g')
    janafData['CrO3(g)'] = db.getphasedata(formula='CrO3',phase='g')
    
    return janafData