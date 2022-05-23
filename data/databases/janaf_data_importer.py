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

    logK = {}

    logK['O2(g)'] = db.getphasedata(formula='O2').logKf
    logK['O(g)'] = db.getphasedata(formula='O',phase='g').logKf

    logK['SiO2(l)'] = db.getphasedata(formula='O2Si',phase='l').logKf
    logK['SiO2(g)'] = db.getphasedata(formula='O2Si',phase='g').logKf
    logK['Si(g)'] = db.getphasedata(formula='Si',phase='g').logKf
    logK['Si(l)'] = db.getphasedata(formula='Si',phase='l').logKf
    logK['Si2(g)'] = db.getphasedata(formula='Si2',phase='g').logKf
    logK['Si3(g)'] = db.getphasedata(formula='Si3',phase='g').logKf
    logK['SiO(g)'] = db.getphasedata(formula='OSi',phase='g').logKf

    logK['MgO(l)'] = db.getphasedata(formula='MgO',phase='l').logKf
    logK['MgO(g)'] = db.getphasedata(formula='MgO',phase='g').logKf
    logK['Mg(g)'] = db.getphasedata(formula='Mg',phase='g').logKf
    logK['Mg(l)'] = db.getphasedata(formula='Mg',phase='l').logKf
    logK['Mg2(g)'] = db.getphasedata(formula='Mg2',phase='g').logKf

    logK['Al2O3(l)'] = db.getphasedata(formula='Al2O3',phase='l').logKf
    logK['Al(g)'] = db.getphasedata(formula='Al',phase='g').logKf
    logK['Al(l)'] = db.getphasedata(formula='Al',phase='l').logKf
    logK['Al2(g)'] = db.getphasedata(formula='Al2',phase='g').logKf
    logK['AlO(g)'] = db.getphasedata(formula='AlO',phase='g').logKf
    logK['Al2O(g)'] = db.getphasedata(formula='Al2O',phase='g').logKf
    logK['AlO2(g)'] = db.getphasedata(formula='AlO2',phase='g').logKf
    logK['Al2O2(g)'] = db.getphasedata(formula='Al2O2',phase='g').logKf

    logK['FeO(l)'] = db.getphasedata(formula='FeO',phase='l').logKf
    logK['FeO(g)'] = db.getphasedata(formula='FeO',phase='g').logKf
    logK['Fe(g)'] = db.getphasedata(formula='Fe',phase='g').logKf
    logK['Fe(l)'] = db.getphasedata(formula='Fe',phase='l').logKf

    logK['CaO(l)'] = db.getphasedata(formula='CaO',phase='l').logKf
    logK['CaO(g)'] = db.getphasedata(formula='CaO',phase='g').logKf
    logK['Ca(g)'] = db.getphasedata(formula='Ca',phase='g').logKf
    logK['Ca(l)'] = db.getphasedata(formula='Ca',phase='l').logKf
    logK['Ca2(g)'] = db.getphasedata(formula='Ca2',phase='g').logKf

    logK['K(g)'] = db.getphasedata(formula='K',phase='g').logKf
    logK['K(l)'] = db.getphasedata(formula='K',phase='l').logKf
    logK['K2(g)'] = db.getphasedata(formula='K2',phase='g').logKf
    logK['KO(g)'] = db.getphasedata(formula='KO',phase='g').logKf

    logK['Na2O(l)'] = db.getphasedata(formula='Na2O',phase='l').logKf 
    logK['Na(g)'] = db.getphasedata(formula='Na',phase='g').logKf
    logK['Na(l)'] = db.getphasedata(formula='Na',phase='l').logKf
    logK['Na2(g)'] = db.getphasedata(formula='Na2',phase='g').logKf
    logK['NaO(g)'] = db.getphasedata(formula='NaO',phase='g').logKf

    logK['TiO2(l)'] = db.getphasedata(formula='O2Ti',phase='l').logKf
    logK['TiO2(g)'] = db.getphasedata(formula='O2Ti',phase='g').logKf
    logK['Ti(g)'] = db.getphasedata(formula='Ti',phase='g').logKf
    logK['Ti(l)'] = db.getphasedata(formula='Ti',phase='l').logKf
    logK['TiO(g)'] = db.getphasedata(formula='OTi',phase='g').logKf

    logK['Cr2O3(l)'] = db.getphasedata(formula='Cr2O3',phase='l').logKf
    logK['Cr(g)'] = db.getphasedata(formula='Cr',phase='g').logKf
    logK['Cr(l)'] = db.getphasedata(formula='Cr',phase='l').logKf
    logK['CrO(g)'] = db.getphasedata(formula='CrO',phase='g').logKf
    logK['CrO2(g)'] = db.getphasedata(formula='CrO2',phase='g').logKf
    logK['CrO3(g)'] = db.getphasedata(formula='CrO3',phase='g').logKf
    
    return logK