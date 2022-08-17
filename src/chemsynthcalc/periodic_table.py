periodic_table = [
    ('H', 1.0080, 2.2, 'H2O'), ('He', 4.0026, 1.2, 'He'),
    ('Li', 6.94, 0.98, 'Li2O'), ('Be', 9.0122, 1.57, 'BeO'), ('B', 10.81, 2.04, 'B2O3'), ('C', 12.011, 2.55, 'CO2'), ('N', 14.007, 3.04, 'NO2'), ('O', 15.999, 3.44, 'O'), ('F', 18.998, 3.98, 'F2O'), ('Ne', 20.180, 1.3, 'Ne'), 
    ('Na', 22.990, 0.93, 'Na2O'), ('Mg', 24.305, 1.31, 'MgO'), ('Al', 26.982, 1.61, 'Al2O3'), ('Si', 28.085, 1.9, 'SiO2'), ('P', 30.974, 2.19, 'P2O3'), ('S', 32.06, 2.58, 'SO3'), ('Cl', 35.45, 3.16, 'ClO2'), ('Ar', 39.95, 2.2, 'Ar'), 
    ('K', 39.098, 0.82, 'K2O'), ('Ca', 40.078, 1.0, 'CaO'), ('Sc', 44.956, 1.36, 'Sc2O3'), ('Ti', 47.867, 1.54, 'TiO2'), ('V', 50.942, 1.63, 'V2O5'), ('Cr', 51.996, 1.66, 'Cr2O3'), ('Mn', 54.938, 1.55, 'MnO2'), ('Fe', 55.845, 1.83, 'Fe2O3'), ('Co', 58.933, 1.88, 'Co2O3'), ('Ni', 58.693, 1.91, 'NiO'), ('Cu', 63.546, 1.9, 'Cu2O'), ('Zn', 65.38, 1.65, 'ZnO'), ('Ga', 69.723, 1.81, 'Ga2O3'), ('Ge', 72.63, 2.01, 'GeO2'), ('As', 74.922, 2.18, 'As2O3'), ('Se', 78.971, 2.55, 'Se3O4'), ('Br', 79.904, 2.96, 'BrO2'), ('Kr', 83.798, 3.0, 'Kr'), 
    ('Rb', 85.468, 0.82, 'Rb2O'), ('Sr', 87.62, 0.95, 'SrO'), ('Y', 88.906, 1.22, 'Y2O3'), ('Zr', 91.224, 1.33, 'ZrO2'), ('Nb', 92.906, 1.6, 'Nb2O5'), ('Mo', 95.95, 2.16, 'MoO3'), ('Tc', 98.906, 1.9, 'TcO2'), ('Ru', 101.07, 2.2, 'RuO2'), ('Rh', 102.91, 2.28, 'Rh2O3'), ('Pd', 106.42, 2.2, 'PdO'), ('Ag', 107.87, 1.93, 'Ag2O'), ('Cd', 112.41, 1.69, 'CdO'), ('In', 114.82, 1.78, 'In2O3'), ('Sn', 118.71, 1.96, 'SnO2'), ('Sb', 121.76, 2.05, 'Sb2O3'), ('Te', 127.60, 2.1, 'TeO3'), ('I', 126.90, 2.66, 'I2O5'), ('Xe', 131.29, 2.6, 'Xe'), 
    ('Cs', 132.91, 0.79, 'Cs2O'), ('Ba', 137.33, 0.89, 'BaO'), ('La', 138.91, 1.1, 'La2O3'), ('Ce', 140.12, 1.12, 'CeO2'), ('Pr', 140.91, 1.13, 'Pr2O3'), ('Nd', 144.24, 1.14, 'Nd2O3'), ('Pm', 146.92, 1.13, 'Pm2O3'), ('Sm', 150.36, 1.17, 'Sm2O3'), ('Eu', 151.96, 1.2, 'Eu2O3'), ('Gd', 157.25, 1.2, 'Gd2O3'), ('Tb', 158.93, 1.1, 'Tb2O3'), ('Dy', 162.5, 1.22, 'Dy2O3'), ('Ho', 164.93, 1.23, 'Ho2O3'), ('Er', 167.26, 1.24, 'Er2O3'), ('Tm', 168.93, 1.25, 'Tm2O3'), ('Yb', 173.05, 1.1, 'Yb2O3'), ('Lu', 174.97, 1.27, 'Lu2O3'), ('Hf', 178.49, 1.3, 'HfO2'), ('Ta', 180.95, 1.5, 'Ta2O5'), ('W', 183.84, 2.36, 'WO3'), ('Re', 186.21, 1.9, 'Re2O7'), ('Os', 190.23, 2.2, 'OsO3'), ('Ir', 192.22, 2.2, 'Ir2O3'), ('Pt', 195.08, 2.28, 'PtO'), ('Au', 196.97, 2.54, 'Au2O3'), ('Hg', 200.59, 2.0, 'HgO2'), ('Tl', 204.38, 1.62, 'Tl2O'), ('Pb', 207.2, 2.33, 'PbO2'), ('Bi', 208.98, 2.02, 'Bi2O3'), ('Po', 208.98, 2.0, 'PoO2'), ('At', 209.99, 2.2, 'At2O'), ('Rn', 222.02, 2.2, 'Rn'), 
    ('Fr', 223.02, 0.79, 'Fr2O'), ('Ra', 226.03, 0.9, 'RaO'), ('Ac', 227.03, 1.1, 'Ac2O3'), ('Th', 232.04, 1.3, 'ThO2'), ('Pa', 231.04, 1.5, 'Pa2O5'), ('U', 238.03, 1.38, 'UO2'), ('Np', 237.05, 1.36, 'NpO2'), ('Pu', 244.06, 1.28, 'PuO2'), ('Am', 243.06, 1.13, 'AmO2'), ('Cm', 247.07, 1.28, 'Cm2O3'), ('Bk', 247.07, 1.3, 'BkO2'), ('Cf', 251.08, 1.3, 'Cf2O3'), ('Es', 252.08, 1.3, 'Es2O3'), ('Fm', 257.10, 1.3, 'Fm2O3'), ('Md', 258.10, 1.3, 'Md2O3'), ('No', 259.10, 1.3, 'No2O3'), ('Lr', 262, 1.3, 'Lr2O3'), ('Rf', 267, 1.3, 'RfO2'), ('Db', 268, 1.3, 'Db2O5'), ('Hs', 269, 1.3, 'SgO3'), ('Bh', 270, 1.3, 'Bh2O7'), ('Sg', 271, 1.3, 'HsO4'), ('Mt', 278, 1.3, 'Mt2O3'), ('Ds', 281, 1.3, 'DsO2'), ('Rg', 281, 1.3, 'RgO'), ('Nh', 284, 1.3, 'CnO2'), ('Cn', 285, 1.3, 'Nh2O3'), ('Fl', 289, 1.3, 'FlO2'), ('Mc', 289, 1.3, 'Mc2O5'), ('Lv', 292, 1.3, 'LvO3'), ('Ts', 294, 1.3, 'Ts2O7'), ('Og', 294, 1.3, 'Og')
                                    ]
'''
Periodic table of elements in the form of (element symbol, atomic weight, electronegativity, oxide formula)
Abridged standard atomic weights are taken from https://doi.org/10.1515/pac-2019-0603. Weights 
of radioactive elements are taken from https://iupac.qmul.ac.uk/AtWt/. Pauling 
electronegativities are taken from https://www.webelements.com/periodicity/eneg_pauling/. Pauling
electronegativities for elements with unknown electronegativities (Bk-Og) is set to 1.3
'''