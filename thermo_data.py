# thermo_data.py
# All thermo databases + helpers (incl. Criss–Cobble) for the Pourbaix workflow.

import numpy as np
from pymatgen.core import Element, Composition

# ============================  CONSTANTS / HELPERS  ==============================================

T_REF = 298.15  # reference T used throughout

def cp_poly(A, B=0, C=0, D=0, E=0):
    """Cp(T) [J/mol/K] = A + B T + C/T^2 + D T^3 + E/T^2."""
    def _cp(T):
        t = max(float(T), 1e-6)
        return A + B*t + C/(t*t) + D*t**3 + E/(t*t)
    return _cp

# Standard states for elements → reference species
STD_STATE = {"H": ("H2", 2), "N": ("N2", 2), "O": ("O2", 2),
             "F": ("F2", 2), "Cl": ("Cl2", 2), "Br": ("Br2", 2), "I": ("I2", 2)}
def standard_species_for_element(el: Element):
    return STD_STATE.get(el.symbol, (el.symbol, 1))

# ============================  THERMO DATABASES  =================================================

# Elemental standard states (S°298 in kJ/mol/K; Cp(T) in J/mol/K).
# >>> Replace placeholders with trusted data for your systems. <<<
element_db = {
    "Al": {"S298": 0.02832, "cp": cp_poly(24.2)},"Li": {"S298": 0.0291, "cp": cp_poly(24.6)},"Be": {"S298": 0.0165, "cp": cp_poly(16.4)},"B": {"S298": 0.0127, "cp": cp_poly(11.1)},
    "Si": {"S298": 0.01881, "cp": cp_poly(22.8)}, "C": {"S298": 0.00574, "cp": cp_poly(4.78,21.16e-3,-2.29e-5)}, "Na": {"S298": 0.0513, "cp": cp_poly(28.2)},
    "O2": {"S298": 0.20515, "cp": cp_poly(29.4, 1.0e-2)}, "Mg": {"S298": 0.0327, "cp": cp_poly(24.87)}, "P": {"S298": 0.0411, "cp": cp_poly(23.82)},
    "H2": {"S298": 0.13068, "cp": cp_poly(28.8)}, "K": {"S298": 0.0647, "cp": cp_poly(29.7)}, "Ca": {"S298": 0.0416, "cp": cp_poly(25.93)},
    "N2": {"S298": 0.19161, "cp": cp_poly(29.1, 0.5e-2)}, "Ti": {"S298": 0.0300, "cp": cp_poly(25.24)}, "V": {"S298": 0.0289, "cp": cp_poly(24.9)},
    "Cl2": {"S298": 0.22308, "cp": cp_poly(33.9)}, "Cr": {"S298": 0.0235, "cp": cp_poly(22.39,9.63e-3,-1.62e5)}, "Mn": {"S298": 0.032, "cp": cp_poly(26.3)},
    # metals (example placeholders)
    "W":  {"S298": 0.032,   "cp": cp_poly(24.0)},"Fe": {"S298": 0.0273,   "cp": cp_poly(17.71,23.77e-3,0.27e5)}, "Co": {"S298": 0.032,   "cp": cp_poly(24.74,6.47e-3,-1.37e5)},
    "Cu": {"S298": 0.0333,  "cp": cp_poly(24.36,3.99e-3,-0.99e5)}, "Ni": {"S298": 0.026,   "cp": cp_poly(25.99)}, "Zn": {"S298": 0.042,   "cp": cp_poly(25.4)},
    "Hf": {"S298": 0.0320,  "cp": cp_poly(25.0)}, "Ga": {"S298": 0.041, "cp": cp_poly(25.86)}, "Ge": {"S298": 0.031,   "cp": cp_poly(23.22)},
    "Zr": {"S298": 0.0330,  "cp": cp_poly(25.0)}, "As": {"S298": 0.036, "cp": cp_poly(23.71,4.85e-3,-0.44e-5)}, "Se": {"S298": 0.042,   "cp": cp_poly(25.0)},
    "Ti": {"S298": 0.031,  "cp": cp_poly(25.24)}, "Rb": {"S298": 0.077, "cp": cp_poly(31.06)}, "Sr": {"S298": 0.056,   "cp": cp_poly(26.79)},
    "Y":  {"S298": 0.032,   "cp": cp_poly(25.0)}, "Zr": {"S298": 0.039,   "cp": cp_poly(25.20)}, "Nb": {"S298": 0.0365,   "cp": cp_poly(24.7)},
    "Mo": {"S298": 0.0288,   "cp": cp_poly(24.0)}, "Tc": {"S298": 0.032,   "cp": cp_poly(24.0)}, "Ru": {"S298": 0.0285,   "cp": cp_poly(20.92,7.1e-3,0.93e-5)},
    "Rh": {"S298": 0.0315,   "cp": cp_poly(23.97,7.93e-3,-1.39e5)}, "Pd": {"S298": 0.0378, "cp": cp_poly(23.97,5.93e-3,-0.36e5)}, "Ag": {"S298": 0.04255,   "cp": cp_poly(22.69,6.92e-3,0.53e5)},
    "Cd": {"S298": 0.0518,  "cp": cp_poly(26.02)}, "In": {"S298": 0.042,   "cp": cp_poly(25.0)}, "Sn": {"S298": 0.051,   "cp": cp_poly(27.0)},
    "Sb": {"S298": 0.0451,   "cp": cp_poly(24.81,10.42e-3,7.87e5)}, "Te": {"S298": 0.052,   "cp": cp_poly(25.0)}, "I2": {"S298": 0.260,   "cp": cp_poly(29.4)},
    "Cs": {"S298": 0.0851,   "cp": cp_poly(32.2)}, "Ba": {"S298": 0.064,   "cp": cp_poly(28.1)}, "La": {"S298": 0.032,   "cp": cp_poly(25.0)},
    "Ta": {"S298": 0.0415,   "cp": cp_poly(25.3)}, "Hf": {"S298": 0.0436,   "cp": cp_poly(25.7)}, "W": {"S298": 0.033,   "cp": cp_poly(24.3)},
    "Re": {"S298": 0.032,   "cp": cp_poly(24.0)}, "Os": {"S298": 0.028,   "cp": cp_poly(20.92,7.1e-3,0.93e5)}, "Ir": {"S298": 0.0315,   "cp": cp_poly(23.97,7.93e-3,-1.39e5)},
    "Pt": {"S298": 0.0378,   "cp": cp_poly(23.97,5.93e-3,-0.36e5)}, "Au": {"S298": 0.04255,   "cp": cp_poly(22.69,6.92e-3,0.53e5)},
    "Hg": {"S298": 0.076,   "cp": cp_poly(27.98)}, "Tl": {"S298": 0.042,   "cp": cp_poly(25.0)}, "Pb": {"S298": 0.065,   "cp": cp_poly(27.0)},
  #  "Bi": {"S298": 0.0451,   "cp": cp_poly(24.81,10.42e-3,7.87e5)}, "Po": {"S298": 0.052,   "cp": cp_poly(25.0)}, "At2": {"S298": 0.260,   "cp": cp_poly(29.4)},
}

# Compounds (reduced formula key). Provide ΔHf°(298) [kJ/mol-fu], S°(298) [kJ/mol/K], Cp(T) callable.
compound_db = {
    "Al2O3": {"dHf298": -1675.7, "S298": 0.0509, "cp": cp_poly(79.0, 3.0e-2)}, "Li2O2": {"dHf298": -632.6, "S298": 0.0572, "cp": cp_poly(70.63)},
    "LiHO":  {"dHf298": -484.1, "S298": 0.0260,  "cp": cp_poly(49.59)}, "LiH":   {"dHf298": -90.5,  "S298": 0.020,  "cp": cp_poly(28.84)},
    "Li2O":  {"dHf298": -597.9, "S298": 0.038,  "cp": cp_poly(55.0)}, "BeO":   {"dHf298": -604.0, "S298": 0.0138,  "cp": cp_poly(25.6)},
    "Be(HO)2": {"dHf298": -902.9, "S298": 0.0536, "cp": cp_poly(65.7)}, "B2O3": {"dHf298": -1272.0, "S298": 0.054, "cp": cp_poly(62.6)},
    "SiO2":  {"dHf298":  -910.9, "S298": 0.0415, "cp": cp_poly(45.0, 2.0e-2)}, "NaH": {"dHf298": -56.5,  "S298": 0.04,  "cp": cp_poly(36.39)},
    "Al":    {"dHf298":     0.0, "S298": 0.02832,"cp": cp_poly(24.2)}, "NaO2": {"dHf298": -260.0, "S298": 0.1159, "cp": cp_poly(72.1)},
    "Si":    {"dHf298":     0.0, "S298": 0.01881,"cp": cp_poly(22.8)}, "Na2O2": {"dHf298": -513.2, "S298": 0.095,  "cp": cp_poly(89.26)},
    "H2O":   {"dHf298":  -285.83, "S298": 0.06995, "cp": cp_poly(75.3)}, "NaHO": {"dHf298": -425.6, "S298": 0.064,  "cp": cp_poly(59.53)},
    "MgO":   {"dHf298":  -601.6, "S298": 0.0269, "cp": cp_poly(37.0, 1.0e-2)}, "MgH2": {"dHf298": -76.15,  "S298": 0.031,  "cp": cp_poly(35.33)},
    "Mg(HO)2": {"dHf298": -924.5, "S298": 0.063,  "cp": cp_poly(77.2)}, "SiO2": {"dHf298": -910.9, "S298": 0.0415, "cp": cp_poly(45.0, 2.0e-2)},
    "KO2":  {"dHf298": -284.5, "S298": 0.122, "cp": cp_poly(77.5)}, "K2O": {"dHf298": -363.32, "S298": 0.094, "cp": cp_poly(83.7)},
    "K2O2": {"dHf298": -495.8, "S298": 0.113, "cp": cp_poly(100.165)}, "CaO": {"dHf298": -635.1, "S298": 0.038, "cp": cp_poly(40.0, 1.0e-2)},
    "CaH2": {"dHf298": -181.5, "S298": 0.0414, "cp": cp_poly(41)}, "ScH2": {"dHf298": -90.0, "S298": 0.032, "cp": cp_poly(30.0)},
    "TiO2": {"dHf298": -944.0, "S298": 0.0503, "cp": cp_poly(50.0, 2.0e-2)}, "TiO": {"dHf298": -542.7, "S298": 0.035, "cp": cp_poly(40.0)},
    "Ti3O5": {"dHf298": -2459.1, "S298": 0.129, "cp": cp_poly(154.8)}, "Ti2O3": {"dHf298": -1520.9, "S298": 0.077, "cp": cp_poly(95.8)},
    "TiH2": {"dHf298": -144.35, "S298": 0.032, "cp": cp_poly(30.0)}, "Ti4O7": {"dHf298": -3060.0, "S298": 0.198, "cp": cp_poly(200.0)},
    "V2O3": {"dHf298": -1120.0, "S298": 0.098, "cp": cp_poly(105)}, "V2O5": {"dHf298": -1550.0, "S298": 0.130, "cp": cp_poly(130.6)},
    "VO2":  {"dHf298": -363.0, "S298": 0.042, "cp": cp_poly(50.0)}, "Cr2O3": {"dHf298": -1139.7, "S298": 0.082, "cp": cp_poly(120.0)},
    "MnO2": {"dHf298": -520.0, "S298": 0.053, "cp": cp_poly(54.1)}, "MnO": {"dHf298": -385.2, "S298": 0.06, "cp": cp_poly(45.4)},
    "Mn2O3": {"dHf298": -958.9, "S298": 0.011, "cp": cp_poly(107.65)}, "Mn3O4": {"dHf298": -1387.8, "S298": 0.155, "cp": cp_poly(139.66)},
    "FeO": {"dHf298": -272.0, "S298": 0.061, "cp": cp_poly(50.0)}, "Fe2O3": {"dHf298": -824.2, "S298": 0.087, "cp": cp_poly(103.8)},
    "Fe3O4": {"dHf298": -1118.4, "S298": 0.145, "cp": cp_poly(120.0)}, "FeO2": {"dHf298": -416.0, "S298": 0.065, "cp": cp_poly(60.0)},
    "CoO": {"dHf298": -237.0, "S298": 0.053, "cp": cp_poly(55.3)}, "Co3O4": {"dHf298": -910, "S298": 0.114, "cp": cp_poly(123.05)},
    "CoHO2": {"dHf298": -493.0, "S298": 0.075, "cp": cp_poly(80.0)}, "NiHO2": {"dHf298": -493.0, "S298": 0.075, "cp": cp_poly(80.0)},
    #"Ni(HO)2": {"dHf298": -475.0, "S298": 0.065, "cp": cp_poly(70.0)}, "NiO": {"dHf298": -240.0, "S298": 0.055, "cp": cp_poly(55.3)},
    "CuO":  {"dHf298": -157.3, "S298": 0.042, "cp": cp_poly(42.6)}, "Cu2O": {"dHf298": -170.7, "S298": 0.092, "cp": cp_poly(62.54)},
    "Cu(HO)2": {"dHf298": -450.37, "S298": 0.108, "cp": cp_poly(95.10)}, "Cu2O3": {"dHf298": -398.0, "S298": 0.11, "cp": cp_poly(80.0)},
    "ZnO": {"dHf298": -348.0, "S298": 0.057, "cp": cp_poly(50.0)}, "Zn(HO)2": {"dHf298": -641, "S298": 0.081, "cp": cp_poly(70.0)},
    "ZnO2": {"dHf298": -490.0, "S298": 0.075, "cp": cp_poly(70.0)}, "Ga2O3": {"S298": -0.256, "cp": cp_poly(-91.08)}, "GeO2": {"dHf298": -580.0, "S298": 0.042, "cp": cp_poly(45.0)},
    "SeO2": {"dHf298": -260.0, "S298": 0.057, "cp": cp_poly(50.0)}, "RbO2": {"dHf298": -284.5, "S298": 0.122, "cp": cp_poly(77.5)}, "SrO": {"dHf298": -592.04, "S298": 0.056, "cp": cp_poly(45.4)},
    "SrO2": {"dHf298": -723.0, "S298": 0.095, "cp": cp_poly(80.0)}, "Y2O3": {"dHf298": -2700.0, "S298": 0.120, "cp": cp_poly(150.0)},
    "ZrO2": {"dHf298": -1080.0, "S298": 0.050, "cp": cp_poly(56.19)}, "ZrO": {"dHf298": -540.0, "S298": 0.045, "cp": cp_poly(50.0)},
    "NbO": {"dHf298": -419.7, "S298": 0.046, "cp": cp_poly(41.11)}, "Nb2O5": {"dHf298": -1899.5, "S298": 0.137, "cp": cp_poly(131.99)},
    "NbO2": {"dHf298": -794.96, "S298": 0.055, "cp": cp_poly(57.5)}, "MoO2": {"dHf298": -587.9, "S298": 0.046, "cp": cp_poly(55.9)},
    "MoO3": {"dHf298": -745.2, "S298": 0.077, "cp": cp_poly(74.9)}, "TcO2": {"dHf298": -540.0, "S298": 0.046, "cp": cp_poly(55.9)},
    "TcO2": {"dHf298": -540.0, "S298": 0.046, "cp": cp_poly(55.9)}, "RuO2": {"dHf298": -520.0, "S298": 0.045, "cp": cp_poly(50.0)},
    "RuO4": {"dHf298": -603.0, "S298": 0.065, "cp": cp_poly(60.0)}, "RhO2": {"dHf298": -500.0, "S298": 0.044, "cp": cp_poly(50.0)},
    "RhO2": {"dHf298": -500.0, "S298": 0.044, "cp": cp_poly(50.0)}, "PdO": {"dHf298": -200.0, "S298": 0.042, "cp": cp_poly(45.0)},
    "PdO2": {"dHf298": -400.0, "S298": 0.065, "cp": cp_poly(60.0)}, "SnO2": {"dHf298": -557.6, "S298": 0.049, "cp": cp_poly(52.59)},
    "SnO": {"dHf298": -284.0, "S298": 0.057, "cp": cp_poly(44.31)}, "Sb2O3": {"dHf298": -540.0, "S298": 0.065, "cp": cp_poly(60.0)},
    "Cs2O": {"dHf298": -345.77, "S298": 0.146, "cp": cp_poly(75.98)}, "BaO": {"dHf298": -550.0, "S298": 0.064, "cp": cp_poly(50.0)},
    "BaO2": {"dHf298": -670.0, "S298": 0.095, "cp": cp_poly(80.0)}, "Ba(HO)2": {"dHf298": -946.3, "S298": 0.107, "cp": cp_poly(101.6)},
    "HfO2": {"dHf298": -1088.0, "S298": 0.059, "cp": cp_poly(60.25)}, "Ta2O5": {"dHf298": -2045.9, "S298": 0.143, "cp": cp_poly(135.03)},
    "WO3":   {"dHf298":  -842.9,  "S298": 0.076,   "cp": cp_poly(73.3)}, "WO2":   {"dHf298":  -589.7,  "S298": 0.051,   "cp": cp_poly(55.74)},
    "W":     {"dHf298":     0.0,  "S298": 0.032,   "cp": cp_poly(24.0)}, "OsO4": {"dHf298": -394.1, "S298": 0.144, "cp": cp_poly(120.0)},
    "IrO2": {"dHf298": -274.0, "S298": 0.044, "cp": cp_poly(57.3)}, "HgO": {"dHf298": -90.0, "S298": 0.069, "cp": cp_poly(44.15)},
    "Cu":    {"dHf298":     0.0,  "S298": 0.0333,  "cp": cp_poly(24.5)},
}

# Aqueous/liquid reference species used in ion reactions
species_db = {
    "H2O(l)": {"S298": 0.06995, "cp": cp_poly(75.3)},
    "H+(aq)": {"S298": 0.0,      "cp": cp_poly(0.0)},
    "e-":     {"S298": 0.0,      "cp": cp_poly(0.0)},
}

# Ion db keyed by signature or formula; examples for Cu+ and Cu2+ (placeholders)
ion_db = {
    ('Ion','Cu', 0.0, 1.0, 0.0): {"S298": 0.041, "cp": cp_poly(58.98)}, ('Ion','Cu', 0.0, 2.0, 0.0): {"S298": -0.0971, "cp": cp_poly(-3.26)}, # Cu+
    ('Ion','Ti', 0.0, -3.0, 0.0): {"S298": -0.1825, "cp": cp_poly(8.49)}, ('Ion','VO4', 0.0, -3.0, 0.0): {"S298": -0.147, "cp": cp_poly(180.9)}, # Ti3+
    ('Ion','VHO4', 0.0, -2.0, 0.0): {"S298": 0.017, "cp": cp_poly(378.2)}, ('Ion','Cr', 0.0, 2.0, 0.0): {"S298": -0.101, "cp": cp_poly(-24.68)},
    ('Ion','Cr', 0.0, 3.0, 0.0): {"S298": -0.322, "cp": cp_poly(-84.57)}, ('Ion','CrHO', 0.0, 2.0, 0.0): {"S298": -0.193, "cp": cp_poly(71.65)},
    ('Ion','Mn', 0.0, 2.0, 0.0): {"S298": -0.068, "cp": cp_poly(-11.41)}, ('Ion','Mn', 0.0, 3.0, 0.0): {"S298": -0.309, "cp": cp_poly(-96.02)},
    ('Ion', 'MnO4', 0.0, -1.0, 0.0): {"S298": 0.195, "cp": cp_poly(-3.74)}, ('Ion', 'MnO4', 0.0, -2.0, 0.0): {"S298": 0.065, "cp": cp_poly(-265.6)},
    ('Ion','Fe', 0.0, 2.0, 0.0): {"S298": -0.101, "cp": cp_poly(-27.16)}, ('Ion','Fe', 0.0, 3.0, 0.0): {"S298": -0.278, "cp": cp_poly(-67.23)},
    ('Ion','FeHO', 0.0, 1.0, 0.0): {"S298": -0.0283, "cp": cp_poly(65.63)}, ('Ion','FeHO', 0.0, 2.0, 0.0): {"S298": -0.13, "cp": cp_poly(-28.41)},
    ('Ion', 'Co', 0.0, 2.0, 0.0): {"S298": -0.113, "cp": cp_poly(-24.68)}, ('Ion', 'Co', 0.0, 3.0, 0.0): {"S298": -0.305, "cp": cp_poly(-84.57)},
    ('Ion', 'Ga', 0.0, 3.0, 0.0): {"S298": 0.201, "cp": cp_poly(96.56)}, ('Ion', 'YHO', 0.0, 2.0, 0.0): {"S298": -0.072, "cp": cp_poly(-68.22)},
    ('Ion', 'RuO4', 0.0, -2.0, 0.0): {"S298": 0.0276, "cp": cp_poly(-188.47)}, ('Ion', 'Rh', 0.0, 2.0, 0.0): {"S298": -0.1176, "cp": cp_poly(-33.23)},
    ('Ion', 'Rh', 0.0, 3.0, 0.0): {"S298": -0.299, "cp": cp_poly(-115.42)}, ('Ion', 'HRhO', 0.0, 1.0, 0.0): {"S298": -0.053, "cp": cp_poly(100.94)},
    ('Ion', 'HRhO', 0.0, 2.0, 0.0): {"S298": -0.135, "cp": cp_poly(5.93)}, ('Ion', 'RhO', 0.0, 1.0, 0.0): {"S298": -0.078, "cp": cp_poly(-161.85)},
    ('Ion', 'RhO', 0.0, 0.0, 0.0): {"S298": -0.081, "cp": cp_poly(26.46)}, ('Ion', 'Ag', 0.0, 1.0, 0.0): {"S298": 0.073, "cp": cp_poly(33.86)},
    ('Ion', 'Ag', 0.0, 2.0, 0.0): {"S298": -0.087, "cp": cp_poly(-16.76)}, ('Ion', 'In', 0.0, 3.0, 0.0): {"S298": -0.263, "cp": cp_poly(-54.88)},
    ('Ion', 'Sn', 0.0, 2.0, 0.0): {"S298": -0.017, "cp": cp_poly(-60.25)}, ('Ion', 'SnHO2', 0.0, -1.0, 0.0): {"S298": 0.0393, "cp": cp_poly(-183.85)},
    ('Ion', 'I', 0.0, -1.0, 0.0): {"S298": 0.107, "cp": cp_poly(-112.82)}, ('Ion', 'Hg', 0.0, 2.0, 0.0): {"S298": -0.036, "cp": cp_poly(13.92)},
    ('Ion', 'Hg', 0.0, 0.0, 0.0): {"S298": -0.0066, "cp": cp_poly(408.04)},
    "LiHO": {"S298": 0.048, "cp": cp_poly(87.1)}, "Be": {"S298": -0.134, "cp": cp_poly(-17.92)}, "BeO2": {"S298": -0.172, "cp": cp_poly(-37.41)},
    "BeHO": {"S298": -0.0749, "cp": cp_poly(126.38)}, "BF4": {"S298": 0.1799, "cp": cp_poly(-19.42)}, "H4C": {"S298": 0.0878, "cp": cp_poly(247.59)},
    "HCO3": {"S298": 0.0984, "cp": cp_poly(-30.17)}, "CO3": {"S298": -0.05, "cp": cp_poly(-276.88)}, "CO2": {"S298": 0.1194, "cp": cp_poly(206.4)},
    "CO": {"S298": 0.105, "cp": cp_poly(260.84)}, "NO3": {"S298": 0.146, "cp": cp_poly(-20.0)}, "N2": {"S298": 0.0958, "cp": cp_poly(225.01)},
    "HNO2": {"S298": 0.136, "cp": cp_poly(22.3)}, "HNO3": {"S298": 0.179, "cp": cp_poly(63.49)}, "H4N": {"S298": 0.111, "cp": cp_poly(67.66)},
    "H3N": {"S298": 0.108, "cp": cp_poly(62.48)}, "NO": {"S298": 0.119, "cp": cp_poly(218.99)}, "NO2": {"S298": 0.123, "cp": cp_poly(-92.63)},
    "F2": {"S298": -0.0138, "cp": cp_poly(-106.56)}, "HF": {"S298": 0.0941, "cp": cp_poly(58.55)}, "HF2": {"S298": 0.09247, "cp": cp_poly(-133.6)},
    "Na": {"S298": 0.05845, "cp": cp_poly(39.33)}, "NaO2": {"S298": 0.075, "cp": cp_poly(39.0)}, "Na2O2": {"S298": 0.087, "cp": cp_poly(50.0)},
    "Mg": {"S298": -0.137, "cp": cp_poly(-16.02)}, "MgHO": {"S298": -0.0799, "cp": cp_poly(132.33)}, "Al": {"S298": -0.338, "cp": cp_poly(-122.65)},
    "AlO2": {"S298": -0.029, "cp": cp_poly(-33.41)}, "Si(HO)4": {"S298": 0.189, "cp": cp_poly(66.04)}, "P2O7": {"S298": -0.117, "cp": cp_poly(-553.66)},
    "PHO3": {"S298": 0.0167, "cp": cp_poly(-155.26)}, "PH3O4": {"S298": 0.162, "cp": cp_poly(90.26)}, "PH2O3": {"S298": 0.133, "cp": cp_poly(27.72)},
    "P(HO2)2": {"S298": 0.0925, "cp": cp_poly(-23.96)}, "PH3": {"S298": 0.115, "cp": cp_poly(203.58)}, "PO4": {"S298": -0.221, "cp": cp_poly(-498.31)},
    "PH2O7": {"S298": 0.046, "cp": cp_poly(-101.02)}, "H2S": {"S298": 0.126, "cp": cp_poly(173.82)}, "H2S2O3": {"S298": 0.188, "cp": cp_poly(101.9)},
    "H(SO2)2": {"S298": 0.2134, "cp": cp_poly(139.82)}, "HS": {"S298": 0.067, "cp": cp_poly(-88.64)}, "HS2O3": {"S298": 0.1276, "cp": cp_poly(19.43)},
    "HSO4": {"S298": 0.132, "cp": cp_poly(26.99)}, "HSO5": {"S298": 0.212, "cp": cp_poly(157.92)}, "HSO3": {"S298": 0.134, "cp": cp_poly(-1.25)},
    "S2O5": {"S298": 0.154, "cp": cp_poly(-190.46)}, "S2O3": {"S298": 0.067, "cp": cp_poly(-227.56)}, "SO3": {"S298": -0.0154, "cp": cp_poly(-268.8)},
    "SO4": {"S298": 0.0185, "cp": cp_poly(-254.54)}, "HClO2": {"S298": 0.188, "cp": cp_poly(101.9)}, "Cl2": {"S298": 0.057, "cp": cp_poly(-117.14)},
    "ClO2": {"S298": 0.101, "cp": cp_poly(-75.57)}, "ClO3": {"S298": 0.162, "cp": cp_poly(-47.19)}, "ClO4": {"S298": 0.184, "cp": cp_poly(-20.44)},
    "HClO": {"S298": 0.142, "cp": cp_poly(31.39)}, "K": {"S298": 0.101, "cp": cp_poly(9.1)}, "Ca": {"S298": -0.0562, "cp": cp_poly(-26.38)},
    "CaHO": {"S298": 0.028, "cp": cp_poly(7.70)}, "CaO2": {"S298": 0.045, "cp": cp_poly(60.0)}, "Ca2O2": {"S298": 0.067, "cp": cp_poly(70.0)},
    "Sc": {"S298": -0.255, "cp": cp_poly(-47.48)}, "ScHO": {"S298": -0.0657, "cp": cp_poly(-75.44)}, "VO": {"S298": 0.021, "cp": cp_poly(-275.85)},
    "VO2": {"S298": -0.042, "cp": cp_poly(-54.45)}, "VO3": {"S298": 0.061, "cp": cp_poly(-100.0)}, "VO4": {"S298": 0.080, "cp": cp_poly(-30.0)},
    "CrHO4": {"S298": 0.195, "cp": cp_poly(10.06)}, "Cr2O7": {"S298": 0.291, "cp": cp_poly(-113.05)}, "MnHO": {"S298": 0.0013, "cp": cp_poly(38.64)},
    "CoHO2": {"S298": -0.107, "cp": cp_poly(150.96)}, "Ni": {"S298": -0.148, "cp": cp_poly(-42.77)}, "NiHO": {"S298": -0.141, "cp": cp_poly(133.08)},
    "NiHO2": {"S298": -0.181, "cp": cp_poly(209.72)}, "CuO": {"S298": -0.052, "cp": cp_poly(22.61)}, "Cu02": {"S298": -0.120, "cp": cp_poly(-209.22)},
    "CuHO": {"S298": -0.010, "cp": cp_poly(69.60)}, "CuHO2": {"S298": -0.0013, "cp": cp_poly(-101.10)}, "Zn": {"S298": -0.1096, "cp": cp_poly(-17.45)},
    "ZnHO2": {"S298": -0.0232, "cp": cp_poly(-40.53)}, "ZnHO": {"S298": -0.0488, "cp": cp_poly(96.39)}, "ZnO2": {"S298": -0.167, "cp": cp_poly(-47.08)},
    "GaHO": {"S298": 0.127, "cp": cp_poly(170.42)}, "Ga(HO)2": {"S298": 0.0503, "cp": cp_poly(77.94)}, "AsO4": {"S298": -0.162, "cp": cp_poly(-463.62)},
    "AsH2O3": {"S298": 0.1297, "cp": cp_poly(191.7)}, "AsH3O4": {"S298": 0.198, "cp": cp_poly(113.15)}, "As(HO2)2": {"S298": 0.117, "cp": cp_poly(2.43)},
    "AsHO4": {"S298": -0.0019, "cp": cp_poly(-185.53)}, "GeO3": {"S298": -0.042, "cp": cp_poly(-54.45)}, "GeHO2": {"S298": 0.0013, "cp": cp_poly(38.64)},
    "HSe": {"S298": 0.079, "cp": cp_poly(-77.32)}, "SeO3": {"S298": 0.012, "cp": cp_poly(-243.91)}, "SeO4": {"S298": 0.053, "cp": cp_poly(-220.25)},
    "HSeO3": {"S298": 0.135, "cp": cp_poly(31.87)}, "H2SeO3": {"S298": 0.208, "cp": cp_poly(131.37)}, "HSeO4": {"S298": 0.149, "cp": cp_poly(55.09)},
    "BrO3": {"S298": 0.162, "cp": cp_poly(-81.91)}, "BrO": {"S298": 0.041, "cp": cp_poly(-11107)}, "Br": {"S298": 0.083, "cp": cp_poly(-121.54)},
    "HBrO": {"S298": 0.142, "cp": cp_poly(31.39)}, "HBrO3": {"S298": 0.2134, "cp": cp_poly(139.82)}, "Rb": {"S298": 0.120, "cp": cp_poly(-11.91)},
    "Sr": {"S298": -0.0315, "cp": cp_poly(-37.39)}, "SrHO": {"S298": 0.061, "cp": cp_poly(-30.44)}, "Y": {"S298": -0.251, "cp": cp_poly(-43.77)},
    "ZrO": {"S298": -0.223, "cp": cp_poly(8.05)}, "NbO3": {"S298": 0.013, "cp": cp_poly(-165.47)}, "MnO4": {"S298": 0.0377, "cp": cp_poly(-186.12)},
    "TcO4": {"S298": 0.1987, "cp": cp_poly(-2.65)}, "Pd": {"S298": -0.088, "cp": cp_poly(-17.37)}, "AgHO2": {"S298": 0.065, "cp": cp_poly(-93.0)},
    "Cd": {"S298": -0.0728, "cp": cp_poly(-9.51)}, "CdO2": {"S298": -0.0849, "cp": cp_poly(-236.19)},"CdHO2": {"S298": -0.0456, "cp": cp_poly(11.13)},
    "CdHO": {"S298": -0.012, "cp": cp_poly(53.92)}, "InO2": {"S298": 0.066, "cp": cp_poly(-275.9)}, "InHO": {"S298": -0.187, "cp": cp_poly(64.86)},
    "SnHO": {"S298": 0.0803, "cp": cp_poly(-52.91)}, "SnO2": {"S298": -0.042, "cp": cp_poly(-54.45)}, "HIO3": {"S298": 0.166, "cp": cp_poly(69.31)},
    "IO3": {"S298": 0.118, "cp": cp_poly(-62.86)}, "HIO": {"S298": 0.095, "cp": cp_poly(-38.7)}, "IO": {"S298": -0.0058, "cp": cp_poly(-137.99)},
    "Cs": {"S298": 0.133, "cp": cp_poly(-25.89)}, "Ba": {"S298": 0.0096, "cp": cp_poly(-47.36)}, "BaHO": {"S298": 0.115, "cp": cp_poly(-93.16)},
    "Hf": {"S298": -0.465, "cp": cp_poly(1.59)}, "HfHO2": {"S298": -0.213, "cp": cp_poly(391.62)}, "HfHO": {"S298": -0.305, "cp": cp_poly(202.61)},
    "HfHO3": {"S298": -0.144, "cp": cp_poly(235.59)}, "HfO2": {"S298": -0.189, "cp": cp_poly(243.19)}, "WO4": {"S298": 0.0405, "cp": cp_poly(-173.61)},
    "ReO4": {"S298": 0.201, "cp": cp_poly(-12.01)}, "Pt": {"S298": -0.099, "cp": cp_poly(-22.77)}, "Au": {"S298": -0.229, "cp": cp_poly(-24.42)},
    "HgHO": {"S298": 0.069, "cp": cp_poly(195.79)}, "Hg(HO)2": {"S298": 0.156, "cp": cp_poly(71.32)}, 


    



}

# Charge map for oxyanions (used by Criss–Cobble fallback)
ion_charge = {
    "WO4": 2, "HWO4": 1, "W2O7": 2, "B4O7": 2, "BH2O3": 1, "B(HO)3": 0, "BO2": 1, "B(HO)4": 1, "BH4O5": 1, "B4HO7": 1, "B4H2O7": 0, "H5C2O": 1,
    "SiO4": 4, "Si(HO)4": 0, "H6C2O": 0, "HCO2": 1, "H3(CO)2": 1, "C2O4": 2, "H(CO2)2": 1, "H4CO": 0, "H5N2": 1, "H2N": 0, "TiO2": 2, "TiHO3": 1,
    "AlHO": 2, "Al(HO)4": 1, "SiH2O3": 0, "P2H3O7": 1, "P2HO7": 3, "P2H2O7": 2, "PHO4": 2, "PH4": 1, "S5O6": 2, "H2SO4": 0, "TiO": 2, "VO4": 1, "V2O7": 4,
    "VHO4": 0, "VH2O3": 3, "V2HO7": 3, "V10HO28": 5, "V(HO2)2": 1, "V2H3O7": 1, "V5HO14": 4, "Mn(HO)3": 1, "FeHO": 4, "Fe(HO)2": 1, "FeHO2": 1,
    "FeO4": 2, "FeO2": 1, "FeO2": 2, "Fe(HO)3": 0, "Co(HO)2": 0, "Ni(HO)2": 0, "Ni(HO)3": 1, "NiO": 0, "Ni(HO)4": 2, "NiO2": 2, "Zn(HO)2": 0,"GaHO3": 2,
    "GaH2O3": 1, "GaO3": 3, "GeHO2": 1, "GeHO3": 1, "AsHO2": 0, "AsO":1, "AsO2": 1, "H2SeO4": 0, "H2Se": 0, "YHO": 4, "Nb(HO)5": 0, "Nb(HO)4": 1,
    "TcHO4": 0, "H2RuO5": 0, "H2RuO2": 2, "HRuO5": 1, "AgHO": 0, "Cd(HO)2": 0, "In(HO)2": 1, "SnO3": 2, "SnHO2": 1, "SbO": 1, "SbO2": 1, "SbHO2": 0,
    "HI2O": 1, "H2IO": 1, "TaO2": 1, "IrO4": 2, "AuO3": 3, "HAuO3": 2, "H2AuO3": 1, "H3AuO3": 0, "HgHO2": 1,
}

# Safe elemental lookup (prevents KeyError & injects a placeholder on the fly)
def get_elem_S298_and_cp(symbol: str):
    if symbol in element_db:
        return element_db[symbol]["S298"], element_db[symbol]["cp"]
    print(f"[WARN] element_db missing '{symbol}'. Using placeholder S298=0.030 kJ/mol/K, Cp≈25 J/mol/K.")
    element_db[symbol] = {"S298": 0.030, "cp": cp_poly(25.0)}
    return element_db[symbol]["S298"], element_db[symbol]["cp"]

# ============================  CRISS–COBBLE DATA & HELPERS  ======================================

# cal → kJ conversion for entropies/heat capacities (per mol·K)
CAL_TO_KJ = 4.184e-3

# α(t), β(t) tables vs temperature (°C); linear interpolation inside the range
_CC_TABLE = {
    "oxyanion": {
        "T_C":  [60,   100,  150,  200],
        "alpha":[-127, -138, -133, -145],   # cal/mol/K
        "beta": [ 1.96,  2.24, 2.27, 2.53], # dimensionless
    },
    "acid_oxyanion": {
        "T_C":  [60,   100,  150,  200],
        "alpha":[-122, -135, -143, -152],   # cal/mol/K
        "beta": [ 3.44,  3.97, 3.95, 4.24], # dimensionless
    },
}

def _cc_alpha_beta(category: str, T2_K: float):
    """Return (alpha_cal_per_molK, beta) by linear interpolation in °C."""
    T2_C = T2_K - 273.15
    tbl = _CC_TABLE[category]
    alpha = float(np.interp(T2_C, tbl["T_C"], tbl["alpha"]))
    beta  = float(np.interp(T2_C, tbl["T_C"], tbl["beta"]))
    return alpha, beta

def _cc_S298_cal(Z: int, n_oxygen: int):
    """Criss–Cobble oxyanion entropy estimate at 298 K in cal/mol/K."""
    return 182.0 - 195.0 * (Z - 0.28 * n_oxygen)
