from typing import NamedTuple


class Atom(NamedTuple):
    """
    Named tuple for representing atomic properties:
    atomic weight and type of oxide that will be used by default.
    """

    atomic_weight: float
    default_oxide: str


_PERIODIC_TABLE_DATA: dict[str, Atom] = {
    "H": Atom(1.008, "H2O"),
    "He": Atom(4.002602, "He"),
    "Li": Atom(6.94, "Li2O"),
    "Be": Atom(9.0121831, "BeO"),
    "B": Atom(10.81, "B2O3"),
    "C": Atom(12.011, "CO2"),
    "N": Atom(14.007, "NO2"),
    "O": Atom(15.999, "O"),
    "F": Atom(18.998403162, "F2O"),
    "Ne": Atom(20.1797, "Ne"),
    "Na": Atom(22.98976928, "Na2O"),
    "Mg": Atom(24.305, "MgO"),
    "Al": Atom(26.9815384, "Al2O3"),
    "Si": Atom(28.085, "SiO2"),
    "P": Atom(30.973761998, "P2O3"),
    "S": Atom(32.06, "SO3"),
    "Cl": Atom(35.45, "ClO2"),
    "Ar": Atom(39.95, "Ar"),
    "K": Atom(39.098, "K2O"),
    "Ca": Atom(40.078, "CaO"),
    "Sc": Atom(44.955907, "Sc2O3"),
    "Ti": Atom(47.867, "TiO2"),
    "V": Atom(50.9415, "V2O5"),
    "Cr": Atom(51.9961, "Cr2O3"),
    "Mn": Atom(54.938043, "MnO2"),
    "Fe": Atom(55.845, "Fe2O3"),
    "Co": Atom(58.933194, "Co2O3"),
    "Ni": Atom(58.6934, "NiO"),
    "Cu": Atom(63.546, "Cu2O"),
    "Zn": Atom(65.38, "ZnO"),
    "Ga": Atom(69.723, "Ga2O3"),
    "Ge": Atom(72.63, "GeO2"),
    "As": Atom(74.921595, "As2O3"),
    "Se": Atom(78.971, "Se3O4"),
    "Br": Atom(79.904, "BrO2"),
    "Kr": Atom(83.798, "Kr"),
    "Rb": Atom(85.4678, "Rb2O"),
    "Sr": Atom(87.62, "SrO"),
    "Y": Atom(88.905838, "Y2O3"),
    "Zr": Atom(91.222, "ZrO2"),
    "Nb": Atom(92.90637, "Nb2O5"),
    "Mo": Atom(95.95, "MoO3"),
    "Tc": Atom(97, "TcO2"),
    "Ru": Atom(101.07, "RuO2"),
    "Rh": Atom(102.90549, "Rh2O3"),
    "Pd": Atom(106.42, "PdO"),
    "Ag": Atom(107.8682, "Ag2O"),
    "Cd": Atom(112.414, "CdO"),
    "In": Atom(114.818, "In2O3"),
    "Sn": Atom(118.71, "SnO2"),
    "Sb": Atom(121.76, "Sb2O3"),
    "Te": Atom(127.6, "TeO3"),
    "I": Atom(126.90447, "I2O5"),
    "Xe": Atom(131.29, "Xe"),
    "Cs": Atom(132.90545196, "Cs2O"),
    "Ba": Atom(137.327, "BaO"),
    "La": Atom(138.90547, "La2O3"),
    "Ce": Atom(140.116, "CeO2"),
    "Pr": Atom(140.90766, "Pr2O3"),
    "Nd": Atom(144.242, "Nd2O3"),
    "Pm": Atom(145, "Pm2O3"),
    "Sm": Atom(150.36, "Sm2O3"),
    "Eu": Atom(151.964, "Eu2O3"),
    "Gd": Atom(157.249, "Gd2O3"),
    "Tb": Atom(158.925354, "Tb2O3"),
    "Dy": Atom(162.5, "Dy2O3"),
    "Ho": Atom(164.930329, "Ho2O3"),
    "Er": Atom(167.259, "Er2O3"),
    "Tm": Atom(168.934219, "Tm2O3"),
    "Yb": Atom(173.045, "Yb2O3"),
    "Lu": Atom(174.96669, "Lu2O3"),
    "Hf": Atom(178.486, "HfO2"),
    "Ta": Atom(180.94788, "Ta2O5"),
    "W": Atom(183.84, "WO3"),
    "Re": Atom(186.207, "Re2O7"),
    "Os": Atom(190.23, "OsO3"),
    "Ir": Atom(192.217, "Ir2O3"),
    "Pt": Atom(195.084, "PtO"),
    "Au": Atom(196.966570, "Au2O3"),
    "Hg": Atom(200.592, "HgO2"),
    "Tl": Atom(204.38, "Tl2O"),
    "Pb": Atom(207.2, "PbO2"),
    "Bi": Atom(208.98040, "Bi2O3"),
    "Po": Atom(209, "PoO2"),
    "At": Atom(210, "At2O"),
    "Rn": Atom(222, "Rn"),
    "Fr": Atom(223, "Fr2O"),
    "Ra": Atom(226, "RaO"),
    "Ac": Atom(227, "Ac2O3"),
    "Th": Atom(232.0377, "ThO2"),
    "Pa": Atom(231.03588, "Pa2O5"),
    "U": Atom(238.02891, "UO2"),
    "Np": Atom(237, "NpO2"),
    "Pu": Atom(244, "PuO2"),
    "Am": Atom(243, "AmO2"),
    "Cm": Atom(247, "Cm2O3"),
    "Bk": Atom(247, "BkO2"),
    "Cf": Atom(251, "Cf2O3"),
    "Es": Atom(252, "Es2O3"),
    "Fm": Atom(257, "Fm2O3"),
    "Md": Atom(258, "Md2O3"),
    "No": Atom(259, "No2O3"),
    "Lr": Atom(262, "Lr2O3"),
    "Rf": Atom(267, "RfO2"),
    "Db": Atom(270, "Db2O5"),
    "Sg": Atom(269, "SgO4"),
    "Bh": Atom(270, "Bh2O7"),
    "Hs": Atom(270, "HsO3"),
    "Mt": Atom(278, "Mt2O3"),
    "Ds": Atom(281, "DsO2"),
    "Rg": Atom(281, "RgO"),
    "Cn": Atom(285, "Cn2O3"),
    "Nh": Atom(286, "NhO2"),
    "Fl": Atom(289, "FlO2"),
    "Mc": Atom(289, "Mc2O5"),
    "Lv": Atom(293, "LvO3"),
    "Ts": Atom(293, "Ts2O7"),
    "Og": Atom(294, "Og"),
}

_ATOMS: set[str] = set(_PERIODIC_TABLE_DATA.keys())


class PeriodicTable:
    """
    Periodic table of elements in the form of "Atom symbol": Atom NamedTuple.
    The standard atomic weights are taken from [IUPAC](https://iupac.qmul.ac.uk/AtWt/).

    Note: For performance, prefer using the module-level constants
    `PERIODIC_TABLE` and `ATOMS` instead of creating instances.
    """

    def __init__(self) -> None:
        self.p_table: dict[str, Atom] = _PERIODIC_TABLE_DATA
        self.atoms: set[str] = _ATOMS


PERIODIC_TABLE: dict[str, Atom] = _PERIODIC_TABLE_DATA
ATOMS: set[str] = _ATOMS
