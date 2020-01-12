import os, json
import numpy as np

el_chem_json_file = str(os.path.join(os.path.dirname(__file__), "Elements.json"))
el_chem_json = open(el_chem_json_file, "r")
data = data = json.load(el_chem_json)
el_chem_json.close()


class Specie(object):
    """
    >>> el = Specie('Al')
    >>> el.Z
    13
    >>> round(el.atomic_mass,2)
    26.98
    >>> el.symbol
    'Al'
    >>> el = Specie('asdfg')
    >>> el.element_property("asdfg")
    nan
    """

    def __init__(self, symbol: str = ""):
        self.symbol = symbol
        self._data = data

    @property
    def Z(self):
        return self.element_property("Z")

    @property
    def atomic_mass(self):
        return self.element_property("atom_mass")

    def element_property(self, key: str = ""):
        val = np.nan
        try:
            keys = [
                "is_halogen",
                "row",
                "GV",
                "nfunfill",
                "C-9",
                "C-8",
                "C-7",
                "C-6",
                "C-5",
                "C-4",
                "C-3",
                "C-2",
                "C-1",
                "C-0",
                "me1",
                "me3",
                "me2",
                "max_oxid_s",
                "npvalence",
                "mp",
                "first_ion_en",
                "ndunfill",
                "op_eg",
                "jv_enp",
                "nfvalence",
                "polzbl",
                "oq_bg",
                "atom_rad",
                "atom_mass",
                "is_alkali",
                "C-13",
                "C-12",
                "C-11",
                "C-10",
                "C-17",
                "C-16",
                "C-15",
                "C-14",
                "C-19",
                "C-18",
                "voro_coord",
                "is_noble_gas",
                "e1",
                "e3",
                "e2",
                "is_lanthanoid",
                "ndvalence",
                "KV",
                "min_oxid_s",
                "nsunfill",
                "C-26",
                "X",
                "is_actinoid",
                "C-28",
                "C-29",
                "C-27",
                "C-24",
                "C-25",
                "C-22",
                "C-23",
                "C-20",
                "C-21",
                "avg_ion_rad",
                "nsvalence",
                "is_metalloid",
                "elec_aff",
                "coulmn",
                "mol_vol",
                "bp",
                "C-31",
                "C-30",
                "C-33",
                "C-32",
                "C-35",
                "C-34",
                "is_transition_metal",
                "block",
                "therm_cond",
                "Z",
                "is_alkaline",
                "npunfill",
                "oq_enp",
                "mop_eg",
                "hfus",
            ]
            val = self._data[self.symbol][key]
        except:
            pass
        return val
