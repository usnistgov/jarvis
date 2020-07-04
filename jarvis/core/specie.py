"""Modules related to chemistry of periodic-table elements."""

import os
import json
import numpy as np

el_chem_json_file = str(os.path.join(
                        os.path.dirname(__file__), "Elements.json"))
el_chem_json = open(el_chem_json_file, "r")
chem_data = json.load(el_chem_json)
el_chem_json.close()

el_chrg_json_file = str(os.path.join(
                        os.path.dirname(__file__), "element_charge.json"))
el_chrg_json = open(el_chrg_json_file, "r")
chrg_data = json.load(el_chrg_json)
el_chrg_json.close()


def get_descrp_arr_name(elm="Al"):
    """
    Get chemical descriptors for an element.

    Can be used in JARVIS-ML.

    Args:

           elm: element name

    Returns:
             arr: array value
    """
    arr = []
    try:
        dat = chem_data
        d = dat[elm]
        arr = []
        for k, v in d.items():
            arr.append(k)
    except Exception:
        pass
    return arr


class Specie(object):
    """
    Specie object for chemistry information.

    Used in defining chemistry of a material.

    >>> el = Specie('Al')
    >>> el.Z
    13
    >>> round(el.atomic_mass,2)
    26.98
    >>> el.symbol
    'Al'
    >>> round(el.get_chgdescrp_arr[1],2)
    12.17
    >>> round(el.get_descrp_arr[1],2)
    2792.11
    >>> el = Specie('asdfg')
    >>> el.element_property("asdfg")
    nan

    """

    def __init__(self, symbol=""):
        """Initialize with periodic table element."""
        self.symbol = symbol
        self._data = chem_data

    @property
    def Z(self):
        """Get atomic number."""
        return self.element_property("Z")

    @property
    def atomic_mass(self):
        """Get atomic mass."""
        return self.element_property("atom_mass")

    @property
    def get_chgdescrp_arr(self):
        """
        Get charge descriptors for an element.

        Gives 378 array data.

        Args:

           elm: element name

        Returns:
             arr: array value
        """
        arr = []

        # try:
        arr = chrg_data[self.symbol][0][1]
        # except:
        #    pass
        return arr

    @property
    def get_descrp_arr(self):
        """
        Get chemical descriptors for an element.

        Gives 438 array data.

        Args:

           elm: element name

        Returns:
             arr: array value
        """
        arr = []

        d = chem_data[self.symbol]
        arr = []
        for k, v in d.items():
            arr.append(v)
        arr = np.array(arr).astype(float)
        return arr

    @property
    def atomic_rad(self):
        """Get atomic radii."""
        return self.element_property("atom_rad")

    def element_property(self, key=""):
        """
        Get element property from the list of keys.

        These 84 keys are:
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
        """
        val = np.nan
        try:
            val = self._data[self.symbol][key]
        except Exception:
            pass
        return val


"""
if __name__ == "__main__":
    el = Specie("Al")
    #print(el.get_chgdescrp_arr)
    print(len(el.get_descrp_arr))
    #print (get_descrp_arr_name())
"""
