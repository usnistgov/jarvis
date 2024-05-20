"""Modules related to chemistry of periodic-table elements."""

import os
import json
import numpy as np
import functools
from jarvis.core.utils import digitize_array
from collections import defaultdict
from collections.abc import Iterable
from jarvis.db.jsonutils import loadjson

el_chem_json_file = str(
    os.path.join(os.path.dirname(__file__), "Elements.json")
)
el_chem_json = open(el_chem_json_file, "r")
chem_data = json.load(el_chem_json)
el_chem_json.close()

chem_data_magpie = loadjson(
    os.path.join(os.path.dirname(__file__), "magpie.json")
)

el_chrg_json_file = str(
    os.path.join(os.path.dirname(__file__), "element_charge.json")
)
el_chrg_json = open(el_chrg_json_file, "r")
chrg_data = json.load(el_chrg_json)
el_chrg_json.close()
cgcnn_feature_json = os.path.join(os.path.dirname(__file__), "atom_init.json")

element_full_name = os.path.join(
    os.path.dirname(__file__), "element_names.json"
)


def get_element_full_names():
    return loadjson(element_full_name)


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

    def __init__(self, symbol="", source="cfid"):
        """Initialize with periodic table element."""
        self.symbol = symbol
        if source == "cfid":
            # Cite reference:
            # https://doi.org/10.1103/PhysRevMaterials.2.083801
            self._data = chem_data
        elif source == "magpie":
            # Cite reference:
            # https://doi.org/10.1038/npjcompumats.2016.28
            self._data = chem_data_magpie
        else:
            raise ValueError("Option not available.", source)

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
        # print ('self._data')
        d = self._data[self.symbol]
        # d = chem_data[self.symbol]
        arr = []
        for k, v in d.items():
            arr.append(v)
        arr = np.array(arr).astype(float)
        return arr

    @property
    def atomic_rad(self):
        """Get atomic radii."""
        return self.element_property("atom_rad")

    @property
    def X(self):
        """Get electronegativity."""
        return self.element_property("X")

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


BASIC_FEATURES = [
    "Z",
    "coulmn",
    "row",
    "X",
    "atom_rad",
    "nsvalence",
    "npvalence",
    "ndvalence",
    "nfvalence",
    "first_ion_en",
    "elec_aff",
]


@functools.lru_cache(maxsize=None)
def get_node_attributes(species, atom_features="atomic_number"):
    """Get specific node features for an element."""
    feature_sets = ("atomic_number", "basic", "cfid", "cgcnn")

    if isinstance(atom_features, str):
        if atom_features not in feature_sets:
            raise NotImplementedError(
                f"atom features must be one of {feature_sets}"
            )
    elif isinstance(atom_features, Iterable):
        # allow custom list of features
        for prop in atom_features:
            if prop not in keys:
                raise NotImplementedError(
                    f"{prop} not supported in custom atom feature list"
                )
        return [
            Specie(species).element_property(prop) for prop in atom_features
        ]

    if atom_features == "cfid":
        return Specie(species).get_descrp_arr
    elif atom_features == "atomic_number":
        return [Specie(species).element_property("Z")]
    elif atom_features == "basic":
        return [
            Specie(species).element_property(prop) for prop in BASIC_FEATURES
        ]
    elif atom_features == "cgcnn":
        # load from json, key by atomic number
        key = str(Specie(species).element_property("Z"))
        with open(cgcnn_feature_json, "r") as f:
            # For alternative features use
            # get_digitized_feats_hot_encoded()
            i = json.load(f)
        try:
            return i[key]
        except KeyError:
            print(f"warning: could not load CGCNN features for {key}")
            print("Setting it to max atomic number available here, 103")
            # TODO Check for the error in oqmd_3d_no_cfid dataset
            # return i['Lr']
            return i["100"]


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


def get_specie_data():
    """Get the json and key data from Specie."""
    return keys, chem_data, chrg_data


def get_digitized_feats_hot_encoded(
    feature_names=keys, filename="feats_encoded.json"
):
    """Get OneHotEncoded features with digitized features."""
    from sklearn.preprocessing import OneHotEncoder
    import pandas as pd

    encoder = OneHotEncoder(categories="auto", sparse_output=False)
    # encoder = OneHotEncoder(categories="auto", sparse=False)
    dat = defaultdict()
    for i, j in chem_data.items():
        tmp = defaultdict()
        for r, s in j.items():
            if r in feature_names:
                tmp[r] = s
        dat[Specie(i).Z] = tmp  # j.values()
    df = pd.DataFrame(dat)
    df = df.T.replace(-9999.0, 0).replace(-0.0, 0).astype("float")

    for i in df.columns:
        df[i] = digitize_array(df[i])
    df = df.T

    vals = []
    for i in range(len(df.values)):
        output = encoder.fit_transform(
            np.array(df.values[i], dtype="float").reshape(-1, 1)
        )  # .toarray()
        vals.extend(output.T)
    vals = np.array(vals, dtype="float").T
    cols = df.columns.tolist()
    new_dat = {}
    for i, j in zip(cols, vals):
        new_dat[int(i)] = list([int(m) for m in j])
    if filename is not None:
        from jarvis.db.jsonutils import dumpjson

        dumpjson(data=new_dat, filename=filename)
    return new_dat


def get_feats_hot_encoded(feature_names=keys, filename="feats_encoded.json"):
    """Get OneHotEncoded features."""
    # Deprecated
    # Kept for reference only
    from sklearn.preprocessing import OneHotEncoder
    import pandas as pd

    encoder = OneHotEncoder(categories="auto", sparse_output=False)
    # encoder = OneHotEncoder(categories="auto", sparse=False)
    dat = {}
    for i, j in chem_data.items():
        tmp = []
        for r, s in j.items():
            if r in feature_names:
                tmp.append(s)
        dat[Specie(i).Z] = tmp  # j.values()
    df = pd.DataFrame(dat)

    vals = []
    for i in range(len(df.values)):
        output = encoder.fit_transform(
            np.array(df.values[i], dtype="float").reshape(-1, 1)
        )  # .toarray()
        vals.extend(output.T)
    vals = np.array(vals, dtype="float").T
    cols = df.columns.tolist()
    new_dat = {}
    for i, j in zip(cols, vals):
        new_dat[i] = list(j)
    if filename is not None:
        from jarvis.db.jsonutils import dumpjson

        dumpjson(data=new_dat, filename=filename)
    return new_dat


x, y, z = get_specie_data()
info_z = {}
for i, j in y.items():
    info_z[j["Z"]] = i


def atomic_numbers_to_symbols(numbers=[1, 2, 3, 4]):
    """Convert atomic number array to atomic symbols."""
    symbs = []
    for i in numbers:
        symbs.append(info_z[i])
    return symbs


# get_digitized_feats_hot_encoded()
"""
if __name__ == "__main__":
    el = Specie("Al")
    #print(el.get_chgdescrp_arr)
    print(len(el.get_descrp_arr))
    #print (get_descrp_arr_name())
"""
