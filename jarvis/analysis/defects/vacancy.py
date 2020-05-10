"""
Modules for making point-defect vacancies
"""
import pprint
import random
from collections import OrderedDict
from jarvis.io.vasp.inputs import Poscar
from jarvis.core.lattice import Lattice
from jarvis.core.atoms import Atoms
import numpy as np
from jarvis.analysis.structure.spacegroup import Spacegroup3D


def rand_select(x=[]):
    uniq = list(set(x))
    info = {}
    for i, ii in enumerate(x):
        info.setdefault(ii, []).append(i)
    selected = {}
    for i, j in info.items():
        chosen = random.choice(j)
        selected.setdefault(i, chosen)
    # print (info)
    # print (selected)
    return selected


class Vacancy(object):
    def __init__(
        self,
        atoms=None,
        defect_structure=None,
        defect_index=None,
        wyckoff_multiplicity=None,
        symbol=None,
    ):
        """
        Get vacancy structures based on Wyckoff positions
        """
        self._atoms = atoms
        self._defect_index = defect_index
        self._defect_structure = defect_structure
        self._wyckoff_multiplicity = wyckoff_multiplicity
        self._symbol = symbol

    def generate_defects(
        self, enforce_c_size=10.0, on_conventional_cell=False, extend=1
    ):
        atoms = self._atoms
        if on_conventional_cell:
            atoms = Spacegroup3D(atoms).conventional_standard_structure
        a = atoms.lattice.lat_lengths()[0]
        b = atoms.lattice.lat_lengths()[1]
        c = atoms.lattice.lat_lengths()[2]
        if enforce_c_size is not None:
            dim1 = int(float(enforce_c_size) / float(a)) + extend
            dim2 = int(float(enforce_c_size) / float(b)) + extend
            dim3 = int(float(enforce_c_size) / float(c)) + extend
            # atoms = atoms.make_supercell([dim1, dim2, dim3])
            supercell_size = [dim1, dim2, dim3]

        element_list = list(set(atoms.elements))
        #print("s_size=", supercell_size, a, b, c, atoms.lattice.lat_lengths())
        #print("s_sizeomg=", atoms.pymatgen_converter().lattice.abc)
        #print(atoms)
        spg = Spacegroup3D(atoms)
        wyckoffs = spg._dataset["wyckoffs"]
        atoms.props = wyckoffs
        supercell = atoms.make_supercell(supercell_size)
        props = rand_select(supercell.props)
        vacs = []
        for i, j in props.items():
            defect_strt = supercell.remove_site_by_index(j)
            vac = Vacancy(
                atoms=supercell,
                defect_structure=defect_strt,
                defect_index=j,
                wyckoff_multiplicity=i,
                symbol=supercell.elements[j],
            )
            vacs.append(vac)
        return vacs

    def to_dict(self):
        d = OrderedDict()
        d["atoms"] = self._atoms
        d["defect_structure"] = self._defect_structure
        d["defect_index"] = self._defect_structure
        d["wyckoff_multiplicity"] = self._wyckoff_multiplicity
        d["symbol"] = self._symbol
        return d

    def __repr__(self, indent=2):
        return pprint.pformat(self.to_dict(), indent=indent)


"""
if __name__ == "__main__":

    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Si = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/2D-bulk/mp-1143_bulk_LDA/MAIN-ELASTIC-bulk@mp_1143/POSCAR"
    ).atoms
    vacs = Vacancy(atoms=Si).generate_defects()
    for i in vacs:
        print(i)
    spg = Spacegroup3D(Si)
    cvn = spg.conventional_standard_structure
    spg = Spacegroup3D(cvn)
    props = spg._dataset["wyckoffs"]
    Si.props = props
    ss = Si.make_supercell([2, 2, 2])
    props = ss.props
    # print (rand_select(props))
"""
