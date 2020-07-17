"""Modules for making point-defect vacancies."""
import pprint
from collections import OrderedDict
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.utils import rand_select
from jarvis.core.atoms import Atoms


class Vacancy(object):
    """Obtain vacancy defects in Atoms class using Wyckoff data."""

    def __init__(
        self,
        atoms=None,
        defect_structure=None,
        defect_index=None,
        wyckoff_multiplicity=None,
        symbol=None,
    ):
        """
        Initialize the method.

        Arguments are given below.
        Args:
            atoms: jarvis.core.Atoms object.

            defect_index: atoms index for defect.

            defect_structure:  Atoms with defect.

            wyckoff_multiplicity: Wyckoff multiplicity.

            symbol: Elemenyt symbol.
        """
        self._atoms = atoms
        self._defect_index = defect_index
        self._defect_structure = defect_structure
        self._wyckoff_multiplicity = wyckoff_multiplicity
        self._symbol = symbol

    @classmethod
    def from_dict(self, d={}):
        """Load from a dictionary."""
        return Vacancy(
            atoms=Atoms.from_dict(d["atoms"]),
            defect_structure=Atoms.from_dict(d["defect_structure"]),
            defect_index=d["defect_index"],
            wyckoff_multiplicity=d["wyckoff_multiplicity"],
            symbol=d["symbol"],
        )

    def generate_defects(
        self, enforce_c_size=10.0, on_conventional_cell=False, extend=1
    ):
        """Provide function to generate defects."""
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
        """Convert to a dictionary."""
        d = OrderedDict()
        d["atoms"] = self._atoms.to_dict()
        if self._defect_structure is not None:
            d["defect_structure"] = self._defect_structure.to_dict()
        else:
            d["defect_structure"] = None
        d["defect_index"] = self._defect_index
        d["wyckoff_multiplicity"] = self._wyckoff_multiplicity
        d["symbol"] = self._symbol
        return d

    def __repr__(self, indent=2):
        """Representation of the class as dict."""
        return pprint.pformat(self.to_dict(), indent=indent)


"""
if __name__ == "__main__":
    from jarvis.io.vasp.inputs import Poscar

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
