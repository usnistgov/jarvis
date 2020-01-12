"""
This module provides classes to specify atomic structure
"""
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie
from jarvis.core.lattice import Lattice
import importlib.util

amu_gm = 1.66054e-24
ang_cm = 1e-8


class Atoms(object):
    def __init__(
        self, lattice_mat=None, coords=None, elements=None, cartesian: bool = False
    ):
        """
        Create atomic structure with lattice, coordinates, atom type and other information
        >>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
        >>> coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
        >>> elements = ["Si", "Si"]
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
        >>> print(round(Si.volume,2))
        40.03
        >>> Si.composition
        {'Si': 2}
        >>> round(Si.density,2)
        2.33
        >>> Si.atomic_numbers
        [14, 14]
        >>> Si.num_atoms
        2
        >>> Si.frac_coords[0][0]
        0
        >>> Si.cart_coords[0][0]
        0.0
        >>> coords = [[0, 0, 0], [1.3575 , 1.22175, 1.22175]]
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements,cartesian=True)
        >>> round(Si.density,2)
        2.33
        >>> Si.spacegroup()
        'C2/m (12)'
        >>> Si.pymatgen_converter()!={}
        True
        """

        self.lattice_mat = lattice_mat
        self.coords = coords
        self.elements = elements
        if cartesian:
            self.cart_coords = self.coords
            self.frac_coords = Lattice(lattice_mat).frac_coords(self.cart_coords)
        else:
            self.frac_coords = self.coords
            self.cart_coords = Lattice(lattice_mat).cart_coords(self.frac_coords)

    @property
    def volume(self):
        m = self.lattice_mat
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol

    @property
    def composition(self):
        comp = {}
        for i in self.elements:
            comp[i] = comp.setdefault(i, 0) + 1
        return Composition(comp)

    @property
    def density(self):
        den = float(self.composition.weight * amu_gm) / (
            float(self.volume) * (ang_cm) ** 3
        )
        return den

    @property
    def atomic_numbers(self):
        numbers = []
        for i in self.elements:
            numbers.append(Specie(i).Z)
        return numbers

    @property
    def num_atoms(self):
        return len(self.coords)

    def pymatgen_converter(self):
        try:
            from pymatgen.core.structure import Structure

            return Structure(
                self.lattice_mat,
                self.elements,
                self.frac_coords,
                coords_are_cartesian=False,
            ) 
        except:pass

    def spacegroup(self, symprec=1e-3):
        #try:
            import spglib
            sg = spglib.get_spacegroup(
                (self.lattice_mat, self.frac_coords, self.atomic_numbers),
                symprec=symprec,
            )
            return sg
        #except:
        #    pass


# if __name__=='__main__':
#    box = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
#    coords = [[0, 0, 0], [0, 0, 0.5]]
#    elements = ["Si", "Si"]
#    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
#    comp=Si.composition
#    #print (comp,Si.density)
#    print (Si.atomic_numbers)
#   print (Si.pymatgen_converter().composition.weight,Si.composition.weight,Si.density)
