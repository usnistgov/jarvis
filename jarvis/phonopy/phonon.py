from __future__ import unicode_literals, print_function

"""
Helper function for running LAMMPS
Used for defects, surface and phonon calculations
"""

from ase import *
import numpy as np
import os

try:
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_CONSTANTS, write_FORCE_CONSTANTS
    from phonopy.structure.atoms import Atoms as PhonopyAtoms
except:
    pass
import glob
from pymatgen.core.structure import Structure


def get_phonopy_atoms(mat=None):
    """
    Helper function to convert pymatgen structure object to phonopy atoms

    Args:
        mat: pymatgen structure object
    Returns:
           phonopy atoms object
    """
    symbols = [str(site.specie.symbol) for site in mat]
    positions = [site.coords for site in mat]
    cell = mat.lattice.matrix
    p = PhonopyAtoms(symbols=symbols, positions=positions, pbc=True, cell=cell)
    return p


if __name__ == "__main__":
    pos = str(
        os.path.join(
            os.path.dirname(__file__),
            "../vasp/examples/SiOptb88/MAIN-MBJ-bulk@mp_149/POSCAR",
        )
    )
    mat = Structure.from_file(pos)
    p = get_phonopy_atoms(mat)
    print(p)
