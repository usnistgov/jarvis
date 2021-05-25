"""
Coulomb matrix for Atoms.

Refer to: 10.1103/PhysRevLett.108.058301
"""

import numpy as np
from jarvis.core.specie import Specie


def coulomb_matrix(atoms="", max_dim=100):
    """Convert Atoms class to max_dim x max_dim matrix.

    Args:

        atoms: atoms object

        max_dim: maximum number of atoms=sqrt(max_dim)

    Returns:
          z: numpy array of 1 x max_dim dimension
    """
    natoms = atoms.num_atoms
    mat = np.zeros((natoms, natoms))

    for ii, i in enumerate(atoms.elements):
        for jj, j in enumerate(atoms.elements):
            if ii == jj:
                mat[ii, jj] = 0.5 * Specie(i).Z ** 2.4
            else:
                a = atoms.cart_coords[ii]
                b = atoms.cart_coords[jj]
                dist = np.linalg.norm(a - b)
                mat[ii, jj] = (Specie(i).Z * Specie(j).Z) / dist
    tmp = mat.ravel()
    if max_dim < len(tmp):
        print("WARNING: Increase max_dim")
    padding = max_dim - len(tmp)
    z = np.pad(tmp, (0, padding), "constant")
    return z


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    z = coulomb_matrix(Si)
    print(z, len(z))
"""
