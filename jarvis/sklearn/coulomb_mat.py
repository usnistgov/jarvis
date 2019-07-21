import numpy as np
from pymatgen.core.structure import Structure

"""
Coulomb matrix for a structure
Refer: 10.1103/PhysRevLett.108.058301
"""


def coulomb_matrix(strt="", max_dim=100):
    natoms = len(strt)
    mat = np.zeros((natoms, natoms))
    for ii, i in enumerate(strt):
        for jj, j in enumerate(strt):
            if ii == jj:
                mat[ii, jj] = 0.5 * i.specie.Z ** 2.4
            else:
                dist = strt.get_distance(ii, jj)
                mat[ii, jj] = (i.specie.Z * j.specie.Z) / dist
    tmp = mat.ravel()
    z = np.pad(tmp, (0, max_dim), "constant")
    return z


if __name__ == "__main__":
    s = Structure.from_file("POSCAR")
    coulomb_matrix(strt=s)
