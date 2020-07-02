"""Modules for handling PDB Protein Data Bank files."""
import numpy as np
from jarvis.core.atoms import Atoms, VacuumPadding


def read_pdb(filename=""):
    """Read PDB file and make Atoms object."""
    f = open(filename, "r")
    lines = f.read().splitlines()
    f.close()
    coords = []
    species = []
    for i in lines:
        tmp = i.split()
        if "ATOM " in i and "REMARK" not in i and "SITE" not in i:
            coord = [float(tmp[6]), float(tmp[7]), float(tmp[8])]
            coords.append(coord)
            species.append(tmp[11])
            # print (coord,tmp[11])
    coords = np.array(coords)

    max_c = np.max(coords, axis=0)
    min_c = np.min(coords, axis=0)
    box = np.zeros((3, 3))
    for j in range(3):
        box[j, j] = abs(max_c[j] - min_c[j])
    pdb = Atoms(lattice_mat=box, elements=species,
                coords=coords, cartesian=True)
    mol = VacuumPadding(pdb, vacuum=20.0).get_effective_molecule()
    return mol


"""
if __name__ == "__main__":
    pdb = read_pdb('/cluster/users/knc6/pdb/pdb101d.ent')
    print (pdb)
    import sys
    sys.exit()

"""
