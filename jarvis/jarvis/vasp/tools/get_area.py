import numpy as np
import pymatgen as mg
structure = mg.Structure.from_file("POSCAR")
m=structure.lattice.matrix
area = np.linalg.norm(np.cross(m[0], m[1]))
print area
