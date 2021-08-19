from jarvis.analysis.defects.substitutions import generate_defect
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os


def test_subs():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Al", "Al"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    v = generate_defect(atoms=Si)


# test_2d()
# test_vacancy()
