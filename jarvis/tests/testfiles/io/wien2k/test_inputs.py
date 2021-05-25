from jarvis.io.wien2k.inputs import get_wien_kpoints
import os
from jarvis.core.atoms import Atoms


def test_kp():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    elements = ["Si", "Si"]
    s = Atoms(lattice_mat=box, coords=coords, elements=elements)
    data = get_wien_kpoints(atoms=s, write_file=True)
    cmd = "rm MyKpoints"
    os.system(cmd)
