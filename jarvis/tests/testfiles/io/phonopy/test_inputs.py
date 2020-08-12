from jarvis.core.atoms import Atoms
from jarvis.io.phonopy.inputs import PhonopyInputs
from jarvis.analysis.structure.spacegroup import Spacegroup3D
import os

def test_input():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(Si)
    Si_cvn = spg.conventional_standard_structure
    p=PhonopyInputs(atoms=Si_cvn).generate_all_files()
    cmd = 'rm *.conf'
    os.system(cmd)

