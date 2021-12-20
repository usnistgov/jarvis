from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
import tempfile
from jarvis.io.qe.inputs import QEinfile


def test_inputs():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    print(Si)
    kp = Kpoints3D().automatic_length_mesh(
        lattice_mat=Si.lattice_mat, length=20
    )
    new_file, filename = tempfile.mkstemp()
    qe = QEinfile(Si, kp)
    qe.write_file(filename)
    kp = Kpoints3D().kpath(atoms=Si)
    qe = QEinfile(Si, kp)
    qe.write_file(filename)
    sp = qe.atomic_species_string()
    sp = qe.atomic_cell_params()
    assert qe.input_params["system"]["nat"] == 2
    
