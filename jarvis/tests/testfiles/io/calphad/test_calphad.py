from jarvis.io.calphad.write_decorated_poscar import get_selective_dyn_decorated_atoms
from jarvis.core.atoms import Atoms


def test_cal():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    elements = ["Si", "Si"]
    atoms = Atoms(lattice_mat=box, coords=coords, elements=elements)
    decorated_atoms, hall_number, wsymbols = get_selective_dyn_decorated_atoms(atoms)
    assert (wsymbols, hall_number) == (["i", "i"], 63)


# test_cal()
