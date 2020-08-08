from jarvis.io.calphad.write_decorated_poscar import get_selective_dyn_decorated_atoms
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar

def test_cal():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    elements = ["Si", "Si"]
    atoms = Atoms(lattice_mat=box, coords=coords, elements=elements)
    decorated_atoms, hall_number, wsymbols = get_selective_dyn_decorated_atoms(atoms)
    print (decorated_atoms)
    p = Poscar(decorated_atoms)
    p.write_file('POSCAR')
    #p2 = Poscar.from_file('POSCAR')
    #print (p2)
    assert (wsymbols, hall_number) == (["i", "i"], 63)


#test_cal()
