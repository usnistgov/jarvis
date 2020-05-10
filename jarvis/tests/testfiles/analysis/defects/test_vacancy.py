from jarvis.analysis.defects.vacancy import Vacancy
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os

C = os.path.join(os.path.dirname(__file__), "POSCAR-667.vasp")


def test_vacancy():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    vacs = Vacancy(atoms=Si).generate_defects(enforce_c_size=10.0)
    # print ( len(vacs),(vacs[0].to_dict()["defect_structure"].num_atoms))
    # print (vacs[0].to_dict()["defect_structure"])
    assert (len(vacs), vacs[0].to_dict()["defect_structure"].num_atoms) == (1, 53)


def test_2d():
    p = Poscar.from_file(C).atoms
    vacs = Vacancy(atoms=p).generate_defects(enforce_c_size=10.0, extend=1)
    # print (vacs[0]._defect_structure)
    assert (vacs[0].to_dict()["defect_structure"]).num_atoms == 49


# test_2d()
# test_vacancy()
