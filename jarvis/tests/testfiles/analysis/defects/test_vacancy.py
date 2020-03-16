from jarvis.analysis.defects.vacancy import Vacancy
from jarvis.core.atoms import Atoms


def test_vacancy():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    vacs = Vacancy(atoms=Si).generate_defects(enforce_c_size=5.0)
    # print ( len(vacs),(vacs[0].to_dict()["defect_structure"].num_atoms))
    assert (len(vacs), vacs[0].to_dict()["defect_structure"].num_atoms) == (1, 127)


# test_vacancy()
