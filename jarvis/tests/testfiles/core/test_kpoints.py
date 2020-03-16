from jarvis.core.kpoints import Kpoints3D
from jarvis.core.atoms import Atoms


def test_kp():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    lattice_mat = Si.lattice_mat
    kp = Kpoints3D().automatic_length_mesh(lattice_mat=lattice_mat, length=20)
    sym = kp.high_symm_path(Si)._path
    x, y = kp.interpolated_points(Si)
    # for i,j in zip(x,y):
    #   print (i,j)
    # print (len(x))
    assert (len(x), x[0][0], y[0]) == (166, 0, "\Gamma")


# test_kp()
