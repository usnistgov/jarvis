from jarvis.core.kpoints import (
    Kpoints3D,
    generate_kpath,
    generate_kgrid,
    HighSymmetryKpoint3DFactory,
)
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os
import tempfile


s1 = Poscar.from_file(
    os.path.join(os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR")
).atoms
s2 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-Aem2"
    )
).atoms
s3 = Poscar.from_file(
    os.path.join(os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-C2m")
).atoms
s4 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-Cmcm"
    )
).atoms
s5 = Poscar.from_file(
    os.path.join(os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-P-1")
).atoms
s6 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-tetragonal"
    )
).atoms
s7 = Poscar.from_file(
    os.path.join(os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-Pc")
).atoms


def test_kp():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    lattice_mat = Si.lattice_mat
    kp = Kpoints3D().automatic_length_mesh(lattice_mat=lattice_mat, length=20)
    sym = kp.high_symm_path(Si)._path
    x, y = kp.interpolated_points(Si)
    kpath = generate_kpath(kpath=[[0, 0, 0], [0, 0.5, 0.5]], num_k=5)
    kps = generate_kgrid(grid=[5, 5, 5])
    new_file, filename = tempfile.mkstemp()
    kp.write_file(filename)
    sym = kp.high_symm_path(s1)._path
    sym = kp.high_symm_path(s2)._path
    sym = kp.high_symm_path(s3)._path
    sym = kp.high_symm_path(s4)._path
    sym = kp.high_symm_path(s5)._path
    sym = kp.high_symm_path(s6)._path
    sym = kp.high_symm_path(s7)._path
    assert (len(x), x[0][0], y[0], kpath[0][0], kps[0][0]) == (
        166,
        0,
        "\Gamma",
        0.0,
        0.0,
    )


def test_highsym():
    kpts = HighSymmetryKpoint3DFactory().cubic()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().fcc()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().bcc()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().tet()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().bctet1(3, 2)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().bctet2(3, 2)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orc()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcf1(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcf2(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcf3(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orci(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcc(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().hex()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().rhl1(47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().rhl2(47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mcl(2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc1(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc2(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc3(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc4(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc5(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().tria()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().trib()._kpoints["\\Gamma"]
    assert kpts[0] == 0


# test_highsym()
# test_kp()
