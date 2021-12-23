from jarvis.core.atoms import (
    Atoms,
    VacuumPadding,
    get_supercell_dims,
    build_xanes_poscar,
)
import os
from jarvis.db.figshare import get_jid_data, data
import tarfile
import tempfile

new_file, filename = tempfile.mkstemp()


example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88.tgz",
)


example_fold = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()


poscar_path = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "SiOptb88",
    "POSCAR",
)

cif_example = os.path.join(
    os.path.dirname(__file__),
    "1000052.cif",
)
cif_example2 = os.path.join(
    os.path.dirname(__file__),
    "Bacomp.cif",
)
cif_example3 = os.path.join(
    os.path.dirname(__file__),
    "mock.cif",
)
cif_example4 = os.path.join(
    os.path.dirname(__file__),
    "exp_000034.cif",
)
cif_example5 = os.path.join(
    os.path.dirname(__file__),
    "1000000.cif",
)


def test_from_cif():
    a = Atoms.from_cif(cif_example)
    a = Atoms.from_cif(cif_example2)
    a = Atoms.from_cif(cif_example3)
    a = Atoms.from_cif(cif_example4)
    a = Atoms.from_cif(cif_example5)
    f = open(cif_example, "r")
    lines = f.read()
    f.close()
    x = Atoms.from_cif(from_string=lines)


def test_basic_atoms():

    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    dim = get_supercell_dims(Si)
    build_xanes_poscar(atoms=Si, filename_with_prefix=True)
    assert dim == [3, 3, 3]
    polar = Si.check_polar
    Si.props = ["a", "a"]
    vac_pad = VacuumPadding(Si)
    den_2d = round(vac_pad.get_effective_2d_slab().density, 2)
    den_0d = round(vac_pad.get_effective_molecule().density, 2)
    den_lll_red = round(Si.get_lll_reduced_structure().density, 2)
    strng = Si.get_string()
    scell_nat = Si.make_supercell([2, 2, 2]).num_atoms
    scell_nat2 = Si.make_supercell_matrix(
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
    ).num_atoms
    # print("scell_nat,scell_nat2", scell_nat, scell_nat2)
    # print(Si.make_supercell([2, 2, 2]))
    # print()
    # print(Si.make_supercell_matrix([[2, 0, 0], [0, 2, 0], [0, 0, 2]]))
    com = round(Si.get_center_of_mass()[0], 3)
    rem = (Si.make_supercell([2, 2, 2]).remove_site_by_index(site=0)).num_atoms
    prim = Si.get_primitive_atoms
    print(prim.cart_coords)
    assert round(prim.cart_coords[0][0], 2) == round(4.37815150, 2)
    # print ('raw_distance_matrix', prim.raw_distance_matrix)
    # print ('raw_distance_matrix', Si.raw_distance_matrix)
    # print ('distance_matrix', Si.pymatgen_converter().distance_matrix)
    assert round(prim.raw_distance_matrix[0][1], 2) == round(
        4.42386329832851, 2
    )
    asee = Si.ase_converter()
    print(prim.raw_angle_matrix)
    d = Si.to_dict()
    new_at = Atoms.from_dict(d)
    angs_a = d["angles"][0]
    Si_2_den = Atoms(
        lattice_mat=d["lattice_mat"],
        coords=d["coords"],
        elements=d["elements"],
    ).density
    Si_xyz = Si.get_xyz_string
    Si.write_xyz(filename="atoms.xyz")
    tmp = Atoms.from_xyz(filename="atoms.xyz")
    cmd = "rm atoms.xyz"
    os.system(cmd)
    Si.center_around_origin()
    # print ('scell_nat', Si_2)
    nb = Si.get_all_neighbors(5)
    assert (
        round(Si.volume, 2),
        Si.atomic_numbers,
        Si.num_atoms,
        Si.frac_coords[0][0],
        Si.cart_coords[0][0],
        round(Si.density, 2),
        Si.spacegroup(),
        Si.pymatgen_converter() != {},
        polar,
        Si.props[0],
        den_2d,
        den_0d,
        round(Si.packing_fraction, 2),
        Si.composition.to_dict(),
        strng != "",
        den_lll_red,
        scell_nat,
        com,
        rem,
        angs_a,
        round(Si_2_den, 2),
    ) == (
        40.03,
        [14, 14],
        2,
        0,
        0.0,
        2.33,
        "C2/m (12)",
        True,
        False,
        "a",
        0.35,
        0.01,
        0.28,
        {"Si": 2},
        True,
        2.33,
        16,
        0.679,
        15,
        60.0,
        2.33,
    )
    cc = Si.center()
    cc = Si.center(axis=[0, 0, 1])
    cc = Si.center(about=[0.5, 0.5, 0.5])

    m1 = Atoms.from_dict(get_jid_data("JVASP-6640")["atoms"])
    assert m1.check_polar == True
    print("Strain test")
    print(m1)
    m1.apply_strain(0.1)
    print(m1)
    assert m1.lattice_mat[2][2] == 32.8717576
    m1.apply_strain([0, 0, 0.1])
    assert m1.lattice_mat[2][2] == 36.158933360000006
    filename = "atoms.cif"
    m1.write_cif(filename)
    a = Atoms.from_cif(filename)
    filename = "POSCAR"
    m1.write_poscar(filename)
    m2 = Atoms.from_poscar(filename)

    filename = "atoms.xyz"
    m1.write_xyz(filename)
    m3 = Atoms.from_xyz(filename)

    cmd = "rm atoms.xyz POSCAR atoms.cif"
    os.system(cmd)


# test_basic_atoms()
# def test_basic_atoms():
