from jarvis.core.atoms import Atoms, VacuumPadding
import os

poscar_path = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "POSCAR",
)


def test_basic_atoms():

    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    polar = Si.check_polar
    Si.props = ["a", "a"]
    vac_pad = VacuumPadding(Si)
    den_2d = round(vac_pad.get_effective_2d_slab().density, 2)
    den_0d = round(vac_pad.get_effective_molecule().density, 2)
    den_lll_red = round(Si.get_lll_reduced_structure().density, 2)
    strng = Si.get_string()
    scell_nat = Si.make_supercell([2, 2, 2]).num_atoms
    com = round(Si.get_center_of_mass()[0], 3)
    rem = (Si.make_supercell([2, 2, 2]).remove_site_by_index(site=0)).num_atoms

    d = Si.to_dict()
    angs_a = d["angles"][0]
    Si_2_den = Atoms(
        lattice_mat=d["lattice_mat"], coords=d["coords"], elements=d["elements"]
    ).density
    # print ('scell_nat', Si_2)
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


# test_basic_atoms()
# def test_basic_atoms():