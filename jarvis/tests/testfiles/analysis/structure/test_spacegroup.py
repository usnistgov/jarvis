from jarvis.analysis.structure.spacegroup import (
    Spacegroup3D,
    symmetrically_distinct_miller_indices,
    get_wyckoff_position_operators,
)
from jarvis.core.atoms import Atoms


def test_spg():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(atoms=Si)  # .spacegroup_data()
    cvn = spg.conventional_standard_structure
    ml = symmetrically_distinct_miller_indices(max_index=3, cvn_atoms=cvn)
    # print ('ml',ml)
    x = get_wyckoff_position_operators(488)
    # print  (x['wyckoff'][0]['letter'])
    assert (
        spg.space_group_number,
        spg.space_group_symbol,
        cvn.num_atoms,
        ml[0][0],
        x["wyckoff"][0]["letter"],
    ) == (227, "Fd-3m", 8, 1, "l")


# test_spg()
