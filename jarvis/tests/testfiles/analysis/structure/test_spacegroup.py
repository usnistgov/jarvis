from jarvis.analysis.structure.spacegroup import (
    Spacegroup3D,
    symmetrically_distinct_miller_indices,
    get_wyckoff_position_operators,
)
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os

s1 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "..","defects","POSCAR-667.vasp")).atoms
s2 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "..","..","io","wannier","POSCAR")).atoms
s3 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-tetragonal")).atoms
s4 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-Cmcm")).atoms

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
    spg = Spacegroup3D(atoms=s1)
    assert spg.space_group_number==191
    print (spg.space_group_number)
    cvn = spg.conventional_standard_structure
    spg = Spacegroup3D(atoms=s2)
    assert spg.space_group_number==166
    print (spg.space_group_number)
    cvn = spg.conventional_standard_structure
    spg = Spacegroup3D(atoms=s3)
    print (spg.space_group_number)
    cvn = spg.conventional_standard_structure
    assert spg.space_group_number==139
    spg = Spacegroup3D(atoms=s4)
    print (spg.space_group_number)
    cvn = spg.conventional_standard_structure
    assert spg.space_group_number==63

   

#test_spg()
