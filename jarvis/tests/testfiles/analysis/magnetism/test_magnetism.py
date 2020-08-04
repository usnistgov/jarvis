from jarvis.core.atoms import Atoms
from jarvis.analysis.magnetism.magmom_setup import get_unique_magnetic_structures, tc_mean_field, get_mag_ions
def test_magnetism_setup():
    from jarvis.db.figshare import get_jid_data
    from jarvis.core.atoms import get_supercell_dims
    atoms = Atoms.from_dict(
        get_jid_data(jid="JVASP-78681", dataset="dft_3d")["atoms"]
    )
    dim = get_supercell_dims(atoms)
    print("dim=", dim)
    dim = [2, 2, 2]
    symm_list, ss = get_unique_magnetic_structures(
        atoms, supercell_dim=dim, magnetic_ions=["Mn"]
    )
    
    print("dim=", dim, len(symm_list))
    assert len(symm_list)==6
    assert ss.num_atoms == 16
    mag_atoms=get_mag_ions(atoms)
    assert mag_atoms==['Mn']
    tc = tc_mean_field(atoms=atoms)
    assert round(tc['Tc'],2)==round(7736.345414497057,2) 
