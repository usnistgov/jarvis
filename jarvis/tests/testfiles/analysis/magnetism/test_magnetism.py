from jarvis.core.atoms import Atoms
from jarvis.analysis.magnetism.magmom_setup import MagneticOrdering


def test_magnetism_setup():
    from jarvis.db.figshare import get_jid_data

    atoms = Atoms.from_dict(
        get_jid_data(jid="JVASP-78681", dataset="dft_3d")["atoms"]
    )
    mag = MagneticOrdering(atoms)
    symm_list, ss = mag.get_minimum_configs(min_configs=3)

    assert len(symm_list) == 3
    assert ss.num_atoms == 8
    mag_atoms = mag.get_mag_ions()
    assert mag_atoms == ["Mn"]
    tc = mag.tc_mean_field()
    assert round(tc["Tc"], 2) == round(3868.17, 2)
