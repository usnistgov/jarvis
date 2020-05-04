import os
from jarvis.analysis.thermodynamics.energetics import form_enp

from jarvis.io.vasp.inputs import Poscar

p = os.path.join(os.path.dirname(__file__), "POSCAR")


def test_form_en():
    atoms = Poscar.from_file(p).atoms

    total_energy = -9.974648
    fen = form_enp(atoms=atoms, total_energy=total_energy)
    assert fen == -0.40429
    # print ('fen', fen)


# test_form_en()
