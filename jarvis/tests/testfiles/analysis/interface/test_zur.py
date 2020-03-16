from jarvis.analysis.interface.zur import ZSLGenerator, mismatch_strts
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os

s1 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-652"))
s2 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-664"))


def test_zur():
    info = mismatch_strts(film=s1.atoms, subs=s2.atoms)
    # print (round(info['mismatch_u'],3), round(info['mismatch_angle'],3))
    assert (round(info["mismatch_u"], 3), round(info["mismatch_angle"], 3)) == (
        0.002,
        0.0,
    )


# test_zur()
