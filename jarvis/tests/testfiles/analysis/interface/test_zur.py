from jarvis.analysis.interface.zur import ZSLGenerator, mismatch_strts, get_hetero
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os

s1 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-652"))
s2 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-664"))
s3 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-688"))


def test_zur():
    info = mismatch_strts(film=s2.atoms, subs=s3.atoms)
    a = info['film_sl'].center_around_origin()
    b = info['subs_sl'].center_around_origin()
    combined  = get_hetero(a,b,seperation=3).center_around_origin()
    print (combined)
    info = mismatch_strts(film=s1.atoms, subs=s2.atoms)
    a = info['film_sl'].center_around_origin()
    b = info['subs_sl'].center_around_origin()
    combined  = get_hetero(a,b,seperation=3).center_around_origin()
    print (combined)
    assert (round(info["mismatch_u"], 3), round(info["mismatch_angle"], 3)) == (
        0.002,
        0.0,
    )


test_zur()
