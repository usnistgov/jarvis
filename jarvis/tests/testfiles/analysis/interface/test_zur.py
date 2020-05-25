from jarvis.analysis.interface.zur import ZSLGenerator,get_hetero_type, mismatch_strts, get_hetero,make_interface
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os

s1 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-652"))
s2 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-664"))
s3 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-688"))
s4 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-75175"))


def test_zur():
    info = mismatch_strts(film=s4.atoms, subs=s2.atoms)
    a = info['film_sl'].center_around_origin()
    b = info['subs_sl'].center_around_origin()
    combined  = get_hetero(a,b,seperation=3).center_around_origin()
    info = make_interface(film=s4.atoms,subs=s2.atoms)
    combined = info['interface']
    print (combined)
    info = mismatch_strts(film=s2.atoms, subs=s3.atoms)
    a = info['film_sl'].center_around_origin()
    b = info['subs_sl'].center_around_origin()
    combined  = get_hetero(a,b,seperation=3).center_around_origin()
    info = make_interface(film=s2.atoms,subs=s3.atoms)
    combined = info['interface']
    print (combined)
    info = mismatch_strts(film=s1.atoms, subs=s2.atoms)
    a = info['film_sl'].center_around_origin()
    b = info['subs_sl'].center_around_origin()
    combined  = get_hetero(a,b,seperation=3).center_around_origin()
    info = make_interface(film=s1.atoms,subs=s2.atoms)
    combined = info['interface']
    print (combined)
    assert (round(info["mismatch_u"], 3), round(info["mismatch_angle"], 3)) == (
        0.002,
        0.0,
    )

def test_type():
    B = {}
    B["scf_vbm"] = -5
    B["scf_cbm"] = -4
    B["avg_max"] = -2
    A = {}
    A["scf_vbm"] = -7
    A["scf_cbm"] = -6
    A["avg_max"] = -2
    int_type,stack = get_hetero_type(A=A,B=B)       
    assert (int_type,stack)==('III','BA')
