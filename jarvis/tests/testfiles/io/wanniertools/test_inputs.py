import os
from jarvis.io.wanniertools.inputs import WTin
from jarvis.io.vasp.inputs import Poscar

pos = os.path.join(os.path.dirname(__file__), "..", "wannier", "POSCAR")
wout = os.path.join(os.path.dirname(__file__), "..", "wannier", "wannier90.wout")
wtin_file = os.path.join(os.path.dirname(__file__), "wt.in")


def test_input():
    p = Poscar.from_file(pos).atoms
    wtin = WTin(wtin=wtin_file, atoms=p, wannierout=wout).write_wt_in()
    assert (os.path.isfile(wtin_file)) == (True)
    os.remove(wtin_file)


# test_input()
