from jarvis.io.wannier.inputs import Wannier90win
from jarvis.core.atoms import Atoms
import os

win = os.path.join(os.path.dirname(__file__), "win.input")


def test_win():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Wannier90win(struct=Si, efermi=0.0).write_win(name=win)
    assert (os.path.isfile(win)) == (True)
    os.remove(win)


test_win()
