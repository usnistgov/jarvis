from jarvis.io.wannier.inputs import Wannier90win
from jarvis.core.atoms import Atoms
import os
from jarvis.io.vasp.inputs import Poscar

win = os.path.join(os.path.dirname(__file__), "win.input")


s1 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "analysis",
        "structure",
        "POSCAR",
    )
).atoms
s2 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "analysis",
        "structure",
        "POSCAR-Cmcm",
    )
).atoms
s3 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "analysis",
        "structure",
        "POSCAR-tetragonal",
    )
).atoms


def test_win():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Wannier90win(struct=Si, efermi=0.0).write_win(name=win)
    w = Wannier90win(struct=Si, efermi=0.0)
    td = w.to_dict()
    fd = Wannier90win.from_dict(td)
    assert (os.path.isfile(win)) == (True)
    os.remove(win)

    Wannier90win(struct=s1, efermi=0.0).write_win(name=win)
    assert (os.path.isfile(win)) == (True)
    Wannier90win(struct=s1, efermi=0.0).write_hr_win(prev_win=win)
    os.remove(win)

    Wannier90win(struct=s2, efermi=0.0).write_win(name=win)
    assert (os.path.isfile(win)) == (True)
    os.remove(win)

    Wannier90win(struct=s3, efermi=0.0).write_win(name=win)
    assert (os.path.isfile(win)) == (True)
    os.remove(win)
    cmd = 'rm *.win'
    os.system(cmd)

