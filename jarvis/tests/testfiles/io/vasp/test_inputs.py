from jarvis.io.vasp.inputs import Poscar, Incar

import os

pos = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "MAIN-RELAX-bulk@mp_149",
    "CONTCAR",
)
inc = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "MAIN-RELAX-bulk@mp_149",
    "INCAR",
)


def test_inputs():
    p = Poscar.from_file(pos)
    i = Incar.from_file(inc)
    assert (round(p.atoms.density, 2), i.to_dict()["ISIF"]) == (2.25, "3")
