from jarvis.io.vasp.inputs import Poscar, Kpoints, Incar

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

kp1 = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "MAIN-RELAX-bulk@mp_149",
    "KPOINTS",
)

kp2 = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "MAIN-BAND-bulk@mp_149",
    "KPOINTS",
)

def test_inputs():
    p = Poscar.from_file(pos)
    i = Incar.from_file(inc)
    assert (round(p.atoms.density, 2), i.to_dict()["ISIF"]) == (2.25, "3")

def test_kpoints():
    kp_mesh = Kpoints(filename=kp1).kpoints    
    kp_bz = Kpoints(filename=kp2).kpoints    
