from jarvis.io.vasp.inputs import (
    Poscar,
    Kpoints,
    Incar,
    Potcar,
    IndividualPotcarData,
)

import tempfile
import os
import tarfile

example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88.tgz",
)


example_fold = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()


pos = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
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
    "SiOptb88",
    "MAIN-BAND-bulk@mp_149",
    "KPOINTS",
)


def test_inputs():
    p = Poscar.from_file(pos)
    print(p)
    new_file, filename = tempfile.mkstemp()
    p.write_file(filename)
    i = Incar.from_file(inc)
    i.update({"ENCUT": 1000})
    assert (round(p.atoms.density, 2), i.to_dict()["ISIF"]) == (2.25, "3")
    f = open(pos, "r")
    lines = f.read()
    f.close()
    p = Poscar.from_string(lines)
    assert (round(p.atoms.density, 2), i.to_dict()["ISIF"]) == (2.25, "3")
    d = i.to_dict()
    ii = Incar.from_dict(d)
    pot = os.path.join(
        os.path.dirname(__file__), "POT_GGA_PAW_PBE", "Xe", "POTCAR",
    )
    potc = IndividualPotcarData.from_file(pot)
    print(potc)
    os.environ["JARVIS_VASP_PSP_DIR"] = os.path.join(os.path.dirname(__file__))
    new_file, filename = tempfile.mkstemp()
    pot = Potcar(elements=["Xe"])
    print(pot)
    pot.write_file(filename)


def test_kpoints():
    new_file, filename = tempfile.mkstemp()
    kp_mesh = Kpoints(filename=kp1).kpoints
    kp_bz = Kpoints(filename=kp2).kpoints
    Kpoints(filename=kp2).kpoints.write_file(filename)
