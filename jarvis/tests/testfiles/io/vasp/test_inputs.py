from jarvis.io.vasp.inputs import (
    Poscar,
    Kpoints,
    Incar,
    Potcar,
    IndividualPotcarData,
    find_ldau_magmom,
)

import tempfile
import os
import tarfile
from jarvis.db.figshare import data
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import get_nelect

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
    pstr = p.to_string()

    td = p.to_dict()
    print("td is:", td)
    fd = Poscar.from_dict(td)
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
    # TODO iu =add_ldau_incar(use_incar_dict=d,atoms=p.atoms)
    ii = Incar.from_dict(d)
    ii.write_file("INCAR")
    print(ii)
    pot = os.path.join(
        os.path.dirname(__file__),
        "POT_GGA_PAW_PBE",
        "Xe",
        "POTCAR",
    )
    potc = IndividualPotcarData.from_file(pot)
    print(potc)
    os.environ["VASP_PSP_DIR"] = os.path.join(os.path.dirname(__file__))
    new_file, filename = tempfile.mkstemp()
    pot = Potcar(elements=["Xe"])
    td = pot.to_dict()
    fd = Potcar.from_dict(td)
    print(pot)
    pot.write_file(filename)
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Xe", "Xe"]
    xe = Atoms(lattice_mat=box, coords=coords, elements=elements)
    p = Potcar.from_atoms(atoms=xe)


def test_kpoints():
    new_file, filename = tempfile.mkstemp()
    kp_mesh = Kpoints(filename=kp1).kpoints
    kp_bz = Kpoints(filename=kp2).kpoints
    Kpoints(filename=kp2).kpoints.write_file(filename)


def test_ldau():
    d = data("dft_2d")
    #d = data("dft_3d")
    for i in d:
        if i["jid"] == "JVASP-76195":
            atoms = Atoms.from_dict(i["atoms"])
            ld = find_ldau_magmom(atoms=atoms, lsorbit=True)
            ld = find_ldau_magmom(atoms=atoms, lsorbit=False)
            nelect = get_nelect(atoms)
        if i["jid"] == "JVASP-76195":
            atoms = Atoms.from_dict(i["atoms"])
            ld = find_ldau_magmom(atoms=atoms, lsorbit=True)
            #assert ld["LDAUU"] == "3.0 0"
            ld = find_ldau_magmom(atoms=atoms, lsorbit=False)
