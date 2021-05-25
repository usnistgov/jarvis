from jarvis.io.boltztrap.inputs import WriteInputs
from jarvis.io.vasp.outputs import Vasprun
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




vrun = os.path.join(
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
    "vasprun.xml",
)
energy = os.path.join(os.path.dirname(__file__), "boltztrap.energyso")
struct = os.path.join(os.path.dirname(__file__), "boltztrap.struct")
intrans = os.path.join(os.path.dirname(__file__), "boltztrap.intrans")


def test_inp():
    inp = WriteInputs(vasprun_path=vrun)
    inp.write_energy(filename=energy)
    inp.write_struct(filename=struct)
    inp.write_intrans(filename=intrans)
    assert (
        os.path.isfile(energy),
        os.path.isfile(struct),
        os.path.isfile(intrans),
    ) == (True, True, True)
    os.remove(energy)
    os.remove(struct)
    os.remove(intrans)


#test_inp()
