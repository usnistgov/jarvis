from jarvis.io.boltztrap.outputs import BoltzTrapOutput
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


bpath = os.path.join(
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
    "boltztrap",
)


def test_out():
    b = BoltzTrapOutput(path=bpath)
    x = b.to_dict()["condtens_fixdoping"]
    td = b.to_dict()
    fd = BoltzTrapOutput.from_dict(td)
    tmp = x["p"][300.0][4.1e-07]["cond"][0]
    assert tmp == 609055450000000.0


# test_out()
