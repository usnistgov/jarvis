import os

import tarfile

example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
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
    "examples",
    "vasp",
    "SiOptb88",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()




vasp_path = os.path.join(
    os.path.dirname(__file__), "..", "..", "..", "examples", "vasp", "SiOptb88", "SiOptb88"
)

cwd = str(os.getcwd())


def test_vasp():
    os.chdir(vasp_path)
    cmd = str("python master.py")
    os.system(cmd)


os.chdir(cwd)
