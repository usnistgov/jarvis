import os

import tarfile

example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist.tgz",
)


example_fold = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()



lammps_path = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold",
)

cwd = str(os.getcwd())


def test_lammps():
    os.chdir(lammps_path)
    cmd = str("python master.py")
    os.system(cmd)


os.chdir(cwd)
