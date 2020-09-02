from jarvis.io.lammps.outputs import analyze_log
import os
from jarvis.io.lammps.outputs import parse_full_ff_folder
import tarfile

example_fold_tgz = os.path.join(
    os.path.dirname(__file__),
    "..",
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
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist",
)

if not os.path.isdir(example_fold):
    tar = tarfile.open(example_fold_tgz)
    tar.extractall(example_fold)
    tar.close()

example_fold_run = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist",
    "Al03.eam.alloy_nist",
)
def test_parse_full_ff_folder():
 abspath=os.path.abspath(example_fold_run)
 info = parse_full_ff_folder(abspath)
 print (info)


log_lammps = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold",
    "mp-134",
    "log.lammps",
)


def test_outputs():
    x = analyze_log(log_lammps)
    assert (x) == (
        -3.36,
        0.0,
        -3.36,
        120.1,
        120.7,
        121.2,
        58.7,
        57.9,
        57.9,
        29.5,
        29.5,
        30.0,
        1.6,
        9.1,
        -6.8,
        -0.0,
        -0.0,
        -0.0,
        -0.0,
        -0.0,
        -0.0,
        -0.0,
        -3.6,
    )

test_parse_full_ff_folder()
# test_outputs()
