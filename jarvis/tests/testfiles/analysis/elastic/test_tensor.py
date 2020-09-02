from jarvis.analysis.elastic.tensor import ElasticTensor
import os, tarfile
from jarvis.io.lammps.outputs import parse_folder
from jarvis.io.vasp.outputs import Vasprun, Outcar

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


outcar = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88",
    "SiOptb88",
    "MAIN-ELASTIC-bulk@mp_149",
    "OUTCAR",
)


def test_vasp_et():
    out = Outcar(outcar)
    et = ElasticTensor(out.elastic_props()["cij"])
    print(et.to_dict())


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


fold_lammps = os.path.join(
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
    "bulk@mp-134",
)


def test_lammps_et():

    data = parse_folder(fold_lammps)
    print(data["elastic_tensor"]["raw_et_tensor"])


# test_lammps_et()
