from jarvis.db.figshare import get_jid_data
import numpy as np
from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.core.atoms import Atoms
import os, tarfile
from jarvis.io.lammps.outputs import parse_folder
from jarvis.io.vasp.outputs import Vasprun, Outcar
from jarvis.core.atoms import Atoms

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

poscar = os.path.join(
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
    "POSCAR",
)

atoms = Atoms.from_poscar(poscar)


def test_vasp_et():
    out = Outcar(outcar)
    et = ElasticTensor(out.elastic_props()["cij"])
    print(et.to_dict())
    theta = et.debye_temperature(atoms=atoms)
    print(
        et.is_brittle, et.is_ductile,et.cauchy_pressure, et.melting_temperature_metals
    )


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


def test_deb():
        x = get_jid_data(jid="JVASP-19821", dataset="dft_3d")
        et = ElasticTensor(et_tensor=np.array(x["elastic_tensor"]))
        atoms = Atoms.from_dict(x["atoms"])
        dd = et.debye_temperature(atoms=atoms)
        assert round(dd,2)==round(1047.547632064132,2)
        # test_lammps_et()
