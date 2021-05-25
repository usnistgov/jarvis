from jarvis.io.lammps.inputs import LammpsData, LammpsInput
from jarvis.core.atoms import Atoms
from jarvis.analysis.structure.spacegroup import Spacegroup3D
import os
import tempfile
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


data = os.path.join(os.path.dirname(__file__), "lammps.data")
init = os.path.join(os.path.dirname(__file__), "init.mod")
inm = os.path.join(os.path.dirname(__file__), "in.main")
pot = os.path.join(
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
    "bulk@cellmax4",
    "potential.mod",
)


def test_inputs():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Al", "Al"]
    Al = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(atoms=Al)
    Al = spg.conventional_standard_structure
    box = LammpsData().atoms_to_lammps(atoms=Al)
    td = box.to_dict()
    fd = LammpsData.from_dict(td)
    box.write_file(filename=data)
    d = box.to_dict()
    # print ('d=',d)
    reload_d = LammpsData.from_dict(d)
    # print ('d=',reload_d)
    atoms = box.lammps_to_atoms()
    lmp = LammpsData().read_data(filename=data, potential_file=pot)
    lmp = LammpsData().atoms_to_lammps(atoms=atoms)
    pair_coeff = "abc"
    inp = LammpsInput(LammpsDataObj=lmp).write_lammps_in(
        lammps_in=init,
        lammps_in1=pot,
        lammps_in2=inm,
        parameters={
            "pair_style": "eam/alloy",
            "pair_coeff": pair_coeff,
            "atom_style": "charge",
            "control_file": "/users/knc6/inelast.mod",
        },
    )
    d=LammpsInput(LammpsDataObj=lmp).to_dict()
    d=LammpsInput.from_dict(d)
    assert (lmp._lammps_box[1], atoms.num_atoms) == (5.43, 8)
    # os.remove(data)
    # os.remove(init)
    # os.remove(inm)
    # os.remove(pot)


# test_inputs()
