from jarvis.io.lammps.inputs import LammpsData, LammpsInput
from jarvis.core.atoms import Atoms
from jarvis.analysis.structure.spacegroup import Spacegroup3D
import os

data = os.path.join(os.path.dirname(__file__), "lammps.data")
init = os.path.join(os.path.dirname(__file__), "init.mod")
pot = os.path.join(os.path.dirname(__file__), "potential.mod")
inm = os.path.join(os.path.dirname(__file__), "in.main")


def test_inputs():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(atoms=Si)
    Si = spg.conventional_standard_structure
    box = LammpsData().atoms_to_lammps(atoms=Si)
    box.write_file(data)
    atoms = box.lammps_to_atoms()
    lmp = LammpsData().read_data(data)
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

    assert (lmp._lammps_box[1], atoms.num_atoms) == (5.43, 8)
    os.remove(data)
    os.remove(init)
    os.remove(inm)
    os.remove(pot)


test_inputs()
