from jarvis.tasks.lammps.lammps import LammpsJob, JobFactory
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
from jarvis.analysis.structure.spacegroup import Spacegroup3D


# atoms = Atoms.from_poscar('POSCAR')
# Get Aluminum FCC
tmp_dict = get_jid_data(jid="JVASP-816", dataset="dft_3d")["atoms"]
atoms = Atoms.from_dict(tmp_dict)
# Get conventional cell
spg = Spacegroup3D(atoms)
cvn_atoms = spg.conventional_standard_structure
ff = "/users/knc6/Software/LAMMPS/lammps-master/potentials/Al_zhou.eam.alloy"
mod = "/users/knc6/Software/Devs/jarvis/jarvis/tasks/lammps/templates/inelast.mod"
cmd = "/users/knc6/Software/LAMMPS/lammps-master/src/lmp_serial<in.main>out"
parameters = {
    "pair_style": "eam/alloy",
    "pair_coeff": ff,
    "atom_style": "charge",
    "control_file": mod,
}
# Test LammpsJob
lmp = LammpsJob(
    atoms=cvn_atoms, parameters=parameters, lammps_cmd=cmd, jobname="Test"
).runjob()

# Test high-throughput
job_fact = JobFactory(pair_style="eam/alloy", name="my_first_lammps_run")
job_fact.all_props_eam_alloy(atoms=cvn_atoms, ff_path=ff, lammps_cmd=cmd)
