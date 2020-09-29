from jarvis.tasks.vasp.vasp import VaspJob, write_jobfact_optb88vdw
from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
from jarvis.db.jsonutils import dumpjson
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
import os
from jarvis.tasks.vasp.vasp import JobFactory, write_jobfact_optb88vdw
from jarvis.db.figshare import get_jid_data


tmp_dict = get_jid_data(jid="JVASP-816", dataset="dft_3d")["atoms"]
atoms = Atoms.from_dict(tmp_dict)


vasp_cmd = "mpirun /users/knc6/VASP/vasp54/src/vasp.5.4.1DobbySOC2/bin/vasp_std"
copy_files = ["/users/knc6/bin/vdw_kernel.bindat"]
submit_cmd = ["qsub",  "submit_job"]
jids = ["JVASP-1002", "JVASP-1067"]

for jid in jids:
    d = get_jid_data(jid=jid, dataset="dft_3d")
    atoms = Atoms.from_dict(d["atoms"]).get_primitive_atoms
    mat = Poscar(atoms)
    mat.comment = "bulk@" + str(jid)
    cwd_home = os.getcwd()
    dir_name = d["jid"] + "_" + str("PBEBO")
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    os.chdir(dir_name)
    job = JobFactory(vasp_cmd=vasp_cmd, poscar=mat, copy_files=copy_files,)
    dumpjson(data=job.to_dict(), filename="job_fact.json")
    write_jobfact_optb88vdw(pyname="job_fact.py", job_json="job_fact.json")

    # Example job commands, need to change based on your cluster
    job_line = (
        "source ~/anaconda2/envs/my_jarvis/bin/activate my_jarvis \n"
        + "python job_fact.py"
    )
    name = jid
    directory = os.getcwd()
    Queue.pbs(
        job_line=job_line,
        jobname=name,
        directory=directory,
        submit_cmd=submit_cmd,
    )
    os.chdir(cwd_home)

