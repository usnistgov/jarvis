"""Run multiple jobs."""
from jarvis.tasks.qe.super import SuperCond
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
from jarvis.db.figshare import data, get_jid_data
from jarvis.io.vasp.inputs import Poscar, Incar, Potcar
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
from jarvis.db.jsonutils import loadjson, dumpjson
import os

jids = ["JVASP-816"]
run_dir = "/home/knc6/Software/qe/jarvis/jarvis/tasks/qe"


def write_qejob(pyname="job.py", job_json=""):
    """Write template job.py with VaspJob.to_dict() job.json."""
    f = open(pyname, "w")
    f.write("from jarvis.tasks.qe.super import SuperCond\n")
    f.write("from jarvis.db.jsonutils import loadjson\n")
    f.write('d=loadjson("' + str(job_json) + '")\n')
    f.write("v=SuperCond.from_dict(d)\n")
    f.write("v.runjob()\n")
    f.close()


for i in jids:
    dir_name = os.path.join(run_dir, i + "_SUPER")
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    os.chdir(dir_name)
    dat = get_jid_data(jid=i, dataset="dft_3d")
    atoms = Atoms.from_dict(dat["atoms"])
    kp = Kpoints3D().automatic_length_mesh(
        lattice_mat=atoms.lattice_mat, length=10
    )
    sup = SuperCond(atoms=atoms, kp=kp).to_dict()
    dumpjson(data=sup, filename="sup.json")
    write_qejob(job_json=os.path.abspath("sup.json"))
    path = (
        "module load intel/2015 openmpi/1.10.2/intel-15"
        + " \nsource activate qe \npython "
        + os.getcwd()
        + "/job.py"
    )
    Queue.pbs(
        job_line=path,
        jobname=i,
        directory=os.getcwd(),
        submit_cmd=["qsub", "submit_job"],
    )

    os.chdir("..")
