"""Run multiple jobs."""
from jarvis.tasks.qe.super import SuperCond
from jarvis.core.utils import get_factors
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import data
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
from jarvis.db.jsonutils import dumpjson
import os
from jarvis.analysis.structure.spacegroup import Spacegroup3D

# import glob
# from jarvis.db.jsonutils import loadjson

my_data = data("dft_3d")


def get_jid_data(jid="JVASP-667", dataset="dft_2d"):
    """Get info for a jid and dataset."""
    # d = data(dataset)
    d = my_data
    for i in d:
        if i["jid"] == jid:
            return i


qe_cmd = "mpirun -np 16 /home/kfg/codes/q-e-qe-6.5/bin/pw.x"
qe_cmd = "/home/kfg/codes/q-e-qe-6.5/bin/pw.x"

run_dir = "/wrk/knc6/Super"
run_dir = "/wrk/knc6/CDVAE_SUP"


def non_prime_kpoints(kpts=[]):
    """Get non prime kpoints."""
    mem = []
    for i in kpts:
        facts = get_factors(i)
        if len(facts) == 1:
            val = i + 1
        else:
            val = i
        mem.append(val)
    return mem


def write_qejob(pyname="job.py", job_json=""):
    """Write template job.py with VaspJob.to_dict() job.json."""
    f = open(pyname, "w")
    f.write("from jarvis.tasks.qe.super import SuperCond\n")
    f.write("from jarvis.db.jsonutils import loadjson\n")
    f.write('d=loadjson("' + str(job_json) + '")\n')
    f.write("v=SuperCond.from_dict(d)\n")
    f.write("v.runjob()\n")
    f.close()


jids = ["JVASP-816", "JVASP-19821"]
jids = [
    "POSCAR-AlN2Sc.vasp",
    "POSCAR-CrV6Pt.vasp",
    "POSCAR-NbNiV6.vasp",
    "POSCAR-RuNb.vasp",
]
submit_job = True
use_preconverged_kpoints = False
for i in jids:
    try:
        print("jid", i)
        if "POSCAR" not in i:
            dat = get_jid_data(jid=i, dataset="dft_3d")
            a_atoms = Atoms.from_dict(dat["atoms"])
        else:
            a_atoms = Atoms.from_poscar(i)
        dir_name = os.path.join(run_dir, i + "_SUPER")
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        os.chdir(dir_name)

        atoms = Spacegroup3D(a_atoms).refined_atoms.get_primitive_atoms
        # print (atoms)
        if use_preconverged_kpoints:

            kp = Kpoints3D().automatic_length_mesh(
                # lattice_mat=atoms.lattice_mat,
                # length=10
                lattice_mat=atoms.lattice_mat,
                length=dat["kpoint_length_unit"],
            )
            kpts = kp._kpoints[0]
            kpts = non_prime_kpoints(kpts)
            kp = Kpoints3D(kpoints=[kpts])
            print("kpts", kpts)

            nq1 = get_factors(kpts[0])[0]
            nq2 = get_factors(kpts[1])[0]
            nq3 = get_factors(kpts[2])[0]
            qp = Kpoints3D(kpoints=[[nq1, nq2, nq3]])
        else:
            kp = Kpoints3D(kpoints=[])
            qp = Kpoints3D(kpoints=[])
        sup = SuperCond(atoms=atoms, kp=kp, qp=qp, qe_cmd=qe_cmd).to_dict()
        dumpjson(data=sup, filename="sup.json")
        write_qejob(job_json=os.path.abspath("sup.json"))
        path = (
            "echo hello"
            + " \n module load intel/2015\nmodule load "
            + "openmpi/2.1.0/intel-15\nsource activate qe \npython "
            + os.getcwd()
            + "/job.py"
        )

        if submit_job:
            Queue.slurm(
                job_line=path,
                jobname=i,
                directory=os.getcwd(),
                # queue="mml",
                walltime="330:0:0",
                submit_cmd=["sbatch", "submit_job"],
            )

        os.chdir("..")
    except Exception as exp:
        print(exp)
        pass
