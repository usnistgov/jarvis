from jarvis.analysis.defects.vacancy import Vacancy
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.db.jsonutils import dumpjson
from jarvis.core.atoms import Atoms
from jarvis.tasks.vasp.vasp import VaspJob, write_vaspjob
from jarvis.tasks.queue_jobs import Queue
from jarvis.io.vasp.inputs import Poscar
from jarvis.db.figshare import data
from jarvis.io.vasp.inputs import Poscar, Incar, Potcar
from jarvis.core.kpoints import Kpoints3D

dat = data("dft_3d")
import os

data0 = dict(
    PREC="Accurate",
    ISMEAR=0,
    GGA="BO",
    PARAM1=0.1833333333,
    PARAM2=0.2200000000,
    LUSE_VDW=".TRUE.",
    AGGAC=0.0000,
    EDIFF="1E-7",
    NELM=400,
    ISPIN=2,
    # NPAR=32,
    POTIM=0.02,
    LCHARG=".FALSE.",
    LVTOT=".FALSE.",
    LVHAR=".FALSE.",
    LWAVE=".TRUE.",
)

data1 = dict(
    PREC="Accurate",
    ISMEAR=0,
    IBRION=6,
    GGA="BO",
    PARAM1=0.1833333333,
    PARAM2=0.2200000000,
    LUSE_VDW=".TRUE.",
    AGGAC=0.0000,
    EDIFF="1E-7",
    NELM=400,
    ISPIN=2,
    # NPAR=32,
    POTIM=0.02,
    LCHARG=".FALSE.",
    LVTOT=".FALSE.",
    LVHAR=".FALSE.",
    LWAVE=".FALSE.",
)

data2 = dict(
    PREC="Accurate",
    ISMEAR=-2,
    IBRION=6,
    GGA="BO",
    PARAM1=0.1833333333,
    PARAM2=0.2200000000,
    LUSE_VDW=".TRUE.",
    AGGAC=0.0000,
    EDIFF="1E-7",
    NELM=400,
    ISPIN=2,
    # NPAR=32,
    POTIM=0.02,
    LCHARG=".FALSE.",
    LVTOT=".FALSE.",
    LVHAR=".FALSE.",
    LWAVE=".FALSE.",
)

inc0 = Incar(data0)
inc1 = Incar(data1)
inc2 = Incar(data2)


def get_symb(atoms):
    new_symb = []
    for i in atoms.elements:
        if i not in new_symb:
            new_symb.append(i)
    return new_symb


# vasp_cmd = "/usr/bin/mpirun /usr/local/src/vasp-5.4.4/build/std/vasp"
# vasp_cmd = "/usr/bin/mpirun /users/knc6/VASP/vasp544/vasp.5.4.4/bin/vasp_std"
# copy_files = ["/users/knc6/bin/vdw_kernel.bindat"]
jids = ["JVASP-19821"]  # ,'JVASP-14615']

jids = [
    "JVASP-93855",
    "JVASP-19821",
    "JVASP-19668",
    "JVASP-19684",
    "JVASP-54933",
]
jids = ["JVASP-1002"]  # ,'JVASP-14615']


def supercon_igor(
    id="JVASP-1002",
    atoms=None,
    copy_files=["/users/knc6/bin/vdw_kernel.bindat"],
    vasp_cmd="/usr/bin/mpirun /users/knc6/VASP/vasp544/vasp.5.4.4/bin/vasp_std",
):
    if "JVASP-" in id:
        for i in dat:
            if i["jid"] == id:
                atoms = Atoms.from_dict(i["atoms"])

    strt = Spacegroup3D(atoms).conventional_standard_structure
    encut = i["encut"]

    # STEP1: Generate WAVECAR

    name = str(id) + "_RUN0"
    wavecar_path = os.getcwd() + "/" + name + "/" + name + "/WAVECAR"
    print("WAVECAR path", wavecar_path)
    pos = Poscar(strt)
    symb = get_symb(strt)
    inc0.update({"ENCUT": i["encut"]})
    pos.comment = name
    pot = Potcar(elements=symb)
    kp = (
        Kpoints3D()
        .automatic_length_mesh(
            lattice_mat=strt.lattice_mat, length=i["kpoint_length_unit"] - 25
        )
        ._kpoints
    )
    kpoints = Kpoints3D(kpoints=kp)

    jobname = name

    job1 = VaspJob(
        poscar=pos,
        incar=inc0,
        potcar=pot,
        kpoints=kpoints,
        vasp_cmd=vasp_cmd,
        copy_files=copy_files,
        jobname=jobname,
    )

    name = str(id) + "_SCREENED"
    # screened
    pos = Poscar(strt)
    symb = get_symb(strt)
    inc1.update({"ENCUT": i["encut"]})
    pos.comment = name
    pot = Potcar(elements=symb)
    kp = (
        Kpoints3D()
        .automatic_length_mesh(
            lattice_mat=strt.lattice_mat, length=i["kpoint_length_unit"] - 25
        )
        ._kpoints
    )
    kpoints = Kpoints3D(kpoints=kp)

    jobname = name
    job2 = VaspJob(
        poscar=pos,
        incar=inc1,
        potcar=pot,
        kpoints=kpoints,
        vasp_cmd=vasp_cmd,
        copy_files=copy_files,
        jobname=jobname,
    )

    # break
    name = str(id) + "_UNSCREENED"
    # if not os.path.exists(name):
    #    os.makedirs(name)
    # os.chdir(name)
    # Unscreened
    pos = Poscar(strt)
    symb = get_symb(strt)
    inc2.update({"ENCUT": i["encut"]})
    pos.comment = name
    pot = Potcar(elements=symb)
    kp = (
        Kpoints3D()
        .automatic_length_mesh(
            lattice_mat=strt.lattice_mat, length=i["kpoint_length_unit"] - 25
        )
        ._kpoints
    )
    kpoints = Kpoints3D(kpoints=kp)

    jobname = name

    # copy_files.append(wavecar_path)
    job3 = VaspJob(
        poscar=pos,
        incar=inc2,
        potcar=pot,
        kpoints=kpoints,
        vasp_cmd=vasp_cmd,
        copy_files=copy_files,
        jobname=jobname,
    )

    jobs = [job1, job2, job3]
    for ii, i in enumerate(jobs):
        if ii > 0:

            copy_files.append(wavecar_path)
            i.copy_files = copy_files
        name = i.jobname
        if not os.path.exists(name):
            os.makedirs(name)
        os.chdir(name)
        i.runjob()
        """
        job_line = "source ~/anaconda2/envs/my_jarvis/bin/activate my_jarvis \npython job.py \n"
        directory = os.getcwd()

        submit_cmd = ["qsub", "-q", "kamal", "submit_job"]
        Queue.pbs(
            job_line=job_line,
            jobname=name,
            directory=directory,
            cores=32,
            submit_cmd=submit_cmd,
        )
        """
        os.chdir("..")
