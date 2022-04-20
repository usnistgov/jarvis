"""Workflow to calculate Curie Temperature for 2D magnets usong Torelli/Olsen method for VASP calculations."""
from jarvis.core.kpoints import Kpoints3D
from collections import defaultdict
from jarvis.tasks.queue_jobs import Queue
from jarvis.io.vasp.inputs import Poscar
from jarvis.tasks.vasp.vasp import VaspJob, write_vaspjob
from jarvis.io.vasp.inputs import Incar
from jarvis.db.jsonutils import dumpjson
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
from jarvis.tasks.vasp.vasp import (
    JobFactory,
    VaspJob,
    GenericIncars,
    write_jobfact,
)
import os

from jarvis.io.vasp.inputs import (
    Poscar,
    Kpoints,
    Incar,
    Potcar,
    IndividualPotcarData,
    find_ldau_magmom,
)


def get_symb(atoms):
    new_symb = []
    for i in atoms.elements:
        if i not in new_symb:
            new_symb.append(i)
    return new_symb
incs = GenericIncars().pbe().incar.to_dict()
dat = get_jid_data(jid="JVASP-76195", dataset="dft_2d")

atoms = Atoms.from_dict(dat["atoms"])
length = dat["kpoint_length_unit"]
encut = dat["encut"]

vasp_cmd = "/cluster/deb9/bin/mpirun /cluster/deb9/bin/vasp_ncl_lmpi"
copy_files = ["/users/dtw2/bin/vdw_kernel.bindat"]
def magnetic_torelli_olsen(
    mat=None,
    encut=None,
    length=20,
    u_val=0.0,
    min_configs=2,
    s_axis=[[1, 0, 0], [0, 0, 1]],
):
    from jarvis.analysis.magnetism.magmom_setup import MagneticOrdering

    count = 0
    mag = MagneticOrdering(mat)
    symm_list, ss = mag.get_minimum_configs(min_configs=min_configs)
    data = defaultdict()
    data["LCHARG"] = ".FALSE."
    data["EDIFF"] = "1E-7"
    data["LORBMOM"] = ".TRUE."
    data["ALGO"] = "All"
    data["AMIX"] = "0.01"
    data["BMIX"] = "0.0001"
    data["AMIX_MAG"] = "0.4"
    data["BMIX_MAG"] = "0.0001"
    data["PREC"] ="High"
    data["LSORBIT"] =".TRUE."
    data["LDAU"] = ".TRUE."
    data["LDAUTYPE"] ="2" #use optimal value of U for transition metal ion
    data["LDAUL"] ="2 -1"
    data["LDAUU"]="2.0 0.0"
    data["LDAUJ"]="0.0 0.0"
    data["LDAUPRINT"]="2"
    data["LMAXMIX"]="4"

    kp = Kpoints3D().automatic_length_mesh(
        length=length, lattice_mat=mat.lattice_mat
    ) #automatic kgrid
    
   # kp = Kpoints3D(kpoints=[[5,5,1]]) manually selected kgrid

    

    comment = "CrI3"
    for ii, i in enumerate(symm_list):
        # print (i)
        if "SAXIS" in data:
            data.pop("SAXIS")
        data["MAGMOM"] = " ".join(map(str, i))    
        job_name = comment + "_" + str(ii)
        chgcar_path = job_name + "/CHGCAR"
        # print(data, chgcar_path, ld)
        for j in s_axis:
            magmom = i
            saxis = j
            data["SAXIS"] = " ".join(map(str, saxis))
            data["MAGMOM"] = " ".join(
                [" ".join(map(str, [0, 0, ii])) for ii in magmom]
            )
            ld = find_ldau_magmom(atoms=mat, lsorbit=True) 
            # print(data, ld)
            # print ("incs",incs)
            data.update(incs)
            # for ii,jj in incs.items():
            #   data[ii]==jj
            inc = Incar(data)
            print("inc", inc)
            print("kp", kp)
            print("mat", mat)
            # print ("poscar",atoms)
            # print ("ld",ld)
            count += 1
            fold_name = "Calc-" + str(count)
            if not os.path.exists(fold_name):
                os.makedirs(fold_name)
            os.chdir(fold_name)
            symb = get_symb(mat)
            pot = Potcar(elements=symb)
            mat.write_poscar("POSCAR")
            incar = Incar(data)  # .write_file("INCAR")
            kp.write_file("KPOINTS")
            pot.write_file("POTCAR")
            job = VaspJob(
                poscar=Poscar(mat),
                incar=incar,
                potcar=pot,
                kpoints=kp,
                vasp_cmd=vasp_cmd,
                copy_files=copy_files,
                jobname=fold_name,
            )
            dumpjson(data=job.to_dict(), filename="job.json")
            write_vaspjob(pyname="job.py", job_json="job.json")

            directory = os.getcwd()
            tmp_name = str(fold_name) + "_submit_job"
            job_line = "source ~/miniconda3/envs/my_jarvis/bin/activate my_jarvis \npython job.py \n"
            submit_cmd = ["sbatch", tmp_name]
            Queue.slurm(
                job_line=job_line,
                filename=tmp_name,
                jobname=fold_name,
                nnodes=1,  # rtx
                cores=16,  # rtx
                directory=directory,
                submit_cmd=submit_cmd,
                queue="fast",
                walltime="30:00:00",
            )

            os.chdir("..")


magnetic_torelli_olsen(mat=atoms, length=length)
