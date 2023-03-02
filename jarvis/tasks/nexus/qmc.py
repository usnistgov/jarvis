"""Module for running QMCPack."""
from jarvis.io.nexus.inputs import (
    get_Zeff,
    get_pseudo_dft_dict,
    # get_pseudo_qmc_dict,
)
from jarvis.io.qe.outputs import QEout
from nexus import generate_physical_system
from nexus import generate_pwscf
from nexus import settings, job, run_project, obj
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
from nexus import read_structure
from jarvis.core.kpoints import Kpoints3D
import numpy as np
import os


def get_energy_nexus(
    atoms=None,
    poscar_path=None,
    length=10,
    ecut=40,
    temp_filename="POSCAR-temp.vasp",
    qe_presub="module load qe/7.0.pw2qmcpack",
    qmcpack_presub="module load qmcpack/3.14.0",
    nodes=2,
    cores=16,
    threads=1,
    app="pw.x",
    net_spin=0,
    path="scf",
    machine="dobby",
    pseudo_path="/rk2/knc6/QMC/pseudopotentials",
    results_path="/users/knc6/Software/jarvis/usnistgov/jarvis/results",
    runs_path="/users/knc6/Software/jarvis/usnistgov/jarvis/runs",
):
    """Get energy from Nexus+qmcpack."""
    results = results_path
    runs = runs_path
    settings(
        results=results,
        pseudo_dir=pseudo_path,  # location of pseudopotential directory
        sleep=1,
        runs=runs,
        machine=machine,
    )
    if atoms is None:
        atoms = Atoms.from_poscar(poscar_path)
        temp_filename = poscar_path

    kp = (
        Kpoints3D()
        .automatic_length_mesh(lattice_mat=atoms.lattice_mat, length=length)
        ._kpoints
    )
    kp = np.array(kp).flatten().tolist()
    atoms.write_poscar(temp_filename)
    structure = read_structure(temp_filename, format="poscar")
    print("structure", structure)
    print("Running job")
    shared_qe = obj(
        occupations="smearing",
        smearing="gaussian",
        degauss=0.005,
        input_dft="PBE",
        ecut=ecut,  # Plane Wave cutoff energy
        conv_thr=1.0e-7,
        mixing_beta=0.2,
        nosym=True,
        use_folded=True,
        spin_polarized=True,
    )

    qe_job = job(
        nodes=nodes, cores=cores, threads=threads, app=app, presub=qe_presub
    )
    Zeff = get_Zeff()
    system = generate_physical_system(
        structure=structure,
        kshift=(0, 0, 0),
        net_spin=net_spin,  # Specify Up - Down Electrons
        **{e: Zeff[e] for e in structure.elem}
    )
    pseudo_dft_dict = get_pseudo_dft_dict()
    arr = []
    for sp in atoms.elements:
        arr.append(pseudo_dft_dict[sp])
    scf = generate_pwscf(
        identifier="scf",  # log output goes to scf.out
        path="scf",  # directory to run in
        job=qe_job,  # pyscf must run w/o mpi
        system=system,
        input_type="scf",
        pseudos=arr,
        kgrid=kp,
        wf_collect=False,
        **shared_qe
    )
    print('scf', scf)
    run_project()
    out_file = os.path.join(runs, "scf", "scf.out")
    print("out_file", out_file)
    qeout = QEout(filename=out_file)
    tot_energy = qeout.get_total_energy()
    print("tot_energy", tot_energy)
    return tot_energy


if __name__ == "__main__":
    jid = "JVASP-32"
    atoms = Atoms.from_dict(
        get_jid_data(jid="JVASP-1002", dataset="dft_3d")["atoms"]
    )
    get_energy_nexus(atoms=atoms)
