"""Module to run Tc calculation."""
from jarvis.io.qe.outputs import DataFileSchema
from jarvis.tasks.qe.qe import QEjob
import numpy as np


def calc_Tc(wlog=300, lamb=1.0, mu=0.1):
    """Calculate Tc."""
    tc = (wlog / 1.2) * np.exp(
        -1.04 * (1 + lamb) / (lamb * (1 - 0.062 * mu) - mu)
    )
    return tc


def parse_lambda(filename="lambda"):
    """Parse lambda file."""
    f = open(filename, "r")
    lines = f.read().splitlines()
    f.close()
    for i in lines:
        if "Broadening" in i:
            tmp = i.split()
            print(i.split())
            wlog = float(tmp[-1])
            lamb = float(tmp[3])
            Tc = calc_Tc(wlog=wlog, lamb=lamb)
            print("Tc", Tc)
            print()


def get_factors(x, start=2):
    """Get factors of a number."""
    facts = []
    for i in range(start, x + 1):
        if x % i == 0:
            facts.append(i)
    return facts


def supercond_workflow(atoms=None, kp=None):
    """Calculate Tc using QE."""
    # Still under development,
    relax = {
        "control": {
            # "calculation": "'scf'",
            "calculation": "'vc-relax'",
            "restart_mode": "'from_scratch'",
            "prefix": "'RELAX'",
            "outdir": "'./'",
            "tstress": ".true.",
            "tprnfor": ".true.",
            "disk_io": "'low'",
            "wf_collect": ".true.",
            "pseudo_dir": None,
            "verbosity": "'high'",
            "nstep": 100,
        },
        "system": {
            "ibrav": 0,
            "nat": None,
            "ntyp": None,
            "ecutwfc": 45,
            "ecutrho": 250,
            "q2sigma": 1,
            "ecfixed": 44.5,
            "qcutz": 800,
            "occupations": "'smearing'",
            "degauss": 0.01,
            "lda_plus_u": ".false.",
        },
        "electrons": {
            "diagonalization": "'david'",
            "mixing_mode": "'local-TF'",
            "mixing_beta": 0.3,
            "conv_thr": "1d-9",
        },
        "ions": {"ion_dynamics": "'bfgs'"},
        "cell": {"cell_dynamics": "'bfgs'", "cell_dofree": "'all'"},
    }
    qejob_relax = QEjob(
        atoms=atoms,
        input_params=relax,
        output_file="relax.out",
        qe_cmd="/home/knc6/Software/qe/q-e/bin/pw.x",
        jobname="relax",
        kpoints=kp,
        input_file="arelax.in",
    )

    info = qejob_relax.runjob()
    xml_path = info["xml_path"]
    xml_dat = DataFileSchema(xml_path)
    final_strt = xml_dat.final_structure

    print("final_strt", final_strt)

    # scf_init =
    scf_init = {
        "control": {
            "calculation": "'scf'",
            "restart_mode": "'from_scratch'",
            "prefix": "'QE'",
            "outdir": "'./'",
            "tstress": ".true.",
            "tprnfor": ".true.",
            "disk_io": "'low'",
            "pseudo_dir": None,
            "verbosity": "'high'",
            "nstep": 100,
            "etot_conv_thr": "1.0d-5",
            "forc_conv_thr": "1.0d-4",
        },
        "system": {
            "ibrav": 0,
            "degauss": 0.01,
            "nat": None,
            "ntyp": None,
            "ecutwfc": 45,
            "ecutrho": 250,
            "occupations": "'smearing'",
            "smearing": "'mp'",
            "la2F ": ".true.",
        },
        "electrons": {
            "diagonalization": "'david'",
            "mixing_mode": "'plain'",
            "mixing_beta": 0.7,
            "conv_thr": "1d-9",
        },
    }

    qejob_scf_init = QEjob(
        atoms=final_strt,
        input_params=scf_init,
        output_file="scf_init.out",
        qe_cmd="/home/knc6/Software/qe/q-e/bin/pw.x",
        jobname="scf_init",
        kpoints=kp,
        input_file="ascf_init.in",
    )

    info_scf = qejob_scf_init.runjob()
    print(info_scf)
    kpts = kp._kpoints[0]
    nq1 = get_factors(kpts[0])[0]
    nq2 = get_factors(kpts[0])[0]
    nq3 = get_factors(kpts[0])[0]
    ph = {
        "inputph": {
            "prefix": "'QE'",
            "fildyn": "'QE.dyn'",
            "outdir": "'./'",
            "ldisp": ".true.",
            "trans": ".true.",
            "fildvscf": "'dvscf'",
            "electron_phonon": "'interpolated'",
            "nq1": nq1,
            "nq2": nq2,
            "nq3": nq3,
            "tr2_ph": "1.0d-12",
        }
    }
    qejob_ph = QEjob(
        atoms=final_strt,
        input_params=ph,
        output_file="ph.out",
        qe_cmd="/home/knc6/Software/qe/q-e/bin/ph.x",
        jobname="ph",
        kpoints=None,
        input_file="aph.in",
    )

    qejob_ph.runjob()
    # import sys
    # sys.exit()
    qr = {
        "input": {
            "zasr": "'simple'",
            "fildyn": "'QE.dyn'",
            "flfrc": "'QE333.fc'",
            "la2F": ".true.",
        }
    }
    qejob_qr = QEjob(
        atoms=final_strt,
        input_params=qr,
        output_file="q2r.out",
        qe_cmd="/home/knc6/Software/qe/q-e/bin/q2r.x",
        jobname="qr",
        kpoints=None,
        input_file="aqr.in",
    )

    qejob_qr.runjob()

    matdyn = {
        "input": {
            "asr": "'simple'",
            "flfrc": "'QE333.fc'",
            "flfrq": "'QE333.freq'",
            "la2F": ".true.",
            "dos": ".true.",
            "fldos": "'phonon.dos'",
            "nk1": 10,
            "nk2": 10,
            "nk3": 10,
            "ndos": 50,
        }
    }

    qejob_matdyn = QEjob(
        atoms=final_strt,
        input_params=matdyn,
        output_file="matdyn.out",
        qe_cmd="/home/knc6/Software/qe/q-e/bin/matdyn.x",
        jobname="matdyn",
        kpoints=None,
        input_file="amatdyn.in",
    )

    qejob_matdyn.runjob()
    parse_lambda()


# """
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    from jarvis.db.figshare import get_jid_data
    from jarvis.core.kpoints import Kpoints3D

    dat = get_jid_data(jid="JVASP-4406", dataset="dft_3d")
    # dat = get_jid_data(jid="JVASP-19821", dataset="dft_3d")
    # dat = get_jid_data(jid="JVASP-816", dataset="dft_3d")
    atoms = Atoms.from_dict(dat["atoms"])
    kp = Kpoints3D().automatic_length_mesh(
        lattice_mat=atoms.lattice_mat, length=10
    )
    supercond_workflow(atoms=atoms, kp=kp)
# """
