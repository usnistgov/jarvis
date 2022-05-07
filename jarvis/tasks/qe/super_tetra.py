"""Module to run Tc calculation."""
from jarvis.io.qe.outputs import DataFileSchema
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.qe.qe import QEjob

# import numpy as np
import os

# from jarvis.core.utils import get_factors


def very_clean():
    """Clean files."""
    cmd = (
        "rm -r elph_dir _ph0 *.dos* *.in *.json  *.fc *.freq* *.save "
        + "*.out *.dyn* *wfc* *.xml  *save lambda *.modes dyn*"
    )
    os.system(cmd)


class SuperCond(object):
    """Module to calculate Tc."""

    def __init__(self, atoms=None, kp=None, qp=None, qe_cmd="pw.x"):
        """Initialize the class."""
        self.atoms = atoms
        self.kp = kp
        self.qp = qp
        self.qe_cmd = qe_cmd

    def to_dict(self):
        """Get dictionary."""
        info = {}
        info["atoms"] = self.atoms.to_dict()
        info["kp"] = self.kp.to_dict()
        info["qp"] = self.qp.to_dict()
        info["qe_cmd"] = self.qe_cmd
        return info

    @classmethod
    def from_dict(self, info={}):
        """Load from a dictionary."""
        return SuperCond(
            atoms=Atoms.from_dict(info["atoms"]),
            kp=Kpoints3D.from_dict(info["kp"]),
            qp=Kpoints3D.from_dict(info["qp"]),
            qe_cmd=info["qe_cmd"],
        )

    def runjob(self):
        """Calculate Tc using QE."""
        # Still under development,
        atoms = self.atoms
        kp = self.kp
        qp = self.qp
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
                "conv_thr": "1d-10",
            },
            "ions": {"ion_dynamics": "'bfgs'"},
            "cell": {"cell_dynamics": "'bfgs'", "cell_dofree": "'all'"},
        }
        qejob_relax = QEjob(
            atoms=atoms,
            input_params=relax,
            output_file="relax.out",
            qe_cmd=self.qe_cmd,
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
                # "degauss": 0.01,
                "nat": None,
                "ntyp": None,
                "ecutwfc": 45,
                "ecutrho": 250,
                "occupations": "'tetrahedra_opt'",
                # "smearing": "'mp'",
                # "la2F ": ".true.",
            },
            "electrons": {
                "diagonalization": "'david'",
                "mixing_mode": "'plain'",
                "mixing_beta": 0.7,
                "conv_thr": "1d-10",
            },
        }

        qejob_scf_init = QEjob(
            atoms=final_strt,
            input_params=scf_init,
            output_file="scf_init.out",
            qe_cmd=self.qe_cmd,
            jobname="scf_init",
            kpoints=kp,
            input_file="ascf_init.in",
        )

        info_scf = qejob_scf_init.runjob()
        print(info_scf)
        kpts = kp._kpoints[0]
        qpts = qp._kpoints[0]
        nq1 = qpts[0]  # get_factors(kpts[0])[0]
        nq2 = qpts[1]  # get_factors(kpts[1])[0]
        nq3 = qpts[2]  # get_factors(kpts[2])[0]
        nk1 = kpts[0]
        nk2 = kpts[1]
        nk3 = kpts[2]
        ph = {
            "inputph": {
                "prefix": "'QE'",
                "fildyn": "'QE.dyn'",
                "outdir": "'./'",
                "ldisp": ".true.",
                "lshift_q": ".true.",
                "fildvscf": "'dvscf'",
                "fildrho": "'dvrho'",
                # "electron_phonon": "'lambda_tetra'",
                # "electron_phonon": "'interpolated'",
                # "el_ph_sigma": 0.005,
                "nq1": nq1,
                "nq2": nq2,
                "nq3": nq3,
                "nk1": nk1,
                "nk2": nk2,
                "nk3": nk3,
                # "tr2_ph": "1.0d-12",
            },
        }
        qejob_ph = QEjob(
            atoms=final_strt,
            input_params=ph,
            output_file="ph.out",
            qe_cmd=self.qe_cmd.replace("pw.x", "ph.x"),
            jobname="ph",
            kpoints=None,
            input_file="aph.in",
        )

        qejob_ph.runjob()

        nq1 = qpts[0]  # get_factors(kpts[0])[0]
        nq2 = qpts[1]  # get_factors(kpts[1])[0]
        nq3 = qpts[2]  # get_factors(kpts[2])[0]
        nk1 = kpts[0]
        nk2 = kpts[1]
        nk3 = kpts[2]
        ph_tetra = {
            "inputph": {
                "prefix": "'QE'",
                "fildyn": "'QE.dyn'",
                "outdir": "'./'",
                "ldisp": ".true.",
                "lshift_q": ".true.",
                "fildvscf": "'dvscf'",
                "fildrho": "'dvrho'",
                "electron_phonon": "'lambda_tetra'",
                # "electron_phonon": "'interpolated'",
                # "el_ph_sigma": 0.005,
                "nq1": nq1,
                "nq2": nq2,
                "nq3": nq3,
                "nk1": nk1,
                "nk2": nk2,
                "nk3": nk3,
                # "tr2_ph": "1.0d-12",
            },
            "inputa2f": {
                "nfreq": 500,
            },
        }
        qejob_ph_tetra = QEjob(
            atoms=final_strt,
            input_params=ph_tetra,
            output_file="ph_tetra.out",
            qe_cmd=self.qe_cmd.replace("pw.x", "ph.x"),
            jobname="ph_tetra",
            kpoints=None,
            input_file="aph_tetra.in",
        )

        qejob_ph_tetra.runjob()
        cmd = (
            self.qe_cmd.replace("pw.x", "alpha2f.x")
            + "<"
            + "aph_tetra.in"
            + ">"
            + "aph_tetra.out"
        )
        os.system(cmd)
