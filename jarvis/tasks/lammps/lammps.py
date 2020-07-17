"""Modules for running LAMMPS calculations."""

from jarvis.analysis.structure.spacegroup import (
    Spacegroup3D,
    symmetrically_distinct_miller_indices,
)
from jarvis.core.atoms import Atoms
from jarvis.io.lammps.inputs import LammpsInput, LammpsData
from jarvis.tasks.lammps.templates.templates import GenericInputs
from jarvis.io.lammps.outputs import analyze_log
from jarvis.analysis.defects.vacancy import Vacancy
from jarvis.analysis.defects.surface import Surface
import shutil
import os
import subprocess
import json
import sys
from collections import OrderedDict


class JobFactory(object):
    """Class for generic LAMMPS calculations."""

    def __init__(self, name="", pair_style="", pair_coeff="", control_file=""):
        """
        Use in defining a LAMMPS job.

        With following arguments.
        Args:
            pair_style :  LAMMPS pair_style, e.g. "eam/alloy"

            pair_coeff : path for pair-coefficients file

            control_file :  control-file with units, include modules
                          for running LAMMPS calculation , see examples

            name : generic name
        """
        self.pair_style = pair_style
        self.pair_coeff = pair_coeff
        self.control_file = control_file
        self.name = name

    def all_props_eam_alloy(
        self,
        atoms=None,
        ff_path="",
        lammps_cmd="",
        enforce_conventional_structure=True,
        enforce_c_size=0,
        extend=1,
    ):
        """
        Provide generic function for LAMMPS calculations using eam/alloy.

        Must provide Atoms class and path to force-field.
        Args:
            atoms :  Atoms object

            ff_path :  inter-atomic potential path

            lammps_cmd : LAMMPS executable path

            enforce_conventional_structure :
            whether to enforce conventional cell

            enforce_c_size : minimum cell-sizes

            extend : used for round-off during making supercells
        """
        if enforce_conventional_structure:
            atoms = Spacegroup3D(atoms).conventional_standard_structure

        a = atoms.lattice.lat_lengths()[0]
        b = atoms.lattice.lat_lengths()[1]
        c = atoms.lattice.lat_lengths()[2]
        if enforce_c_size is not None:
            dim1 = int(float(enforce_c_size) / float(a)) + extend
            dim2 = int(float(enforce_c_size) / float(b)) + extend
            dim3 = int(float(enforce_c_size) / float(c)) + extend
            atoms = atoms.make_supercell([dim1, dim2, dim3])

        self.pair_style = "eam/alloy"
        self.pair_coeff = ff_path
        parameters = {
            "pair_style": self.pair_style,
            "atom_style": "charge",
            "pair_coeff": self.pair_coeff,
        }
        parameters["control_file"] = "inelast.mod"
        en, final_str, forces = LammpsJob(
            atoms=atoms,
            jobname="ELASTIC",
            parameters=parameters,
            lammps_cmd=lammps_cmd,
        ).runjob()
        print("en, final_str, forces", en, final_str, forces)

        indices = symmetrically_distinct_miller_indices(
            max_index=1, cvn_atoms=atoms
        )
        for i in indices:
            surf = Surface(atoms=final_str, indices=i).make_surface()
            jobname = str("Surf-") + str("_".join(map(str, i)))
            en2, final_str2, forces2 = LammpsJob(
                atoms=surf,
                jobname=jobname,
                parameters=parameters,
                lammps_cmd=lammps_cmd,
            ).runjob()

        sys.exit()

        v = Vacancy(atoms=final_str).generate_defects(enforce_c_size=5)
        print("vacs=", v)
        for i, ii in enumerate(v):
            jobname = (
                str("symbol-")
                + str(ii._symbol)
                + str("-")
                + str("Wycoff-")
                + str(ii._wyckoff_multiplicity)
            )
            print("ii._defect_structure", ii._atoms)
            en2, final_str2, forces2 = LammpsJob(
                atoms=ii._defect_structure,
                jobname=jobname,
                parameters=parameters,
                lammps_cmd=lammps_cmd,
            ).runjob()

    # def optimize_and_elastic(self):
    #     pass

    # def surface_energy(self):
    #     pass

    # def vacancy(self):
    #     pass

    # def phonon(self):
    #     pass


class LammpsJob(object):
    """Construct a class representing a LAMMPS job."""

    def __init__(
        self,
        atoms=None,
        parameters={
            "pair_style": "eam/alloy",
            "pair_coeff": "abc.alloy",
            "atom_style": "charge",
            "control_file": "inelast.mod",
        },
        lammps_cmd=None,
        output_file="lammps.out",
        stderr_file="std_err.txt",
        jobname="testt",
        attempts=5,
        copy_files=[],
        element_order=[],
    ):
        """
        Use for defining a LAMMPS job.

        Provide follwoing arguments.

        Args:
            atoms :  Atoms object

            element_order :
            element order used in accessing force-field parameters

            parameters :  LAMMPS input parameter dictionary

            lammps_cmd : LAMMPS executable path

            output_file :  standard output file

            stderr_file :  standard error file

            jobname :  Job name

            attempts :  number of attempts before crashing the job, TODO

            copy_files : copy certain files before a job
        """
        self.atoms = atoms
        self.element_order = element_order
        self.parameters = parameters
        self.lammps_cmd = lammps_cmd
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.jobname = jobname
        self.attempts = attempts
        self.copy_files = copy_files

    def to_dict(self):
        """Convert to a dictionary."""
        d = OrderedDict()
        d["atoms"] = self.atoms.to_dict()
        d["element_order"] = self.element_order
        d["parameters"] = self.parameters
        d["lammps_cmd"] = self.lammps_cmd
        d["output_file"] = self.output_file
        d["stderr_file"] = self.stderr_file
        d["jobname"] = self.jobname
        d["attempts"] = self.attempts
        d["copy_files"] = self.copy_files
        return d

    @classmethod
    def from_dict(self, d={}):
        """Load from a dictionary."""
        return LammpsJob(
            atoms=Atoms.from_dict(d["atoms"]),
            element_order=self.element_order,
            parameters=d["parameters"],
            lammps_cmd=d["lammps_cmd"],
            output_file=d["output_file"],
            stderr_file=d["stderr_file"],
            jobname=d["jobname"],
            attempts=d["attempts"],
            copy_files=d["copy_files"],
        )

    def write_input(self):
        """Write LAMMPS input files."""
        lmp = LammpsData().atoms_to_lammps(atoms=self.atoms)
        self.element_order = lmp._element_order
        lmp.write_file("data")
        LammpsInput(LammpsDataObj=lmp).write_lammps_in(
            parameters=self.parameters
        )
        for i in self.copy_files:
            shutil.copy2(i, ".")
        if "control_file" in self.parameters:
            if self.parameters["control_file"] == "inelast.mod":
                GenericInputs().elastic_general(path=".")

    def run(self):
        """Run a job with subprocess."""
        with open(self.output_file, "w") as f_std, open(
            self.stderr_file, "w", buffering=1
        ) as f_err:
            p = subprocess.Popen(
                self.lammps_cmd, shell=True, stdout=f_std, stderr=f_err
            )
            p.wait()
        return p

    def runjob(self):
        """Constrct  a generic LAMMPS job submission."""
        attempt = 0
        wait = False
        while not wait:
            attempt = attempt + 1
            if attempt == self.attempts:
                wait = True
            folder = str(os.getcwd()) + str("/") + str(self.jobname)
            if not os.path.exists(folder):
                os.makedirs(folder)
            os.chdir(folder)

            if os.path.isfile("./log.lammps"):
                try:

                    (
                        en,
                        press,
                        toten,
                        c11,
                        c22,
                        c33,
                        c12,
                        c13,
                        c23,
                        c44,
                        c55,
                        c66,
                        c14,
                        c16,
                        c24,
                        c25,
                        c26,
                        c34,
                        c35,
                        c36,
                        c45,
                        c46,
                        c56,
                    ) = analyze_log("./log.lammps")
                    wait = True
                except Exception:
                    pass
                # print ('toten',toten)
            else:
                self.write_input()
                self.run()
                try:

                    (
                        en,
                        press,
                        toten,
                        c11,
                        c22,
                        c33,
                        c12,
                        c13,
                        c23,
                        c44,
                        c55,
                        c66,
                        c14,
                        c16,
                        c24,
                        c25,
                        c26,
                        c34,
                        c35,
                        c36,
                        c45,
                        c46,
                        c56,
                    ) = analyze_log("./log.lammps")
                    wait = True
                except Exception:
                    pass
            pot = os.path.join(os.getcwd(), "potential.mod")
            # print ('toten2',toten,pot)
            initial_str = LammpsData().read_data(
                filename="data",
                element_order=self.element_order,
                potential_file=pot,
            )
            final_str = LammpsData().read_data(
                potential_file=pot,
                filename="data0",
                element_order=self.element_order,
            )
            forces = []

            data_cal = []
            data_cal.append(
                {
                    "jobname": self.jobname,
                    "initial_pos": initial_str.to_dict(),
                    "pair_style": str(self.parameters["pair_style"]),
                    "pair_coeff": str(self.parameters["pair_coeff"]),
                    "final_energy": float(toten),
                    "en": en,
                    "press": press,
                    "final_str": final_str.to_dict(),
                }
            )

            json_file = str(self.jobname) + str(".json")
            os.chdir("../")
            f_json = open(json_file, "w")
            f_json.write(json.dumps(data_cal))
            f_json.close()
            return en, final_str, forces
