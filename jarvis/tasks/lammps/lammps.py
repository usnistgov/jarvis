from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
from jarvis.io.lammps.inputs import LammpsInput, LammpsData
from jarvis.tasks.lammps.templates.templates  import LammpsJobFactory
from jarvis.io.lammps.outputs import analyze_log
import os
import subprocess
import json

class LammpsJob(object):
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

        self.atoms = atoms
        self.element_order = element_order
        self.parameters = parameters
        self.lammps_cmd = lammps_cmd
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.jobname = jobname
        self.attempts = attempts
        self.copy_files = copy_files

    def write_input(self):
        lmp = LammpsData().atoms_to_lammps(atoms=self.atoms)
        self.element_order=lmp._element_order 
        lmp.write_file('data')
        LammpsInput(LammpsDataObj=lmp).write_lammps_in(parameters=self.parameters)
        for i in self.copy_files:
            shutil.copy2(i, ".")
        if "control_file" in self.parameters:
             if self.parameters["control_file"]=="inelast.mod":
                    LammpsJobFactory().elastic_general(path=".")

    def run(self):
        with open(self.output_file, "w") as f_std, open(
            self.stderr_file, "w", buffering=1
        ) as f_err:
            p = subprocess.Popen(
                self.lammps_cmd, shell=True, stdout=f_std, stderr=f_err
            )
            p.wait()
        return p

    def runjob(self):
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
            except:
                pass
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
            except:
                pass
            initial_str = LammpsData().read_data(
                filename="data", element_order=self.element_order,potential_file='potential.mod'
            )
            final_str = LammpsData().read_data(
                filename="data0", element_order=self.element_order
            )
            forces = []

            data_cal = []
            data_cal.append(
                {
                    "jobname": self.jobname,
                    "initial_pos": initial_str.to_dict(),
                    "pair_style": str(parameters["pair_style"]),
                    "pair_coeff": str(parameters["pair_coeff"]),
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


if __name__ == "__main__":
    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-FF/COMB/ffield.comb3.NiAlO_nist/bulk@mp-1143_fold/bulk@mp-1143/new_pymatgen_slab.vasp"
    )
    atoms = Poscar.from_file(
        "/rk2/knc6/JARVIS-FF/FS/Al1.eam.fs_nist/bulk@mp-134_fold/mp-134/new_pymatgen_slab.vasp"
    ).atoms
    from jarvis.analysis.structure.spacegroup import Spacegroup3D
    cvn=Spacegroup3D(atoms).conventional_standard_structure
    parameters={
            "pair_style": "eam/alloy",
            "pair_coeff": "/data/knc6/JARVIS-FF-NEW/FS/Al1.eam.fs",
            "atom_style": "charge",
            "control_file": "inelast.mod",
        }

    cmd='/users/knc6/Software/LAMMPS/lammps-master/src/lmp_serial<in.main'
    LammpsJob(atoms=cvn,parameters=parameters,lammps_cmd=cmd).runjob()

    lmp = LammpsData().atoms_to_lammps(atoms=cvn)
    LammpsInput(LammpsDataObj=lmp).write_lammps_in(
        parameters={
            "pair_style": "eam/alloy",
            "pair_coeff": "/data/knc6/JARVIS-FF-NEW/FS/Al1.eam.fs",
            "atom_style": "charge",
            "control_file": "inelast.mod",
        }
    )
