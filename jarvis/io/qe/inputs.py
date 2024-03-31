"""Modules for generating input files for QE."""

import os
import requests
import tarfile
from jarvis.core.specie import Specie
import numpy as np

# from jarvis.analysis.structure.spacegroup import Spacegroup3D


class QEinfile(object):
    """
    Gnererate quantum espresso input files.

    jarvis.core.Atoms and Kpoints are required as inputs.
    For other input parameters,
    see https://www.quantum-espresso.org/Doc/INPUT_PW.html .
    """

    def __init__(
        self,
        atoms=None,
        kpoints=None,
        psp_dir=None,
        input_params={},
        url=None,
        sanitize=True,
        sanitize_tol=2e-4,
        psp_temp_name=None,
    ):
        """Initialize input parameters for qunatum espresso."""
        if psp_temp_name is None:
            psp_temp_name = "QE_PSPs"
        if input_params == {}:
            input_params = GenericInputs().geometry_optimization()
        if psp_dir is None:
            psp_dir = str(
                os.path.join(os.path.dirname(__file__), psp_temp_name)
            )
            # Download GBRV PSPs by default
            if url is None:
                url = (
                    "https://www.physics.rutgers.edu/"
                    "gbrv/all_pbesol_UPF_v1.5.tar.gz"
                )
                print("Please cite for PSPs:")
                print("https://doi.org/10.1016/j.commatsci.2013.08.053")
            if not os.path.exists(psp_dir):
                print("Downloading PSPs")
                tar_file_name = str(
                    os.path.join(os.path.dirname(__file__), "psp.tar.gz")
                )
                r = requests.get(url)
                f = open(tar_file_name, "wb")
                f.write(r.content)
                f.close()
                tar = tarfile.open(tar_file_name)
                tar.extractall(psp_dir)
                tar.close()
                os.remove(tar_file_name)
        self.species = []
        self.psp_dir = psp_dir
        self.input_params = input_params
        self.atoms = atoms
        self.kpoints = kpoints
        self.sanitize = sanitize
        self.sanitize_tol = sanitize_tol
        if "system" in input_params:
            self.system_params = input_params["system"]
            if (
                input_params["system"]["nat"] is None
                or input_params["system"]["nat"] == ""
            ):
                input_params["system"]["nat"] = self.atoms.num_atoms
            if (
                input_params["system"]["ntyp"] is None
                or input_params["system"]["ntyp"] == ""
            ):
                self.species = self.atoms.uniq_species
                input_params["system"]["ntyp"] = len(self.species)
            if (
                "nspin" in input_params["system"]
                and input_params["system"]["nspin"] == 2
            ):

                for ii in range(input_params["system"]["nat"]):
                    tmp = "starting_magnetization(" + str(ii + 1) + str(")")
                    input_params["system"][tmp] = 1.0

                # 'starting_magnetization' in input_params["system"] and
                # input_params["system"]['starting_magnetization'] is None:

        else:
            self.system_params = {}

        if "electrons" in input_params:
            self.electron_params = input_params["electrons"]
        else:
            self.electron_params = {}

        if "control" in input_params:
            self.control_params = input_params["control"]
            self.control_params["pseudo_dir"] = str("'") + self.psp_dir + "'"
            if (
                "prefix" in self.control_params
                and self.control_params["prefix"] is None
            ):
                self.control_params[
                    "prefix"
                ] = self.atoms.composition.reduced_formula
        else:
            self.control_params = {}

        if "inputa2f" in input_params:
            self.inputa2f = input_params["inputa2f"]
        else:
            self.inputa2f = {}

        if "ions" in input_params:
            self.ion_params = input_params["ions"]
        else:
            self.ion_params = {}

        if "cell" in input_params:
            self.cell_params = input_params["cell"]
        else:
            self.cell_params = {}
        if "input" in input_params:
            self.input = input_params["input"]
            if (
                "amass(1)" not in input_params["input"]
                and "zasr" not in input_params["input"]
            ):
                for ii, jj in enumerate(self.atoms.uniq_species):
                    tmp = "amass(" + str(ii + 1) + ")"
                    input_params["input"][tmp] = str(
                        round(Specie(jj).atomic_mass, 2)
                    )

        else:
            self.input = {}

        if "inputph" in input_params:
            self.inputph = input_params["inputph"]
            if "amass(1)" not in input_params["inputph"]:
                for ii, jj in enumerate(self.atoms.uniq_species):
                    tmp = "amass(" + str(ii + 1) + ")"
                    input_params["inputph"][tmp] = str(
                        round(Specie(jj).atomic_mass, 2)
                    )
        else:
            self.inputph = {}

    def dictionary_to_string(self, tags={}):
        """Convert a dictionary to string with '=' sign."""
        lines = ""
        for i, j in tags.items():
            if j is not None:
                lines = lines + str(i) + str(" = ") + str(j) + "\n"
        return lines

    def kpoints_to_string(self):
        """Convert a jarvis.core.Kpoints3D to string."""
        kp = ""
        if self.kpoints:
            kpoint_mode = self.kpoints._kpoint_mode
            if kpoint_mode == "automatic":
                kp = (
                    kp
                    + "K_POINTS automatic\n"
                    + (" ".join(map(str, self.kpoints.kpts[0])) + " 0 0 0\n")
                )
            elif kpoint_mode == "linemode":
                points = ""
                for i in self.kpoints.kpts:
                    points = points + " ".join(map(str, i)) + " 1.0" + "\n"
                kp = (
                    kp
                    + "K_POINTS crystal\n"
                    + str(len(self.kpoints.kpts))
                    + "\n"
                    + points
                )

            else:
                print(
                    "Kpoint scheme not implemented except linemode&automatic"
                )
        return kp

    def get_psp(self, element):
        """Obtain psuedopotential for an element."""
        element = str(element).lower()
        for i in os.listdir(self.psp_dir):
            el = str(i.split(".")[0].split("_")[0]).lower()
            if el == element:
                return i

    def atomic_species_string(self):
        """Obtain string for QE atomic species."""
        line = ""
        for i in self.species:
            line = line + (
                i
                + " "
                + str(round(Specie(i).atomic_mass, 2))
                + " "
                + str(self.get_psp(i))
                + "\n"
            )
        return line

    def check_frac(self, n):
        """Check fractional coordinates or lattice parameters."""
        items = [
            0.0,
            0.3333333333333333,
            0.25,
            0.5,
            0.75,
            0.6666666666666667,
            1.0,
            1.5,
            2.0,
            -0.5,
            -2.0,
            -1.5,
            -1.0,
            1.0 / 2.0 ** 0.5,
            -1.0 / 2.0 ** 0.5,
            3.0 ** 0.5 / 2.0,
            -(3.0 ** 0.5) / 2.0,
            1.0 / 3.0 ** 0.5,
            -1.0 / 3.0 ** 0.5,
            1.0 / 2.0 / 3 ** 0.5,
            -1.0 / 2.0 / 3 ** 0.5,
            1 / 6,
            5 / 6,
        ]
        items = items + [(-1) * i for i in items]
        for f in items:
            if abs(f - n) < self.sanitize_tol:
                return f
        return n

    def atomic_pos(self):
        """Obtain string for QE atomic positions."""
        line = ""
        # self.atoms = Spacegroup3D(self.atoms).refined_atoms
        coords = np.array(self.atoms.frac_coords)
        ntot = self.atoms.num_atoms

        if self.sanitize:
            for i in range(ntot):
                for j in range(3):  # neatin
                    coords[i, j] = self.check_frac(coords[i, j])

        for i, j in zip(self.atoms.elements, coords):
            line = line + str(i) + " " + " ".join(map(str, j)) + "\n"
        return line

    def atomic_cell_params(self):
        """Obtain string for QE atomic lattice parameters."""
        lat_mat = np.array(self.atoms.lattice_mat)
        if self.sanitize:
            print("Sanitizing Atoms.")
            a_lat = np.linalg.norm(lat_mat[0, :])
            at = lat_mat / a_lat
            for i in range(3):
                for j in range(3):
                    at[i, j] = self.check_frac(at[i, j])

            lat_mat = at * a_lat

        line = (
            str(lat_mat[0][0])
            + " "
            + str(lat_mat[0][1])
            + " "
            + str(lat_mat[0][2])
            + "\n"
            + str(lat_mat[1][0])
            + " "
            + str(lat_mat[1][1])
            + " "
            + str(lat_mat[1][2])
            + "\n"
            + str(lat_mat[2][0])
            + " "
            + str(lat_mat[2][1])
            + " "
            + str(lat_mat[2][2])
        )
        return line

    def to_string(self):
        """Convert inputs to a string to write in file."""
        control = ""
        system = ""
        electrons = ""
        ions = ""
        cell = ""
        input = ""
        inputph = ""
        inputa2f = ""
        spec = ""
        if self.control_params:
            control = (
                "&control\n\n"
                + self.dictionary_to_string(self.control_params)
                + "/"
                + "\n"
            )
        if self.system_params:
            system = (
                "\n&system\n\n"
                + self.dictionary_to_string(self.system_params)
                + "/"
                + "\n"
            )
        if self.electron_params:
            electrons = (
                "\n&electrons\n\n"
                + self.dictionary_to_string(self.electron_params)
                + "/"
                + "\n"
            )
        if self.ion_params:
            ions = (
                "\n&ions\n\n"
                + self.dictionary_to_string(self.ion_params)
                + "/"
                + "\n"
            )
        if self.cell_params:
            cell = (
                "\n&cell\n\n"
                + self.dictionary_to_string(self.cell_params)
                + "/"
                + "\n"
            )
        if self.input:
            input = (
                "\n&input\n\n"
                + self.dictionary_to_string(self.input)
                + "/"
                + "\n"
            )
        if self.inputph:
            inputph = (
                "\n&inputph\n\n"
                + self.dictionary_to_string(self.inputph)
                + "/"
                + "\n"
            )
        if self.inputa2f:
            inputa2f = (
                "\n&inputa2F\n\n"
                + self.dictionary_to_string(self.inputa2f)
                + "/"
                + "\n"
            )
        if self.species:
            spec = "ATOMIC_SPECIES\n\n" + self.atomic_species_string() + "\n"
        line = (
            control
            + system
            + electrons
            + ions
            + cell
            + input
            + inputph
            + inputa2f
            + spec
            # + "ATOMIC_SPECIES\n\n"
            # + self.atomic_species_string()
            # + "\n"
            + "ATOMIC_POSITIONS crystal\n\n"
            + self.atomic_pos()
            + "\n"
            + "CELL_PARAMETERS angstrom\n\n"
            + self.atomic_cell_params()
            + "\n\n"
            + self.kpoints_to_string()
        )
        return line

    def write_file(self, filename="qe.in"):
        """Write input file."""
        line = self.to_string()
        f = open(filename, "w")
        f.write(line)
        f.close()


class GenericInputs(object):
    """Obtain set of minimal QE input parameters."""

    def __init__(self):
        """Initialize with minimum parameters that can be updated."""
        self.sample_qe_inputs = {
            "control": {
                "calculation": "'scf'",
                "restart_mode": "'from_scratch'",
                "prefix": "'QE'",
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
                "force_symmorphic": ".true.",
                "nosym": ".false.",
                "noinv": ".false.",
            },
            "electrons": {
                "diagonalization": "'david'",
                "mixing_mode": "'local-TF'",
                "mixing_beta": 0.3,
                "conv_thr": "1d-9",
            },
            "ions": {"ion_dynamics": "'bfgs'"},
            "cell": {
                "cell_dynamics": "'bfgs'",
                "cell_dofree": "'all'",
                "cell_factor": 1.0,
            },
        }

    def geometry_optimization(self):
        """Obtain QE inputs for geometry optimization."""
        input = self.sample_qe_inputs
        input["control"]["calculation"] = "'vc-relax'"
        return input

    # def simple_scf(self):
    #    """Obtain QE inputs for single SCF."""
    #    input = self.sample_qe_inputs
    #    return input

    # def electron_band_structure(self):
    #    """Obtain QE inputs for bandstructure calculations optimization."""
    #    input = self.sample_qe_inputs
    #    # TODO: increase nbands
    #    input["control_params"]["calculation"] = "'nscf'"
    #    input["control_params"]["restart_mode"] = None
    #    return input


"""
if __name__ == "__main__":
    from jarvis.core.kpoints import Kpoints3D
    from jarvis.core.atoms import Atoms
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    print(Si)
    kp = Kpoints3D().automatic_length_mesh(
        lattice_mat=Si.lattice_mat, length=20
    )
    qe = QEinfile(Si, kp)
    qe.write_file()
    kp = Kpoints3D().kpath(atoms=Si)
    qe = QEinfile(Si, kp)
    qe.write_file("qe.in2")
    sp = qe.atomic_species_string()
    sp = qe.atomic_cell_params()
    print("sp", sp)
    print(qe.input_params['system_params']['nat'])
#  /users/knc6/Software/QE/q-e/bin/pw.x -i qe.in
"""
