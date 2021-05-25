"""Modules for generating input files for QE."""

import os
import requests
import tarfile
from jarvis.core.specie import Specie


class QEinfile(object):
    """
    Gnererate quantum espresso input files.

    jarvis.core.Atoms and Kpoints are required as inputs.
    For other input parameters,
    see https://www.quantum-espresso.org/Doc/INPUT_PW.html .
    """

    def __init__(
        self, atoms, kpoints, psp_dir=None, input_params={}, url=None
    ):
        """Initialize input parameters for qunatum espresso."""
        if input_params == {}:
            input_params = GenericInputs().geometry_optimization()
        if psp_dir is None:
            psp_dir = str(os.path.join(os.path.dirname(__file__), "QE_PSPs"))
            # Download GBRV PSPs by default
            if url is None:
                url = ("http://www.physics.rutgers.edu/"
                       "gbrv/all_pbesol_UPF_v1.5.tar.gz")
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

        self.psp_dir = psp_dir

        self.atoms = atoms
        input_params["system_params"]["nat"] = self.atoms.num_atoms
        self.species = self.atoms.uniq_species
        input_params["system_params"]["ntyp"] = len(self.species)

        self.kpoints = kpoints

        self.control_params = input_params["control_params"]
        self.control_params["pseudo_dir"] = str("'") + self.psp_dir + "'"
        self.system_params = input_params["system_params"]
        self.electron_params = input_params["electron_params"]
        self.ion_params = input_params["ion_params"]
        self.cell_params = input_params["cell_params"]
        self.input_params = input_params

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
            print("Kpoint scheme not implemented except linemode, & automatic")
        return kp

    def get_psp(self, element):
        """Obtain psuedopotential for an element."""
        element = str(element).lower()
        for i in os.listdir(self.psp_dir):
            el = str(i.split("_")[0]).lower()
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

    def atomic_pos(self):
        """Obtain string for QE atomic positions."""
        line = ""
        for i, j in zip(self.atoms.elements, self.atoms.frac_coords):
            line = line + str(i) + " " + " ".join(map(str, j)) + "\n"
        return line

    def atomic_cell_params(self):
        """Obtain string for QE atomic lattice parameters."""
        line = (
            str(self.atoms.lattice_mat[0][0])
            + " "
            + str(self.atoms.lattice_mat[0][1])
            + " "
            + str(self.atoms.lattice_mat[0][2])
            + "\n"
            + str(self.atoms.lattice_mat[1][0])
            + " "
            + str(self.atoms.lattice_mat[1][1])
            + " "
            + str(self.atoms.lattice_mat[1][2])
            + "\n"
            + str(self.atoms.lattice_mat[2][0])
            + " "
            + str(self.atoms.lattice_mat[2][1])
            + " "
            + str(self.atoms.lattice_mat[2][2])
        )
        return line

    def to_string(self):
        """Convert inputs to a string to write in file."""
        control = self.dictionary_to_string(self.control_params)
        system = self.dictionary_to_string(self.system_params)
        electrons = self.dictionary_to_string(self.electron_params)
        ions = self.dictionary_to_string(self.ion_params)
        cell = self.dictionary_to_string(self.cell_params)
        line = (
            "&control\n\n"
            + control
            + "/"
            + "\n&system\n\n"
            + system
            + "/"
            + "\n&electrons\n\n"
            + electrons
            + "/"
            + "\n&ions\n\n"
            + ions
            + "/"
            + "\n&cell\n\n"
            + cell
            + "/"
            + "\n"
            + "ATOMIC_SPECIES\n\n"
            + self.atomic_species_string()
            + "\n"
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
            "control_params": {
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
            "system_params": {
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
            "electron_params": {
                "diagonalization": "'david'",
                "mixing_mode": "'local-TF'",
                "mixing_beta": 0.3,
                "conv_thr": "1d-9",
            },
            "ion_params": {"ion_dynamics": "'bfgs'"},
            "cell_params": {
                "cell_dynamics": "'bfgs'",
                "cell_dofree": "'all'",
                "cell_factor": 1.0,
            },
        }

    def geometry_optimization(self):
        """Obtain QE inputs for geometry optimization."""
        input = self.sample_qe_inputs
        input["control_params"]["calculation"] = "'vc-relax'"
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
