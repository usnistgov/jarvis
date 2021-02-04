"""Modules handling input files for VASP calculations."""
import pprint
import json
import os
import numpy as np
from jarvis.core.atoms import Atoms
from collections import OrderedDict
from jarvis.core.kpoints import generate_kgrid, Kpoints3D
from jarvis.core.utils import get_counts
from jarvis.core.specie import Specie
from jarvis.core.utils import update_dict


class Poscar(object):
    """
    Class defining Poscar object.

    Constructued from the Atoms object.

    Args:
        atoms : Atoms object

        comment : Header of Poscar file
    """

    def __init__(self, atoms, comment="System"):
        """Initialize the Poscar object."""
        self.atoms = atoms
        self.comment = comment

    @staticmethod
    def from_file(filename="POSCAR"):
        """Read simple POSCAR file from the path."""
        with open(filename, "r") as f:
            return Poscar.from_string(f.read())

    def to_dict(self):
        """Convert Poscar object to a dictionary."""
        d = OrderedDict()
        d["atoms"] = self.atoms.to_dict()
        d["comment"] = self.comment
        return d

    @classmethod
    def from_dict(self, d={}):
        """Construct Poscar object from a dictionary."""
        return Poscar(atoms=Atoms.from_dict(d["atoms"]), comment=d["comment"])

    def write_file(self, filename):
        """Write the Poscar object to a file."""
        f = open(filename, "w")
        header = (
            str(self.comment)
            + str("\n1.0\n")
            + str(self.atoms.lattice_mat[0][0])
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
            + "\n"
        )
        # order = np.argsort(self.atoms.elements)
        coords = self.atoms.frac_coords
        # DO NOT USE ORDER
        coords_ordered = np.array(coords)  # [order]
        elements_ordered = np.array(self.atoms.elements)  # [order]
        props_ordered = np.array(self.atoms.props)  # [order]
        # check_selective_dynamics = False
        counts = get_counts(elements_ordered)
        if "T" in "".join(map(str, self.atoms.props[0])):
            middle = (
                " ".join(map(str, counts.keys()))
                + "\n"
                + " ".join(map(str, counts.values()))
                + "\nSelective dynamics\n"
                + "Direct\n"
            )
        else:
            middle = (
                " ".join(map(str, counts.keys()))
                + "\n"
                + " ".join(map(str, counts.values()))
                + "\ndirect\n"
            )
        rest = ""
        # print ('repr',self.frac_coords, self.frac_coords.shape)
        for ii, i in enumerate(coords_ordered):
            p_ordered = str(props_ordered[ii])
            rest = rest + " ".join(map(str, i)) + " " + str(p_ordered) + "\n"

        result = header + middle + rest

        f.write(result)
        f.close()

    @staticmethod
    def from_string(lines):
        """Read Poscar from strings, useful in reading files/streams."""
        text = lines.splitlines()
        comment = text[0]
        scale = float(text[1])
        lattice_mat = []
        lattice_mat.append([float(i) for i in text[2].split()])
        lattice_mat.append([float(i) for i in text[3].split()])
        lattice_mat.append([float(i) for i in text[4].split()])
        lattice_mat = scale * np.array(lattice_mat)

        uniq_elements = text[5].split()
        element_count = np.array([int(i) for i in text[6].split()])
        elements = []
        for i, ii in enumerate(element_count):
            for j in range(ii):
                elements.append(uniq_elements[i])
        cartesian = True
        if "d" in text[7] or "D" in text[7]:
            cartesian = False
        # print ('cartesian poscar=',cartesian,text[7])
        num_atoms = int(np.sum(element_count))
        coords = []
        for i in range(num_atoms):
            coords.append([float(i) for i in text[8 + i].split()[0:3]])
        coords = np.array(coords)
        atoms = Atoms(
            lattice_mat=lattice_mat,
            coords=coords,
            elements=elements,
            cartesian=cartesian,
        )
        # print (atoms)
        # formula = atoms.composition.formula

        return Poscar(atoms, comment=comment)

    def __repr__(self):
        """Represent the Poscar class."""
        header = (
            str(self.comment)
            + str("\n1.0\n")
            + str(self.atoms.lattice_mat[0][0])
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
            + "\n"
        )
        order = np.argsort(self.atoms.elements)
        coords = self.atoms.frac_coords
        coords_ordered = np.array(coords)[order]
        elements_ordered = np.array(self.atoms.elements)[order]
        props_ordered = np.array(self.atoms.props)[order]
        # check_selective_dynamics = False
        counts = get_counts(elements_ordered)
        if "T" in "".join(map(str, self.atoms.props[0])):
            middle = (
                " ".join(map(str, counts.keys()))
                + "\n"
                + " ".join(map(str, counts.values()))
                + "\nSelective dynamics\n"
                + "Direct\n"
            )
        else:
            middle = (
                " ".join(map(str, counts.keys()))
                + "\n"
                + " ".join(map(str, counts.values()))
                + "\ndirect\n"
            )
        rest = ""
        # print ('repr',self.frac_coords, self.frac_coords.shape)
        for ii, i in enumerate(coords_ordered):
            p_ordered = str(props_ordered[ii])
            rest = rest + " ".join(map(str, i)) + " " + str(p_ordered) + "\n"
        result = header + middle + rest
        return result


class Incar(object):
    """Make VASP INCAR files as python dictionary of keys and values."""

    def __init__(self, tags={}):
        """Initialize with a dictionary."""
        self._tags = tags

    @staticmethod
    def from_file(filename="INCAR"):
        """Read INCAR file."""
        with open(filename, "r") as f:
            return Incar.from_string(f.read())

    def update(self, d={}):
        """Provide the new tags as a dictionary to update Incar object."""
        # print("selftags1=", self._tags)
        if self._tags != {}:
            for i, j in self._tags.items():
                self._tags[i] = j  # .strip(' ')
        for i, j in d.items():
            self._tags[i] = j
        # print("selftags2=", self._tags)
        return Incar(self._tags)

    # def get(self, key="POTIM", temp=0.5):
    #    """Get the value of a key."""
    #    if key in list(self._tags.keys()):
    #        return self._tags[key]
    #    else:
    #        self._tags[key] = temp
    #        print("Setting the temp value to the key", temp, key)
    #        return self._tags[key]

    def to_dict(self):
        """Convert into dictionary."""
        return self._tags

    @classmethod
    def from_dict(self, data={}):
        """Construct from dictionary."""
        return Incar(tags=data)

    @staticmethod
    def from_string(lines):
        """Construct from string."""
        text = lines.splitlines()
        tags = OrderedDict()
        for i in text:
            if "=" in i:
                tmp = i.split("=")
                tags.setdefault(tmp[0].strip(" "), tmp[1].strip(" "))

        return Incar(tags=tags)

    def __repr__(self):
        """Representation of the class."""
        return str(self._tags)

    def write_file(self, filename):
        """Write Incar to a file."""
        tags = self._tags
        lines = ""
        for i, j in tags.items():
            lines = lines + str(i) + str("=") + str(j) + "\n"
        f = open(filename, "w")
        f.write(lines)
        f.close()


class IndividualPotcarData(object):
    """Class for individual POTCAR file handling."""

    def __init__(self, data={}):
        """Intialize with some key info."""
        self._data = data

    def from_string(lines):
        """Generate some of the contents in the POTCAR."""
        text = lines.splitlines()
        individual_data = OrderedDict()
        individual_data["header1"] = text[0]
        individual_data["header2"] = text[2]
        individual_data["VRHFIN"] = text[3].split("=")[1]
        individual_data["element"] = text[3].split("=")[1].split(":")[0]
        individual_data["LEXCH"] = text[4].split("=")[1].replace(" ", "")
        tmp = text[5].split("=")[1].split()[0].replace(" ", "")
        individual_data["EATOM"] = tmp
        individual_data["TITEL"] = text[7].split("=")[1].replace(" ", "")
        return IndividualPotcarData(data=individual_data)

    def from_file(filename="POTCAR"):
        """Read from file."""
        with open(filename, "r") as f:
            return IndividualPotcarData.from_string(f.read())

    def __repr__(self):
        """Reprent the class."""
        return pprint.pformat(self._data, indent=4)


class Potcar(object):
    """Construct muti-atoms Postcar."""

    def __init__(
        self,
        elements=[],
        pot_type="POT_GGA_PAW_PBE",
        pot_json_path="",
        potcar_strings=[],
    ):
        """
        Initialize the Potcar class.

        POTCAR file contains the pseudopotential for each
        atomic species used in the calculation.
        Args:
            elements: atomic elements.

            pot_type: Type of pseudopotential.
            Look for VASP provided PSPs.
            VASP_PSP_DIR should be in the path.

            pot_json_path: Path to dictionary of chemical elements along with
            their orbitals taken into conisideration, e.g. V_pv.

            potcar_strings: One can directly provide
            the above mentioned strings.
        """
        self._elements = elements
        self._pot_type = pot_type
        self._potcar_strings = potcar_strings
        self._pot_json_path = pot_json_path

        if self._potcar_strings == []:
            pot_json_file = str(
                os.path.join(os.path.dirname(__file__), "default_potcars.json")
            )
            self._pot_json_path = pot_json_file
            pot_json = open(pot_json_file, "r")
            pot_json_selected = json.load(pot_json)
            pot_json.close()
            potcar_strings = []  # OrderedDict()
            for i in self._elements:
                potcar_strings.append(pot_json_selected[i])
                # for j, k in pot_json_selected.items():
                #    if i == j:
                #        potcar_strings.setdefault(i, k)
            self._potcar_strings = potcar_strings
            #  print("self._elements", self._elements)
            #  print("self._potcar_strings", self._potcar_strings)
            if len(self._elements) != len(self._potcar_strings):
                raise ValueError(
                    "Number of elements is not same as potcar_strings",
                    self._elements,
                    self._potcar_strings.keys(),
                )
        else:
            pot_json = open(self._pot_json_path, "r")
            pot_json_selected = json.load(pot_json)
            pot_json.close()
            self._potcar_strings = potcar_strings
            if len(self._elements) != len(self._potcar_strings):
                msg = "Number of elements not same as potcar_strings"
                raise ValueError(msg)

    @classmethod
    def from_dict(self, d={}):
        """Build class from a dictionary."""
        return Potcar(
            elements=d["elements"],
            pot_type=d["pot_type"],
            pot_json_path=d["pot_json_path"],
            potcar_strings=d["potcar_strings"],
        )

    def to_dict(self):
        """Convert to a dictionary."""
        d = OrderedDict()
        d["elements"] = np.array(self._elements).tolist()
        d["pot_type"] = self._pot_type
        d["pot_json_path"] = self._pot_json_path
        d["potcar_strings"] = np.array(self._potcar_strings).tolist()
        return d

    def catenate_potcar_files(
        self, destination_filename="POTCAR", filenames=[]
    ):
        """Catenate potcars of sifferent elements."""
        with open(destination_filename, "w") as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def list_potcar_files(self):
        """List POTCAR files."""
        pot_files = []
        vasp_dir = str(os.environ["VASP_PSP_DIR"])
        vasp_psp_dir = str(os.path.join(vasp_dir, self._pot_type))
        potcar_strings = self._potcar_strings
        for j in potcar_strings:
            tmp = os.path.join(vasp_psp_dir, j, "POTCAR")
            pot_files.append(tmp)
        return pot_files

    def write_file(self, filename="POTCAR"):
        """Write POTCAR file."""
        pot_files = self.list_potcar_files()
        self.catenate_potcar_files(
            destination_filename=filename, filenames=pot_files
        )

    def __repr__(self):
        """Represent Potcar."""
        return str(str(self._pot_type) + "\n" + str(self._potcar_strings))


class Kpoints(object):
    """Make VASP KPOINTS as object."""

    def __init__(self, filename=""):
        """Initialize Kpoints from filename else read from file-stream."""
        self.filename = filename
        if filename != "":
            f = open(self.filename, "r")
            self.lines = f.read().splitlines()
            f.close()
            self.kpoints = self.read(self.lines)

    @classmethod
    def read(self, lines):
        """Read from an open file."""
        if lines[1] == "0":
            return self.get_mesh_kp(lines=lines)
        elif "Reciprocal" in lines[2]:
            return self.get_ibz_kp(lines=lines)
        else:
            ValueError("K-point scheme is not implemented")

    def get_mesh_kp(lines=""):
        """Read Kpoints as grid."""
        grid = [int(i) for i in lines[3].split()]
        kpts = generate_kgrid(grid)
        kpts_cls = Kpoints3D(kpoints=np.array(kpts))
        return kpts_cls

    def get_ibz_kp(lines=""):
        """Read the Kpoints in the line-mode."""
        kp_labels = []
        all_kp = []
        kp_labels_points = []
        for ii, i in enumerate(lines):
            if ii > 2:
                tmp = i.split()
                k = [tmp[0], tmp[1], tmp[2]]
                all_kp.append(k)
                if len(tmp) == 5:
                    tmp = str("$") + str(tmp[4]) + str("$")
                    if len(kp_labels) == 0:
                        kp_labels.append(tmp)
                        kp_labels_points.append(ii - 3)
                    elif tmp != kp_labels[-1]:
                        kp_labels.append(tmp)
                        kp_labels_points.append(ii - 3)
        labels = []
        for i in range(len(all_kp)):
            labels.append("")
        for i, j in zip(kp_labels, kp_labels_points):
            labels[j] = i
        all_kp = np.array(all_kp, dtype="float")
        kpts_cls = Kpoints3D(
            kpoints=all_kp, labels=labels, kpoint_mode="linemode"
        )
        return kpts_cls


def find_ldau_magmom(
    atoms="",
    default_magmom=True,
    U=3.0,
    mag=5.0,
    amix=0.2,
    bmix=0.00001,
    amixmag=0.8,
    bmixmag=0.00001,
    lsorbit=False,
):
    """Get necessary INCAR tags for DFT+U calculations."""
    sps = atoms.uniq_species
    LDAUL = []
    LDAUU = []
    LDAUTYPE = 2
    lmix = 4
    for i in sps:
        el = Specie(i)
        el_u = 0
        el_l = -1
        if el.element_property("is_transition_metal"):
            el_u = U
            el_l = 2
        if el.element_property("is_actinoid") or el.element_property(
            "is_lanthanoid"
        ):
            el_u = U
            el_l = 3
            lmix = 6
        LDAUL.append(el_l)
        LDAUU.append(el_u)
    if 3 in LDAUL:
        LDAUTYPE = 3
    nat = atoms.num_atoms
    magmom = str(nat) + str("*") + str(mag)
    if lsorbit:
        magmom = ""
        tmp = " 0 0 " + str(mag)
        for i in range(0, nat):
            magmom = magmom + tmp
    info = {}
    info["LDAU"] = ".TRUE."
    info["LDAUTYPE"] = LDAUTYPE
    info["LDAUL"] = " ".join(str(m) for m in LDAUL)
    info["LDAUU"] = " ".join(str(m) for m in LDAUU)
    info["LDAUPRINT"] = 2
    info["LMAMIX"] = lmix
    if not default_magmom:
        info["MAGMOM"] = magmom
    info["AMIX"] = amix
    info["BMIX"] = bmix
    info["AMIX_MAG"] = amixmag
    info["BMIX_MAG"] = bmixmag
    return info


def add_ldau_incar(
    use_incar_dict={}, atoms=None, Uval=2, lsorbit=False, default_magmom=True
):
    """Add LDAU in incase, especially made for spillage calcs."""
    info_ldau = find_ldau_magmom(
        U=Uval, atoms=atoms, lsorbit=lsorbit, default_magmom=default_magmom
    )
    tmp = update_dict(use_incar_dict, info_ldau)
    use_incar_dict = tmp
    return use_incar_dict
