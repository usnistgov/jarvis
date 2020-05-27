"""
Modules handling input files for VASP calculations
"""
from collections import Counter
import pprint
import json
import os
import numpy as np
from jarvis.core.atoms import Atoms
from collections import OrderedDict
from jarvis.core.kpoints import generate_kgrid, Kpoints3D
from jarvis.core.utils import get_counts


class Poscar(object):
    """
    Class defining Poscar object 

    Args:
    
        atoms : Atoms object
        
        comment : Header of Poscar file
    """

    def __init__(self, atoms, comment="System"):
        self.atoms = atoms
        self.comment = comment

    @staticmethod
    def from_file(filename="POSCAR"):
        """
        Read simple POSCAR file from the path
        """
        with open(filename, "r") as f:
            return Poscar.from_string(f.read())

    def write_file(self, filename):
        """
        Write the Poscar object to a file
        """
        
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
        order = np.argsort(self.atoms.elements)
        coords = self.atoms.frac_coords
        coords_ordered = np.array(coords)#[order]
        elements_ordered = np.array(self.atoms.elements)#[order]
        props_ordered = np.array(self.atoms.props)#[order]
        check_selective_dynamics = False
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
            rest = rest + " ".join(map(str, i)) + " " + str(props_ordered[ii]) + "\n"

        result = header + middle + rest

        f.write(result)
        f.close()

    @staticmethod
    def from_string(lines):
        """
        Read Poscar from strings, useful in reading files/streams
        """
        
        text = lines.splitlines()
        comment = text[0]
        scale = float(text[1])
        lattice_mat = []
        lattice_mat.append([float(i) for i in text[2].split()])
        lattice_mat.append([float(i) for i in text[3].split()])
        lattice_mat.append([float(i) for i in text[4].split()])
        lattice_mat = np.array(lattice_mat)

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
        formula = atoms.composition.formula

        return Poscar(atoms, comment=comment)

    def __repr__(self):
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
        check_selective_dynamics = False
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
            rest = rest + " ".join(map(str, i)) + " " + str(props_ordered[ii]) + "\n"
        result = header + middle + rest
        return result


class Incar(object):
    """
    VASP INCAR files as python dictionary of keys and values
    """

    def __init__(self, tags={}):
        self._tags = tags

    @staticmethod
    def from_file(filename="INCAR"):
        """
        Read INCAR file
        """
        
        with open(filename, "r") as f:
            return Incar.from_string(f.read())

    def update(self, d={}):
        """
        Provide the new tags as a dictionary to update Incar object
        """
        
        #print("selftags1=", self._tags)
        if self._tags != {}:
            for i, j in self._tags.items():
                self._tags[i] = j  # .strip(' ')
        for i, j in d.items():
            self._tags[i] = j
        print("selftags2=", self._tags)
        return Incar(self._tags)

    def get(self, key="POTIM", temp=0.5):
        """
        Get the value of a key
        """
        
        if key in list(self._tags.keys()):
            return self._tags[key]
        else:
            self._tags[key] = temp
            print("Setting the temp value to the key", temp, key)
            return self._tags[key]

    def to_dict(self):
        """
        Convert into dictionary
        """
        
        return self._tags

    def from_dict(self, data={}):
        return Incar(tags=data)

    @staticmethod
    def from_string(lines):
        text = lines.splitlines()
        tags = OrderedDict()
        for i in text:
            if "=" in i:
                tmp = i.split("=")
                tags.setdefault(tmp[0].strip(" "), tmp[1].strip(" "))

        return Incar(tags=tags)

    def __repr__(self):
        return str(self._tags)

    def write_file(self, filename):
        """
        Write Incar to a file
        """
        
        tags = self._tags
        lines = ""
        for i, j in tags.items():
            lines = lines + str(i) + str("=") + str(j) + "\n"
        f = open(filename, "w")
        f.write(lines)
        f.close()


class IndividualPotcarData(object):
    """
    Class for individual POTCAR file handling
    """

    def __init__(self, data={}):
        self._data = data

    def from_string(lines):
        """
        Some of the contents in the POTCAR
        """
        
        text = lines.splitlines()
        individual_data = OrderedDict()
        individual_data["header1"] = text[0]
        individual_data["header2"] = text[2]
        individual_data["VRHFIN"] = text[3].split("=")[1]
        individual_data["element"] = text[3].split("=")[1].split(":")[0]
        individual_data["LEXCH"] = text[4].split("=")[1].replace(" ", "")
        individual_data["EATOM"] = text[5].split("=")[1].split()[0].replace(" ", "")
        individual_data["TITEL"] = text[7].split("=")[1].replace(" ", "")

        return IndividualPotcarData(data=individual_data)

    def from_file(filename="POTCAR"):
        with open(filename, "r") as f:
            return IndividualPotcarData.from_string(f.read())

    def __repr__(self):
        return pprint.pformat(self._data, indent=4)


class Potcar(object):
    def __init__(
        self,
        elements=[],
        pot_type="POT_GGA_PAW_PBE",
        pot_json_path="",
        potcar_strings={},
    ):
        self._elements = elements
        self._pot_type = pot_type
        self._potcar_strings = potcar_strings
        self._pot_json_path = pot_json_path

        if self._potcar_strings == {}:
            pot_json_file = str(
                os.path.join(os.path.dirname(__file__), "default_potcars.json")
            )
            pot_json = open(pot_json_file, "r")
            pot_json_selected = json.load(pot_json)
            pot_json.close()
            potcar_strings = OrderedDict()
            for i in self._elements:
                for j, k in pot_json_selected.items():
                    if i == j:
                        potcar_strings.setdefault(i, k)
            self._potcar_strings = potcar_strings
            if len(self._elements) != len(self._potcar_strings.keys()):
                raise ValueError(
                    "Number of elements is not same as potcar_strings",
                    self._elements,
                    self._potcar_strings.keys(),
                )
        else:
            pot_json = open(self._pot_json_path, "r")
            pot_json_selected = json.load(pot_json)
            pot_json.close()
            potcar_strings = OrderedDict()
            for i, j in pot_json_selected.items():
                if i in self._elements:
                    potcar_strings.setdefault(i, j)
            self._potcar_strings = potcar_strings
            if len(self._elements) != len(self._potcar_strings.keys()):
                raise ValueError("Number of elements is not same as potcar_strings")

    def catenate_potcar_files(self, destination_filename="POTCAR", filenames=[]):
        with open(destination_filename, "w") as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def list_potcar_files(self):
        pot_files = []
        vasp_dir = str(os.environ["JARVIS_VASP_PSP_DIR"])
        vasp_psp_dir = str(os.path.join(vasp_dir, self._pot_type))
        potcar_strings = self._potcar_strings
        for i, j in potcar_strings.items():
            tmp = os.path.join(vasp_psp_dir, j, "POTCAR")
            pot_files.append(tmp)
        return pot_files

    def write_file(self, filename="POTCAR"):
        pot_files = self.list_potcar_files()
        self.catenate_potcar_files(destination_filename=filename, filenames=pot_files)

    def __repr__(self):
        return str(str(self._pot_type) + "\n" + str(self._potcar_strings))


class Kpoints(object):
    """
    VASP KPOINTS as object
    """

    def __init__(self, filename=""):
        self.filename = filename
        if filename != "":
            f = open(self.filename, "r")
            self.lines = f.read().splitlines()
            f.close()
            self.kpoints = self.read(self.lines)

    @classmethod
    def read(self, lines):

        if lines[1] == "0":
            return self.get_mesh_kp(lines=lines)
        elif "Reciprocal" in lines[2]:
            return self.get_ibz_kp(lines=lines)
        else:
            ValueError("K-point scheme is not implemented")

    def get_mesh_kp(lines=""):
        #print("lines", lines)
        """
        Read Kpoints as grid
        """
        
        grid = [int(i) for i in lines[3].split()]
        # print (grid)
        kpts = generate_kgrid(grid)
        kpoints = Kpoints3D(kpoints=np.array(kpts))
        return kpoints

    def get_ibz_kp(lines=""):
        """
        Read the Kpoints in the line-mode
        """
        
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
        kpoints = Kpoints3D(kpoints=all_kp, labels=labels, kpoint_mode="linemode")
        return kpoints
        # print ('kp_labels',labels,len(all_kp))
        # print ('kp2',kp_labels_points, kp_labels,np.array(all_kp,dtype='float'))

        # return kp_labels_points, kp_labels,np.array(all_kp,dtype='float')


"""
if __name__ == "__main__":

    kp = open(
        "../../examples/vasp/SiOptb88/MAIN-RELAX-bulk@mp_149/KPOINTS", "r"
    )  # .read_file()
    lines = kp.read().splitlines()
    kp.close()
    print('lbl',Kpoints.read(lines).labels)
    # print (kp.kpoints)
    import sys

    sys.exit()
    kp = Kpoints(filename="../../examples/vasp/SiOptb88/MAIN-RELAX-bulk@mp_149/KPOINTS")

    p = Poscar.from_file(
        filename="/rk2/knc6/JARVIS-DFT/2D-1L/POSCAR-mp-2815-1L.vasp_PBEBO/MAIN-RELAX-Surf-mp-2815/POSCAR"
    )
    p.write_file("POSCAR")
    from jarvis.db.jsonutils import dumpjson

    dumpjson(p.atoms.to_dict(), "pp.json")
    # print (p.atoms*[3,3,3])
    # p = Incar.from_file(filename='/rk2/knc6/JARVIS-DFT/2D-1L/POSCAR-mp-2815-1L.vasp_PBEBO/MAIN-RELAX-Surf-mp-2815/INCAR')
    # print (p)
    # p.write_file('pp')
    # elements=['S','P']
    # p=Potcar(elements=elements)
    # print (p)
    # p.write_file('POTCAR')
    # inp=IndividualPotcarData.from_file('POTCAR')
    # print (inp)
"""
