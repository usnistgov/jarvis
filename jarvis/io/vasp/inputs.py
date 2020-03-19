from collections import Counter
import pprint
import json
import os
import numpy as np
from jarvis.core.atoms import Atoms
from collections import OrderedDict


class Poscar(object):
    def __init__(self, atoms, comment="System"):
        self.atoms = atoms
        self.comment = comment

    @staticmethod
    def from_file(filename="POSCAR"):
        with open(filename, "r") as f:
            return Poscar.from_string(f.read())

    def write_file(self, filename):
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
        coords_ordered = np.array(coords)[order]
        elements_ordered = np.array(self.atoms.elements)[order]
        props_ordered = np.array(self.atoms.props)[order]
        check_selective_dynamics = False
        if "T" in "".join(map(str, self.atoms.props[0])):
            middle = (
                " ".join(map(str, Counter(elements_ordered).keys()))
                + "\n"
                + " ".join(map(str, Counter(elements_ordered).values()))
                + "\nSelective dynamics\n"
                + "Direct\n"
            )
        else:
            middle = (
                " ".join(map(str, Counter(elements_ordered).keys()))
                + "\n"
                + " ".join(map(str, Counter(elements_ordered).values()))
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
        if "T" in "".join(map(str, self.atoms.props[0])):
            middle = (
                " ".join(map(str, Counter(elements_ordered).keys()))
                + "\n"
                + " ".join(map(str, Counter(elements_ordered).values()))
                + "\nSelective dynamics\n"
                + "Direct\n"
            )
        else:
            middle = (
                " ".join(map(str, Counter(elements_ordered).keys()))
                + "\n"
                + " ".join(map(str, Counter(elements_ordered).values()))
                + "\ndirect\n"
            )
        rest = ""
        # print ('repr',self.frac_coords, self.frac_coords.shape)
        for ii, i in enumerate(coords_ordered):
            rest = rest + " ".join(map(str, i)) + " " + str(props_ordered[ii]) + "\n"
        result = header + middle + rest
        return result


class Incar(object):
    def __init__(self, tags={}):
        self._tags = tags

    @staticmethod
    def from_file(filename="INCAR"):
        with open(filename, "r") as f:
            return Incar.from_string(f.read())

    def update(self, d={}):
        print("selftags1=", self._tags)
        if self._tags != {}:
            for i, j in self._tags.items():
                self._tags[i] = j  # .strip(' ')
        for i, j in d.items():
            self._tags[i] = j
        print("selftags2=", self._tags)
        return Incar(self._tags)

    def get(self, key="POTIM", temp=0.5):
        if key in list(self._tags.keys()):
            return self._tags[key]
        else:
            self._tags[key] = temp
            print("Setting the temp value to the key", temp, key)
            return self._tags[key]

    def to_dict(self):
        return self._tags

    def from_dict(self, data={}):
        print("data", data)
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
        tags = self._tags
        lines = ""
        for i, j in tags.items():
            lines = lines + str(i) + str("=") + str(j) + "\n"
        f = open(filename, "w")
        f.write(lines)
        f.close()


class IndividualPotcarData(object):
    def __init__(self, data={}):
        self._data = data

    def from_string(lines):
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

"""
if __name__ == "__main__":
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
