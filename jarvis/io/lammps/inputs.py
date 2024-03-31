"""Module for writing LAMMPS input files."""

import numpy as np
import pprint
from collections import OrderedDict
from jarvis.core.atoms import Atoms
from jarvis.core.specie import Specie


class LammpsData(object):
    """Construct Lammps data file."""

    def __init__(
        self,
        lammps_box=[],
        species=[],
        charges=[],
        cart_coords=[],
        element_order=[],
    ):
        """Provide following information."""
        self._lammps_box = lammps_box
        self._species = species
        self._charges = charges
        self._cart_coords = cart_coords
        self._element_order = element_order

    def lammps_to_atoms(self):
        """Convert Lammps data to Atoms object."""
        # n_atoms = len(self._species)
        # n_atom_types = len(self._element_order)
        box = self._lammps_box
        xhi = box[1]
        xlo = box[0]
        yhi = box[3]
        ylo = box[2]
        zhi = box[5]
        zlo = box[4]
        xy = box[6]
        xz = box[7]
        yz = box[8]
        lat = np.array(
            [[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]]
        )
        elements = [self._element_order[i - 1] for i in self._species]
        atoms = Atoms(
            lattice_mat=lat,
            elements=elements,
            coords=self._cart_coords,
            cartesian=True,
        )
        return atoms

    def atoms_to_lammps(self, atoms, origin=(0, 0, 0)):
        """Convert Atoms object to Lammps data."""
        # n_atoms = len(self._species)
        elements = atoms.elements
        # num_atoms = atoms.num_atoms
        a, b, c = atoms.lattice.abc
        xlo, ylo, zlo = origin
        xhi = a + xlo
        m = atoms.lattice._lat
        xy = np.dot(m[1], m[0] / a)
        yhi = np.sqrt(b**2 - xy**2) + ylo
        xz = np.dot(m[2], m[0] / a)
        yz = (np.dot(m[1], m[2]) - xy * xz) / (yhi - ylo)
        zhi = np.sqrt(c**2 - xz**2 - yz**2) + zlo
        rot_matrix = np.linalg.solve(
            [[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]], m
        )
        # box = np.array([[xlo, xhi], [ylo, yhi], [zlo, zhi], [xy, xz, yz]])
        box = np.array([xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz])
        new_coords = (np.dot(rot_matrix, atoms.cart_coords.T)).T
        unique_elements = list(set(elements))
        # n_atom_types = len(unique_elements)
        n_atoms = len(elements)
        if self._charges == []:
            self._charges = np.zeros(n_atoms)
        lmp_species = []
        for i, ii in enumerate(elements):
            s = unique_elements.index(ii) + 1
            lmp_species.append(s)
        lmp_species = np.array(lmp_species, dtype="int")
        return LammpsData(
            lammps_box=box,
            species=lmp_species,
            charges=self._charges,
            cart_coords=new_coords,
            element_order=unique_elements,
        )

    def read_data(
        self,
        filename="lammps.data",
        element_order=[],
        potential_file="pot.mod",
        verbose=False,
        has_charges=True,
    ):
        """Read Lammps data file."""
        # n_atoms = len(self._species)
        if element_order == []:
            # Reading potential file for element order
            if verbose:
                print("element_order is empty, reading from", potential_file)
            pot_file = open(potential_file, "r")
            lines = pot_file.read().splitlines()
            pot_file.close()
            symb = []
            # count = 0

            for i, line in enumerate(lines):
                if "pair_coeff" in line.split():
                    sp = line.split()
                    # print("spsplit", sp, os.getcwd())
                    for el in sp:
                        try:
                            if str(Specie(el).Z) != "nan":
                                # if el=='M':
                                #    el='Mo'
                                # count=count+1
                                # if count >4:
                                symb.append(Specie(el).symbol)
                        except Exception:
                            pass
        else:
            symb = [Specie(i).symbol for i in element_order]

        # print("symb=", symb)
        f = open(filename, "r")
        lines = f.read().splitlines()
        xy = 0
        xz = 0
        yz = 0
        for i, line in enumerate(lines):
            if "atoms" in line.split():
                natoms = int(line.split()[0])
            if "types" in line.split():
                # print(line)
                ntypes = int(line.split()[0])
            if "xlo" in line.split():
                xlo = float(line.split()[0])
                xhi = float(line.split()[1])
            if "ylo" in line.split():
                ylo = float(line.split()[0])
                yhi = float(line.split()[1])
            if "zlo" in line.split():
                zlo = float(line.split()[0])
                zhi = float(line.split()[1])
            if "xy" in line.split():
                xy = float(line.split()[0])
                xz = float(line.split()[1])
                yz = float(line.split()[2])
        if len(symb) != ntypes:
            ValueError(
                "Something wrong in atom type assignment", len(symb), ntypes
            )
        lat = np.array(
            [[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]]
        )
        typ = np.empty((natoms), dtype="S20")
        x = np.zeros((natoms))
        y = np.zeros((natoms))
        z = np.zeros((natoms))
        q = np.zeros((natoms))
        coords = list()
        if has_charges:
            it = 2
        else:
            it = 1
        for i, line in enumerate(lines):
            if "Atoms" in line.split():
                for j in range(0, natoms):
                    # print int(((lines[j+2]).split()[1]))-1
                    typ[j] = symb[int(((lines[i + j + 2]).split()[1])) - 1]
                    if has_charges:
                        q[j] = float((lines[i + j + 2]).split()[2])
                    x[j] = float((lines[i + j + 2]).split()[it + 1])
                    y[j] = float((lines[i + j + 2]).split()[it + 2])
                    z[j] = float((lines[i + j + 2]).split()[it + 3])
                    coords.append([x[j], y[j], z[j]])
                    # print(coords[-1])
        f.close()
        # print ("info",(typ),'coo',(coords),'latt',lat)
        typ_sp = [str(i, "utf-8") for i in typ]
        # print ('typ_sp',typ_sp)
        atoms = Atoms(
            lattice_mat=lat,
            elements=typ_sp,
            coords=np.array(coords),
            cartesian=True,
        )
        return atoms

    def write_file(self, filename="lammps.data"):
        """Write Lammps data input file."""
        # n_atoms = len(self._species)
        n_atoms = len(self._species)
        n_atom_types = len(self._element_order)
        box = self._lammps_box
        xhi = box[1]
        xlo = box[0]
        yhi = box[3]
        ylo = box[2]
        zhi = box[5]
        zlo = box[4]
        xy = box[6]
        xz = box[7]
        yz = box[8]

        f = open(filename, "w")
        f.write("datafile (written by JARVIS-Tools) \n\n")
        f.write("%d \t atoms \n" % n_atoms)
        f.write("%d  atom types\n" % n_atom_types)
        f.write("%s %s  xlo xhi\n" % (xlo, xhi))
        f.write("%s %s  ylo yhi\n" % (ylo, yhi))
        f.write("%s %s  zlo zhi\n" % (zlo, zhi))
        f.write("%s %s %s  xy xz yz\n" % (xy, xz, yz))
        f.write("\n\n")
        f.write("Atoms \n\n")
        for i in range(n_atoms):
            s = self._species[i]
            r = self._cart_coords[i]
            charge = self._charges[i]
            f.write(
                "%6d %3d %6f %s %s %s\n" % (i + 1, s, charge, r[0], r[1], r[2])
            )
        f.close()

    def to_dict(self):
        """Convert the infor to a dictionary."""
        d = OrderedDict()
        d["box"] = self._lammps_box
        d["species"] = self._species
        d["charges"] = self._charges
        d["cart_coords"] = self._cart_coords
        d["element_order"] = self._element_order
        return d

    @classmethod
    def from_dict(self, d={}):
        """Construct from a dictionary."""
        return LammpsData(
            lammps_box=d["box"],
            species=d["species"],
            charges=d["charges"],
            cart_coords=d["cart_coords"],
            element_order=d["element_order"],
        )

    def __repr__(self, indent=4):
        """Print method."""
        return pprint.pformat(self.to_dict(), indent=indent)


class LammpsInput(object):
    """Construct LAMMPS input."""

    def __init__(self, LammpsDataObj=None, pbc=["p", "p", "p"]):
        """Require LammpsData class and periodic boundary conditions."""
        self.LammpsDataObj = LammpsDataObj
        self.pbc = pbc

    def to_dict(self):
        """Convert to a dictionary."""
        d = OrderedDict()
        d["LammpsDataObj"] = self.LammpsDataObj.to_dict()
        d["pbc"] = self.pbc
        return d

    @classmethod
    def from_dict(self, d={}):
        """Costruct class from a dictionary."""
        return LammpsInput(
            LammpsDataObj=LammpsData.from_dict(d["LammpsDataObj"]),
            pbc=d["pbc"],
        )

    def write_lammps_in(
        self,
        lammps_in="init.mod",
        lammps_in1="potential.mod",
        lammps_in2="in.main",
        lammps_trj=None,
        lammps_data=None,
        parameters={},
    ):
        """
        Write lammps input file.

        From ase with custom modifications
        LAMMPS input is devided into three parts

        Args:
            lammps_in: generally"init.mod",
            with unit and conversion factor information

            lammps_in1: generally "potential.mod",
            with force-field/potential style and
            element type information

            lammps_in2: generally "in.elastic",
            a generic main input file to be fed
            in LAMMPS usin lmp_*<...,parameters['exec']

            parameters: input parameters

        """
        f = open(lammps_in, "w")
        f1 = open(lammps_in1, "w")  # potential.mod
        f2 = open(lammps_in2, "w")
        f.write(
            ('variable dump_file string "%s"\n' % lammps_trj)
            + ("variable up  equal 1.0e-6\n")
            + ("variable cfac  equal 1.0e-4\n")
            + ("variable cunits  string GPa\n")
            + ("variable etol  equal 0\n")
            + ("variable ftol  equal 1.0e-10\n")
            + ("variable maxiter equal 1000\n")
            + ("variable maxeval equal 10000\n")
            + ("variable dmax equal 1.0e-2\n")
            + ('variable data_file string "%s"\n' % "data")
        )
        if "control_file" in parameters:
            f2.write("include %s \n" % parameters["control_file"])
        if "units" in parameters:
            f.write("units %s \n" % parameters["units"])
        else:
            f.write("units metal \n")
        if "atom_style" in parameters:
            f.write("atom_style %s \n" % parameters["atom_style"])
        else:
            f.write("atom_style atomic \n")
        if "boundary" in parameters:
            f.write("boundary %s \n" % parameters["boundary"])
        else:
            pbc = self.pbc
            f.write("boundary %s %s %s \n" % (pbc[0], pbc[1], pbc[2]))
        f.write("atom_modify sort 0 0.0 \n")
        for key in ("neighbor", "newton"):
            if key in parameters:
                f.write("%s %s \n" % (key, parameters[key]))
        f.write("\n")
        f.write("read_data %s\n" % "data")
        f.write("\n### interactions \n")
        if "lib" in parameters:
            lib = parameters["lib"]
            f1.write("%s \n" % lib)
        if ("pair_style" in parameters) and ("pair_coeff" in parameters):
            pair_style = parameters["pair_style"]
            f1.write("pair_style %s \n" % pair_style)
            symbols = self.LammpsDataObj._element_order
            tag = ""
            for i in symbols:
                tag = tag + " " + i
            pair_coef = "* * " + str(parameters["pair_coeff"]) + " " + tag
            f1.write("pair_coeff %s \n" % pair_coef)

            masses = []
            for i in symbols:
                m = Specie(i).atomic_mass
                if m not in masses:
                    masses.append(m)
            count = 0
            for i in masses:
                count = count + 1
                f.write("mass" + " " + str(count) + " " + str(i) + "\n")

        else:

            # default interaction
            f.write(
                "pair_style lj/cut 2.5 \n"
                + "pair_coeff * * 1 1 \n"
                + "mass * 1.0 \n"
            )
        f1.write("neighbor 1.0 nsq\n")
        f1.write("neigh_modify once no every 1 delay 0 check yes\n")
        if "min" not in parameters:
            f1.write("min_style  cg\n")
            f1.write("min_modify           dmax ${dmax} line quadratic\n")
        f1.write("thermo          1\n")
        f1.write(
            "thermo_style custom step temp press cpu pxx pyy pzz\
             pxy pxz pyz ke pe etotal vol lx ly lz atoms\n"
        )
        f1.write("thermo_modify norm no\n")
        if "fix" in parameters:
            if parameters["fix"]:
                for i in parameters["fix"]:
                    f1.write("fix %s\n" % i)
        f.close()
        f1.close()
        f2.close()


"""
if __name__ == "__main__":
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(atoms=Si)
    Si = spg.conventional_standard_structure
    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-FF/FS/Al1.eam.fs_nist/bulk@mp-134_fold/mp-134/new_pymatgen_slab.vasp"
    )
    p = Poscar.from_file(
        "/rk2/knc6/JARVIS-FF/COMB/ffield.comb3.NiAlO_nist/bulk@mp-1143_fold/bulk@mp-1143/new_pymatgen_slab.vasp"
    )
    # print(lmp)
    lmp = LammpsData().atoms_to_lammps(atoms=p.atoms)
    lmp.write_file()
    x = lmp.lammps_to_atoms()
    # print (x)
    lmp = LammpsData().read_data("lammps.data")
    # print(lmp)
    lmp = LammpsData().atoms_to_lammps(atoms=p.atoms)
    pair_coeff = "abc"
    LammpsInput(LammpsDataObj=lmp).write_lammps_in(
        parameters={
            "pair_style": "eam/alloy",
            "pair_coeff": pair_coeff,
            "atom_style": "charge",
            "control_file": "/users/knc6/inelast.mod",
        }
    )
"""
