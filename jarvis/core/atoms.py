"""This module provides classes to specify atomic structure."""
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie, atomic_numbers_to_symbols
from jarvis.core.lattice import Lattice, lattice_coords_transformer
from collections import OrderedDict
from jarvis.core.utils import get_counts
import itertools
from jarvis.core.utils import get_angle
from jarvis.core.utils import (
    check_duplicate_coords,
    get_new_coord_for_xyz_sym,
)
import os
import math
import tempfile

amu_gm = 1.66054e-24
ang_cm = 1e-8


class Atoms(object):
    """
    Generate Atoms python object.

    >>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    >>> coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    >>> elements = ["Si", "Si"]
    >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    >>> print(round(Si.volume,2))
    40.03
    >>> Si.composition
    {'Si': 2}
    >>> round(Si.density,2)
    2.33
    >>> round(Si.packing_fraction,2)
    0.28
    >>> Si.atomic_numbers
    [14, 14]
    >>> Si.num_atoms
    2
    >>> Si.frac_coords[0][0]
    0
    >>> Si.cart_coords[0][0]
    0.0
    >>> coords = [[0, 0, 0], [1.3575 , 1.22175, 1.22175]]
    >>> round(Si.density,2)
    2.33
    >>> Si.spacegroup()
    'C2/m (12)'
    >>> Si.pymatgen_converter()!={}
    True
    """

    def __init__(
        self,
        lattice_mat=None,
        coords=None,
        elements=None,
        props=None,
        cartesian=False,
        show_props=False,
    ):
        """
        Create atomic structure.

        Requires lattice, coordinates, atom type  information.
        """
        self.lattice_mat = np.array(lattice_mat)
        self.show_props = show_props
        self.lattice = Lattice(lattice_mat)
        self.coords = np.array(coords)
        self.elements = elements
        self.cartesian = cartesian
        self.props = props
        if self.props is None:
            self.props = ["" for i in range(len(self.elements))]
        if self.cartesian:
            self.cart_coords = self.coords
            self.frac_coords = np.array(self.lattice.frac_coords(self.coords))
        else:
            self.frac_coords = self.coords
            self.cart_coords = np.array(self.lattice.cart_coords(self.coords))

    def write_cif(
        self, filename="atoms.cif", comment=None, with_spg_info=True
    ):
        """
        Write CIF format file from Atoms object.

        Caution: can't handle fractional occupancies right now
        """
        if comment is None:
            comment = "CIF file written using JARVIS-Tools package."
        comment = comment + "\n"
        f = open(filename, "w")
        f.write(comment)
        composition = self.composition
        line = "data_" + str(composition.reduced_formula) + "\n"
        f.write(line)
        from jarvis.analysis.structure.spacegroup import Spacegroup3D

        if with_spg_info:
            spg = Spacegroup3D(self)
            line = (
                "_symmetry_space_group_name_H-M  "
                + str("'")
                + str(spg.space_group_symbol)
                + str("'")
                + "\n"
            )
        else:
            line = "_symmetry_space_group_name_H-M  " + str("'P 1'") + "\n"

        f.write(line)
        a, b, c, alpha, beta, gamma = self.lattice.parameters
        f.write("_cell_length_a       %g\n" % a)
        f.write("_cell_length_b       %g\n" % b)
        f.write("_cell_length_c       %g\n" % c)
        f.write("_cell_angle_alpha    %g\n" % alpha)
        f.write("_cell_angle_beta     %g\n" % beta)
        f.write("_cell_angle_gamma    %g\n" % gamma)
        f.write("\n")
        if with_spg_info:
            line = (
                "_symmetry_Int_Tables_number  "
                + str(spg.space_group_number)
                + "\n"
            )
        else:
            line = "_symmetry_Int_Tables_number  " + str(1) + "\n"
        f.write(line)

        line = (
            "_chemical_formula_structural  "
            + str(composition.reduced_formula)
            + "\n"
        )
        f.write(line)
        line = "_chemical_formula_sum  " + str(composition.formula) + "\n"
        f.write(line)
        line = "_cell_volume  " + str(self.volume) + "\n"
        f.write(line)
        reduced, repeat = composition.reduce()
        line = "_cell_formula_units_Z  " + str(repeat) + "\n"
        f.write(line)
        f.write("loop_\n")
        f.write("  _symmetry_equiv_pos_site_id\n")
        f.write(" _symmetry_equiv_pos_as_xyz\n")
        f.write(" 1  'x, y, z'\n")
        f.write("loop_\n")
        f.write(" _atom_site_type_symbol\n")
        f.write(" _atom_site_label\n")
        f.write(" _atom_site_symmetry_multiplicity\n")
        f.write(" _atom_site_fract_x\n")
        f.write(" _atom_site_fract_y\n")
        f.write(" _atom_site_fract_z\n")
        f.write(" _atom_site_fract_occupancy\n")
        order = np.argsort(self.elements)
        coords_ordered = np.array(self.frac_coords)[order]
        elements_ordered = np.array(self.elements)[order]
        occ = 1
        extra = 1
        element_types = []
        # count = 0
        for ii, i in enumerate(elements_ordered):
            if i not in element_types:
                element_types.append(i)
                count = 0
            symbol = i
            count = count + 1
            label = str(i) + str(count)
            element_types.append(i)
            coords = coords_ordered[ii]
            f.write(
                " %s  %s  %s  %7.5f  %7.5f  %7.5f  %s\n"
                % (symbol, label, occ, coords[0], coords[1], coords[2], extra)
            )
        f.close()

    @staticmethod
    def read_with_cif2cell(filename="1000000.cif", get_primitive_atoms=False):
        """Use cif2cell package to read cif files."""
        # https://pypi.org/project/cif2cell/
        # tested on version 2.0.0a3
        try:
            new_file, fname = tempfile.mkstemp(text=True)
            cmd = (
                "cif2cell "
                + filename
                + " -p vasp --vasp-cartesian-positions --vca -o "
                + fname
            )
            os.system(cmd)

        except Exception as exp:
            print(exp)
            pass
        f = open(fname, "r")
        text = f.read().splitlines()
        f.close()
        os.close(new_file)
        elements = text[0].split("Species order:")[1].split()
        scale = float(text[1])
        lattice_mat = []
        lattice_mat.append([float(i) for i in text[2].split()])
        lattice_mat.append([float(i) for i in text[3].split()])
        lattice_mat.append([float(i) for i in text[4].split()])
        lattice_mat = scale * np.array(lattice_mat)
        uniq_elements = elements
        element_count = np.array([int(i) for i in text[5].split()])
        elements = []
        for i, ii in enumerate(element_count):
            for j in range(ii):
                elements.append(uniq_elements[i])
        cartesian = True
        if "d" in text[6] or "D" in text[6]:
            cartesian = False
        # print ('cartesian poscar=',cartesian,text[7])
        num_atoms = int(np.sum(element_count))
        coords = []
        for i in range(num_atoms):
            coords.append([float(i) for i in text[7 + i].split()[0:3]])
        coords = np.array(coords)
        atoms = Atoms(
            lattice_mat=lattice_mat,
            coords=coords,
            elements=elements,
            cartesian=cartesian,
        )
        if get_primitive_atoms:
            atoms = atoms.get_primitive_atoms
        return atoms

    @staticmethod
    def from_pdb_old(filename="abc.pdb"):
        """Read PDB file, kept of checking, use from_pdb instead."""
        f = open(filename, "r")
        lines = f.read().splitlines()
        f.close()
        coords = []
        species = []
        for i in lines:
            tmp = i.split()
            if "ATOM " in i and "REMARK" not in i and "SITE" not in i:
                coord = [float(tmp[6]), float(tmp[7]), float(tmp[8])]
                coords.append(coord)
                species.append(tmp[11])
                # print (coord,tmp[11])
        coords = np.array(coords)

        max_c = np.max(coords, axis=0)
        min_c = np.min(coords, axis=0)
        box = np.zeros((3, 3))
        for j in range(3):
            box[j, j] = abs(max_c[j] - min_c[j])
        pdb = Atoms(
            lattice_mat=box, elements=species, coords=coords, cartesian=True
        )
        mol = VacuumPadding(pdb, vacuum=20.0).get_effective_molecule()
        return mol

    @staticmethod
    def from_pdb(filename="abc.pdb", max_lat=200):
        """Read pdb/sdf/mol2 etc. file and make Atoms object, using pytraj."""
        try:
            import pytraj as pt
        except Exception:
            print("Pytraj not installed, trying native version.")
            return Atoms.from_pdb_old(filename)
        x = pt.load(filename)
        coords = x.xyz
        at_numbs = [int(i.atomic_number) for i in list(x.top.atoms)]
        elements = atomic_numbers_to_symbols(at_numbs)
        lattice_mat = [[max_lat, 0, 0], [0, max_lat, 0], [0, 0, max_lat]]
        a = Atoms(
            lattice_mat=lattice_mat,
            elements=elements,
            coords=coords[0],
            cartesian=True,
        )
        return a

    @staticmethod
    def from_cif(
        filename="atoms.cif",
        from_string="",
        get_primitive_atoms=True,
        use_cif2cell=True,
    ):
        """Read .cif format file."""
        # Warnings:
        # May not work for:
        # system with partial occupancy
        # cif file with multiple blocks
        # _atom_site_U_iso, instead of fractn_x, cartn_x
        # with non-zero _atom_site_attached_hydrogens
        try:

            if use_cif2cell:
                # https://pypi.org/project/cif2cell/
                # tested on version 2.0.0a3
                if from_string != "":
                    new_file, filename = tempfile.mkstemp(text=True)
                    f = open(filename, "w")
                    f.write(from_string)
                    f.close()
                atoms = Atoms.read_with_cif2cell(
                    filename=filename, get_primitive_atoms=get_primitive_atoms
                )

                return atoms
        except Exception as exp:
            print(exp)
            pass

        if from_string == "":
            f = open(filename, "r")
            lines = f.read().splitlines()
            # lines = [ii.encode('utf-8') for ii in f.read().splitlines()]
            f.close()
        else:
            lines = from_string.splitlines()
        lat_a = ""
        lat_b = ""
        lat_c = ""
        lat_alpha = ""
        lat_beta = ""
        lat_gamma = ""
        # TODO: check if chemical_formula_sum
        # matches Atoms.compsotion.reduced_formula

        # chemical_formula_structural = ""
        # chemical_formula_sum = ""
        # chemical_name_mineral = ""
        sym_xyz_line = ""
        for ii, i in enumerate(lines):
            if "_cell_length_a" in i:
                lat_a = float(i.split()[1].split("(")[0])
            if "_cell_length_b" in i:
                lat_b = float(i.split()[1].split("(")[0])
            if "_cell_length_c" in i:
                lat_c = float(i.split()[1].split("(")[0])
            if "_cell_angle_alpha" in i:
                lat_alpha = float(i.split()[1].split("(")[0])
            if "_cell_angle_beta" in i:
                lat_beta = float(i.split()[1].split("(")[0])
            if "_cell_angle_gamma" in i:
                lat_gamma = float(i.split()[1].split("(")[0])
            # if "_chemical_formula_structural" in i:
            #    chemical_formula_structural = i.split()[1]
            # if "_chemical_formula_sum" in i:
            #    chemical_formula_sum = i.split()[1]
            # if "_chemical_name_mineral" in i:
            #    chemical_name_mineral = i.split()[1]
            if "_symmetry_equiv_pos_as_xyz" in i:
                sym_xyz_line = ii
            if "_symmetry_equiv_pos_as_xyz_" in i:
                sym_xyz_line = ii
            if "_symmetry_equiv_pos_as_xyz_" in i:
                sym_xyz_line = ii
            if "_space_group_symop_operation_xyz_" in i:
                sym_xyz_line = ii
            if "_space_group_symop_operation_xyz" in i:
                sym_xyz_line = ii
        symm_ops = []
        terminate = False
        count = 0
        while not terminate:
            # print("sym_xyz_line", sym_xyz_line)
            tmp = lines[sym_xyz_line + count + 1]
            if "x" in tmp and "y" in tmp and "z" in tmp:
                # print("tmp", tmp)
                symm_ops.append(tmp)
                count += 1
            else:
                terminate = True
        tmp_arr = [lat_a, lat_b, lat_c, lat_alpha, lat_beta, lat_gamma]
        if any(ele == "" for ele in tmp_arr):
            raise ValueError("Lattice information is incomplete.", tmp_arr)
        lat = Lattice.from_parameters(
            lat_a, lat_b, lat_c, lat_alpha, lat_beta, lat_gamma
        )
        terminate = False
        atom_features = []
        count = 0
        beginning_atom_info_line = 0
        for ii, i in enumerate(lines):
            if "loop_" in i and "_atom_site" in lines[ii + count + 1]:
                beginning_atom_info_line = ii
        while not terminate:
            if "_atom" in lines[beginning_atom_info_line + count + 1]:
                atom_features.append(
                    lines[beginning_atom_info_line + count + 1]
                )
            count += 1
            if "_atom" not in lines[beginning_atom_info_line + count]:
                terminate = True
        terminate = False
        count = 1
        atom_liines = []
        while not terminate:
            number = beginning_atom_info_line + len(atom_features) + count
            if number == len(lines):
                terminate = True
                break
            line = lines[number]
            # print ('tis line',line)
            if len(line.split()) == len(atom_features):
                atom_liines.append(line)
                count += 1
            else:
                terminate = True
        label_index = ""
        fract_x_index = ""
        fract_y_index = ""
        fract_z_index = ""
        cartn_x_index = ""
        cartn_y_index = ""
        cartn_z_index = ""
        occupancy_index = ""
        for ii, i in enumerate(atom_features):
            if "_atom_site_label" in i:
                label_index = ii
            if "fract_x" in i:
                fract_x_index = ii
            if "fract_y" in i:
                fract_y_index = ii
            if "fract_z" in i:
                fract_z_index = ii
            if "cartn_x" in i:
                cartn_x_index = ii
            if "cartn_y" in i:
                cartn_y_index = ii
            if "cartn_z" in i:
                cartn_z_index = ii
            if "occupancy" in i:
                occupancy_index = ii
        if fract_x_index == "" and cartn_x_index == "":
            raise ValueError("Cannot find atomic coordinate info.")
        elements = []
        coords = []
        cif_atoms = None
        if fract_x_index != "":
            for ii, i in enumerate(atom_liines):
                tmp = i.split()
                tmp_lbl = list(
                    Composition.from_string(tmp[label_index]).to_dict().keys()
                )
                elem = tmp_lbl[0]
                coord = [
                    float(tmp[fract_x_index].split("(")[0]),
                    float(tmp[fract_y_index].split("(")[0]),
                    float(tmp[fract_z_index].split("(")[0]),
                ]
                if len(tmp_lbl) > 1:
                    raise ValueError("Check if labesl are correct.", tmp_lbl)
                if (
                    occupancy_index != ""
                    and not float(
                        tmp[occupancy_index].split("(")[0]
                    ).is_integer()
                ):
                    raise ValueError(
                        "Fractional occupancy is not supported.",
                        float(tmp[occupancy_index].split("(")[0]),
                        elem,
                    )

                elements.append(elem)
                coords.append(coord)
            cif_atoms = Atoms(
                lattice_mat=lat.matrix,
                elements=elements,
                coords=coords,
                cartesian=False,
            )
        elif cartn_x_index != "":
            for ii, i in enumerate(atom_liines):
                tmp = i.split()
                tmp_lbl = list(
                    Composition.from_string(tmp[label_index]).to_dict().keys()
                )
                elem = tmp_lbl[0]
                coord = [
                    float(tmp[cartn_x_index].split("(")[0]),
                    float(tmp[cartn_y_index].split("(")[0]),
                    float(tmp[cartn_z_index].split("(")[0]),
                ]
                if len(tmp_lbl) > 1:
                    raise ValueError("Check if labesl are correct.", tmp_lbl)
                if (
                    occupancy_index != ""
                    and not float(
                        tmp[occupancy_index].split("(")[0]
                    ).is_integer()
                ):
                    raise ValueError(
                        "Fractional occupancy is not supported.",
                        float(tmp[occupancy_index].split("(")[0]),
                        elem,
                    )
                elements.append(elem)
                coords.append(coord)
            cif_atoms = Atoms(
                lattice_mat=lat.matrix,
                elements=elements,
                coords=coords,
                cartesian=True,
            )
        else:
            raise ValueError(
                "Cannot find atomic coordinate info from cart or frac."
            )
        # frac_coords=list(cif_atoms.frac_coords)
        cif_elements = cif_atoms.elements
        lat = cif_atoms.lattice.matrix
        if len(symm_ops) > 1:
            frac_coords = list(cif_atoms.frac_coords)
            for i in symm_ops:
                for jj, j in enumerate(frac_coords):
                    new_c_coord = get_new_coord_for_xyz_sym(
                        xyz_string=i, frac_coord=j
                    )
                    new_frac_coord = [new_c_coord][0]
                    if not check_duplicate_coords(frac_coords, new_frac_coord):
                        frac_coords.append(new_frac_coord)
                        cif_elements.append(cif_elements[jj])
            new_atoms = Atoms(
                lattice_mat=lat,
                coords=frac_coords,
                elements=cif_elements,
                cartesian=False,
            )
            cif_atoms = new_atoms
        if get_primitive_atoms:
            cif_atoms = cif_atoms.get_primitive_atoms
        return cif_atoms

    def write_poscar(self, filename="POSCAR"):
        """Write POSCAR format file from Atoms object."""
        from jarvis.io.vasp.inputs import Poscar

        pos = Poscar(self)
        pos.write_file(filename)

    @property
    def get_xyz_string(self):
        """Get xyz string for atoms."""
        line = str(self.num_atoms) + "\n"
        line += " ".join(map(str, np.array(self.lattice_mat).flatten())) + "\n"
        for i, j in zip(self.elements, self.cart_coords):
            line += (
                str(i)
                + str(" ")
                + str(round(j[0], 4))
                + str(" ")
                + str(round(j[1], 4))
                + str(" ")
                + str(round(j[2], 4))
                + "\n"
            )
        return line

    def write_xyz(self, filename="atoms.xyz"):
        """Write XYZ format file."""
        f = open(filename, "w")
        line = str(self.num_atoms) + "\n"
        f.write(line)
        line = ",".join(map(str, np.array(self.lattice_mat).flatten())) + "\n"
        f.write(line)
        for i, j in zip(self.elements, self.cart_coords):
            f.write("%s %7.5f %7.5f %7.5f\n" % (i, j[0], j[1], j[2]))
        f.close()

    @classmethod
    def from_xyz(self, filename="dsgdb9nsd_057387.xyz", box_size=40):
        """Read XYZ file from to make Atoms object."""
        lattice_mat = [[box_size, 0, 0], [0, box_size, 0], [0, 0, box_size]]
        f = open(filename, "r")
        lines = f.read().splitlines()
        f.close()
        coords = []
        species = []
        natoms = int(lines[0])
        for i in range(natoms):
            tmp = (lines[i + 2]).split()
            coord = [(tmp[1]), (tmp[2]), (tmp[3])]
            coord = [
                0 if "*" in ii else float(ii) for ii in coord
            ]  # dsgdb9nsd_000212.xyz
            coords.append(coord)
            species.append(tmp[0])
        coords = np.array(coords)
        atoms = Atoms(
            lattice_mat=lattice_mat,
            coords=coords,
            elements=species,
            cartesian=True,
        ).center_around_origin(new_origin=[0.5, 0.5, 0.5])
        # print (atoms)
        return atoms

    @classmethod
    def from_poscar(self, filename="POSCAR"):
        """Read POSCAR/CONTCAR file from to make Atoms object."""
        from jarvis.io.vasp.inputs import Poscar

        return Poscar.from_file(filename).atoms

    @property
    def check_polar(self):
        """
        Check if the surface structure is polar.

        Comparing atom types at top and bottom.
        Applicable for sufcae with vaccums only.

        Args:

             file:atoms object (surface with vacuum)

        Returns:
               polar:True/False
        """
        up = 0
        dn = 0
        coords = np.array(self.frac_coords)
        z_max = max(coords[:, 2])
        z_min = min(coords[:, 2])
        for site, element in zip(self.frac_coords, self.elements):
            if site[2] == z_max:
                up = up + Specie(element).Z
            if site[2] == z_min:
                dn = dn + Specie(element).Z
        polar = False
        if up != dn:
            print("Seems like a polar materials.")
            polar = True
        if up == dn:
            print("Non-polar")
            polar = False
        return polar

    def strain_atoms(self, strain):
        """Apply volumetric strain to atoms."""
        # strains=np.arange(0.80,1.2,0.02)
        s = np.eye(3) + np.array(strain) * np.eye(3)
        lattice_mat = np.array(self.lattice_mat)
        lattice_mat = np.dot(lattice_mat, s)
        return Atoms(
            lattice_mat=lattice_mat,
            elements=self.elements,
            coords=self.coords,
            cartesian=self.cartesian,
        )

    def apply_strain(self, strain):
        """Apply a strain(e.g. 0.01, [0,0,.01]) to the lattice."""
        print("Use strain_atoms instead.")
        s = (1 + np.array(strain)) * np.eye(3)
        self.lattice_mat = np.dot(self.lattice_mat.T, s).T

    def to_dict(self):
        """Provide dictionary representation of the atoms object."""
        d = OrderedDict()
        d["lattice_mat"] = self.lattice_mat.tolist()
        d["coords"] = np.array(self.coords).tolist()
        d["elements"] = np.array(self.elements).tolist()
        d["abc"] = self.lattice.lat_lengths()
        d["angles"] = self.lattice.lat_angles()
        d["cartesian"] = self.cartesian
        d["props"] = self.props
        return d

    @classmethod
    def from_dict(self, d={}):
        """Form atoms object from the dictionary."""
        return Atoms(
            lattice_mat=d["lattice_mat"],
            elements=d["elements"],
            props=d["props"],
            coords=d["coords"],
            cartesian=d["cartesian"],
        )

    def remove_site_by_index(self, site=0):
        """Remove an atom by its index number."""
        new_els = []
        new_coords = []
        new_props = []
        for ii, i in enumerate(self.frac_coords):
            if ii != site:
                new_els.append(self.elements[ii])
                new_coords.append(self.frac_coords[ii])
                new_props.append(self.props[ii])
        return Atoms(
            lattice_mat=self.lattice_mat,
            elements=new_els,
            coords=np.array(new_coords),
            props=new_props,
            cartesian=False,
        )

    def add_site(
        self, element="Si", coords=[0.1, 0.1, 0.1], props=[], index=0
    ):
        """Ad an atom, coords in fractional coordinates."""
        new_els = list(self.elements)
        new_coords = list(self.frac_coords)
        new_props = list(self.props)
        new_els.insert(index, element)
        new_coords.insert(index, coords)
        new_props.insert(index, props)
        # new_els.append(element)
        # new_coords.append(coords)
        # new_props.append(props)
        new_els = np.array(new_els)
        new_coords = np.array(new_coords)
        new_props = np.array(new_props, dtype=object)
        return Atoms(
            lattice_mat=self.lattice_mat,
            elements=new_els,
            coords=np.array(new_coords),
            props=new_props,
            cartesian=False,
        )

    @property
    def get_primitive_atoms(self):
        """Get primitive Atoms using spacegroup information."""
        from jarvis.analysis.structure.spacegroup import Spacegroup3D

        return Spacegroup3D(self).primitive_atoms

    def get_all_neighbors(self, r=5, bond_tol=0.15):
        """
        Get neighbors for each atom in the unit cell, out to a distance r.

        Contains [index_i, index_j, distance, image] array.
        Adapted from pymatgen.
        """
        recp_len = np.array(self.lattice.reciprocal_lattice().abc)
        maxr = np.ceil((r + bond_tol) * recp_len / (2 * math.pi))
        nmin = np.floor(np.min(self.frac_coords, axis=0)) - maxr
        nmax = np.ceil(np.max(self.frac_coords, axis=0)) + maxr
        all_ranges = [np.arange(x, y) for x, y in zip(nmin, nmax)]
        matrix = self.lattice_mat
        neighbors = [list() for _ in range(len(self.cart_coords))]
        # all_fcoords = np.mod(self.frac_coords, 1)
        coords_in_cell = self.cart_coords  # np.dot(all_fcoords, matrix)
        site_coords = self.cart_coords
        indices = np.arange(len(site_coords))
        for image in itertools.product(*all_ranges):
            coords = np.dot(image, matrix) + coords_in_cell
            z = (coords[:, None, :] - site_coords[None, :, :]) ** 2
            all_dists = np.sum(z, axis=-1) ** 0.5
            all_within_r = np.bitwise_and(all_dists <= r, all_dists > 1e-8)
            for (j, d, within_r) in zip(indices, all_dists, all_within_r):
                for i in indices[within_r]:
                    if d[i] > bond_tol:
                        # if d[i] > bond_tol and i!=j:
                        neighbors[i].append([i, j, d[i], image])
        return np.array(neighbors, dtype="object")

    def get_neighbors_cutoffs(self, max_cut=10, r=5, bond_tol=0.15):
        neighbors = self.get_all_neighbors(r=r, bond_tol=bond_tol)
        dists = np.hstack(([[xx[2] for xx in yy] for yy in neighbors]))
        hist, bins = np.histogram(dists, bins=np.arange(0.1, 10.2, 0.1))
        # shell_vol = (
        #    4.0
        #    / 3.0
        #    * np.pi
        #    * (np.power(bins[1:], 3) - np.power(bins[:-1], 3))
        # )
        arr = []
        for i, j in zip(bins[:-1], hist):
            if j > 0:
                arr.append(i)
        rcut_buffer = 0.11
        io1 = 0
        io2 = 1
        io3 = 2
        try:
            delta = arr[io2] - arr[io1]
            while delta < rcut_buffer and arr[io2] < max_cut:
                io1 = io1 + 1
                io2 = io2 + 1
                io3 = io3 + 1
                delta = arr[io2] - arr[io1]

            rcut1 = (arr[io2] + arr[io1]) / float(2.0)
        except Exception:
            print("Warning:Setting first nbr cut-off as minimum bond-dist")
            rcut1 = arr[0]
            pass
        try:
            delta = arr[io3] - arr[io2]
            while (
                delta < rcut_buffer
                and arr[io3] < max_cut
                and arr[io2] < max_cut
            ):
                io2 = io2 + 1
                io3 = io3 + 1
                delta = arr[io3] - arr[io2]
            rcut2 = float(arr[io3] + arr[io2]) / float(2.0)
        except Exception:
            print("Warning:Setting first and second nbr cut-off equal")
            print("You might consider increasing max_n parameter")
            rcut2 = rcut1
            pass
        return rcut1, rcut2, neighbors

    def atomwise_angle_and_radial_distribution(
        self, r=5, bond_tol=0.15, c_size=10, verbose=False
    ):
        rcut1, rcut2, neighbors = self.get_neighbors_cutoffs(
            r=r, bond_tol=bond_tol
        )
        from jarvis.analysis.structure.neighbors import NeighborsAnalysis
        from jarvis.core.utils import bond_angle as angle

        nbor_info = NeighborsAnalysis(
            self, rcut1=rcut1, rcut2=rcut2
        ).nbor_list(rcut=rcut1, c_size=c_size)
        nat = nbor_info["nat"]
        dist = nbor_info["dist"]
        atom_rdfs = []
        nbins = 180
        actual_prdf = []
        # actual_pangs = []
        actual_pangs = np.zeros((nat, 380))
        for i in range(nat):
            hist, bins = np.histogram(
                dist[:, i], bins=np.arange(0.1, rcut1 + 0.2, 0.1)
            )
            actual_prdf.append(dist[:, i])
            atom_rdfs.append(hist.tolist())
            if verbose:
                exact_dists = np.arange(0.1, rcut1 + 0.2, 0.1)[hist.nonzero()]
                print("exact_dists", exact_dists)
        # prdf_arr = np.array(atom_rdfs)[0 : self.num_atoms]

        atom_angles = []
        for i in range(nbor_info["nat"]):
            angles = [
                angle(
                    nbor_info["dist"][in1][i],
                    nbor_info["dist"][in2][i],
                    nbor_info["bondx"][in1][i],
                    nbor_info["bondx"][in2][i],
                    nbor_info["bondy"][in1][i],
                    nbor_info["bondy"][in2][i],
                    nbor_info["bondz"][in1][i],
                    nbor_info["bondz"][in2][i],
                )
                for in1 in range(nbor_info["nn"][i])
                for in2 in range(nbor_info["nn"][i])
                if in2 > in1
                and nbor_info["dist"][in1][i] * nbor_info["dist"][in2][i] != 0
            ]
            ang_hist, ang_bins = np.histogram(
                angles, bins=np.arange(1, nbins + 2, 1), density=False,
            )
            for jj, j in enumerate(angles):
                actual_pangs[i, jj] = j
            # actual_pangs.append(angles)
            atom_angles.append(ang_hist)
            if verbose:
                exact_angles = np.arange(1, nbins + 2, 1)[ang_hist.nonzero()]
                print("exact_angles", exact_angles)
        # return (atom_angles)#/nbor_info['nat']

        # pangle_arr = np.array(atom_angles)[0 : self.num_atoms]
        return (
            neighbors,
            np.array(atom_rdfs)[0 : self.num_atoms],
            np.array(atom_angles)[0 : self.num_atoms],
            np.array(actual_prdf[0 : self.num_atoms]),
            np.array(actual_pangs[0 : self.num_atoms]),
            nbor_info,
        )

    @property
    def raw_distance_matrix(self):
        """Provide distance matrix."""
        coords = np.array(self.cart_coords)
        z = (coords[:, None, :] - coords[None, :, :]) ** 2
        return np.sum(z, axis=-1) ** 0.5

    @property
    def raw_angle_matrix(self, cut_off=5.0):
        """Provide distance matrix."""
        coords = np.array(self.cart_coords)
        angles = []
        for a, b, c in itertools.product(coords, coords, coords):
            if (
                not np.array_equal(a, b)
                and not np.array_equal(b, c)
                and np.linalg.norm((a - b)) < cut_off
                and np.linalg.norm((c - b)) < cut_off
            ):
                angle = get_angle(a, b, c)
                angles.append(angle)
        return angles

    def center(self, axis=2, vacuum=18.0, about=None):
        """
        Center structure with vacuum padding.

        Args:
          vacuum:vacuum size

          axis: direction
        """
        cell = self.lattice_mat
        p = self.cart_coords

        dirs = np.zeros_like(cell)
        for i in range(3):
            dirs[i] = np.cross(cell[i - 1], cell[i - 2])
            dirs[i] /= np.sqrt(np.dot(dirs[i], dirs[i]))  # normalize
            if np.dot(dirs[i], cell[i]) < 0.0:
                dirs[i] *= -1

        if isinstance(axis, int):
            axes = (axis,)
        else:
            axes = axis

        # if vacuum and any(self.pbc[x] for x in axes):
        #     warnings.warn(
        #         'You are adding vacuum along a periodic direction!')

        # Now, decide how much each basis vector should be made longer
        longer = np.zeros(3)
        shift = np.zeros(3)
        for i in axes:
            p0 = np.dot(p, dirs[i]).min() if len(p) else 0
            p1 = np.dot(p, dirs[i]).max() if len(p) else 0
            height = np.dot(cell[i], dirs[i])
            if vacuum is not None:
                lng = (p1 - p0 + 2 * vacuum) - height
            else:
                lng = 0.0  # Do not change unit cell size!
            top = lng + height - p1
            shf = 0.5 * (top - p0)
            cosphi = np.dot(cell[i], dirs[i]) / np.sqrt(
                np.dot(cell[i], cell[i])
            )
            longer[i] = lng / cosphi
            shift[i] = shf / cosphi

        # Now, do it!
        translation = np.zeros(3)
        for i in axes:
            nowlen = np.sqrt(np.dot(cell[i], cell[i]))
            if vacuum is not None or cell[i].any():
                cell[i] = cell[i] * (1 + longer[i] / nowlen)
                translation += shift[i] * cell[i] / nowlen

        new_coords = p + translation
        if about is not None:
            for vector in cell:
                new_coords -= vector / 2.0
            new_coords += about
        atoms = Atoms(
            lattice_mat=cell,
            elements=self.elements,
            coords=new_coords,
            cartesian=True,
        )
        return atoms

    @property
    def volume(self):
        """Get volume of the atoms object."""
        m = self.lattice_mat
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol

    @property
    def composition(self):
        """Get composition of the atoms object."""
        comp = {}
        for i in self.elements:
            comp[i] = comp.setdefault(i, 0) + 1
        return Composition(OrderedDict(comp))

    @property
    def density(self):
        """Get density in g/cm3 of the atoms object."""
        den = float(self.composition.weight * amu_gm) / (
            float(self.volume) * (ang_cm) ** 3
        )
        return den

    @property
    def atomic_numbers(self):
        """Get list of atomic numbers of atoms in the atoms object."""
        numbers = []
        for i in self.elements:
            numbers.append(Specie(i).Z)
        return numbers

    @property
    def num_atoms(self):
        """Get number of atoms."""
        return len(self.coords)

    @property
    def uniq_species(self):
        """Get unique elements."""
        uniq = set()
        return [x for x in self.elements if not (x in uniq or uniq.add(x))]

    def get_center_of_mass(self):
        """Get center of mass of the atoms object."""
        # atomic_mass
        m = []
        for i in self.elements:
            m.append(Specie(i).atomic_mass)
        m = np.array(m)
        com = np.dot(m, self.cart_coords) / m.sum()
        # com = np.linalg.solve(self.lattice_mat.T, com)
        return com

    def get_origin(self):
        """Get center of mass of the atoms object."""
        # atomic_mass
        return self.frac_coords.mean(axis=0)

    def center_around_origin(self, new_origin=[0.0, 0.0, 0.5]):
        """Center around given origin."""
        lat = self.lattice_mat
        typ_sp = self.elements
        natoms = self.num_atoms
        # abc = self.lattice.lat_lengths()
        COM = self.get_origin()
        # COM = self.get_center_of_mass()
        x = np.zeros((natoms))
        y = np.zeros((natoms))
        z = np.zeros((natoms))
        coords = list()
        for i in range(0, natoms):
            # cart_coords[i]=self.cart_coords[i]-COM
            x[i] = self.frac_coords[i][0] - COM[0] + new_origin[0]
            y[i] = self.frac_coords[i][1] - COM[1] + new_origin[1]
            z[i] = self.frac_coords[i][2] - COM[2] + new_origin[2]
            coords.append([x[i], y[i], z[i]])
        struct = Atoms(
            lattice_mat=lat, elements=typ_sp, coords=coords, cartesian=False
        )
        return struct

    def pymatgen_converter(self):
        """Get pymatgen representation of the atoms object."""
        try:
            from pymatgen.core.structure import Structure

            return Structure(
                self.lattice_mat,
                self.elements,
                self.frac_coords,
                coords_are_cartesian=False,
            )
        except Exception:
            print("Requires pymatgen for this functionality.")
            pass

    def phonopy_converter(self, pbc=True):
        """Get phonopy representation of the atoms object."""
        try:
            from phonopy.structure.atoms import Atoms as PhonopyAtoms

            return PhonopyAtoms(
                symbols=self.elements,
                positions=self.cart_coords,
                pbc=pbc,
                cell=self.lattice_mat,
            )
        except Exception:
            print("Requires phonopy for this functionality.")
            pass

    def ase_converter(self, pbc=True):
        """Get ASE representation of the atoms object."""
        try:
            from ase import Atoms as AseAtoms

            return AseAtoms(
                symbols=self.elements,
                positions=self.cart_coords,
                pbc=pbc,
                cell=self.lattice_mat,
            )
        except Exception:
            print("Requires ASE for this functionality.")
            pass

    def spacegroup(self, symprec=1e-3):
        """Get spacegroup of the atoms object."""
        import spglib

        sg = spglib.get_spacegroup(
            (self.lattice_mat, self.frac_coords, self.atomic_numbers),
            symprec=symprec,
        )
        return sg

    @property
    def packing_fraction(self):
        """Get packing fraction of the atoms object."""
        total_rad = 0
        for i in self.elements:
            total_rad = total_rad + Specie(i).atomic_rad ** 3
        pf = np.array([4 * np.pi * total_rad / (3 * self.volume)])
        return round(pf[0], 5)

    def lattice_points_in_supercell(self, supercell_matrix):
        """
        Adapted from Pymatgen.

        Returns the list of points on the original lattice contained in the
        supercell in fractional coordinates (with the supercell basis).
        e.g. [[2,0,0],[0,1,0],[0,0,1]] returns [[0,0,0],[0.5,0,0]]

        Args:

            supercell_matrix: 3x3 matrix describing the supercell

        Returns:
            numpy array of the fractional coordinates
        """
        diagonals = np.array(
            [
                [0, 0, 0],
                [0, 0, 1],
                [0, 1, 0],
                [0, 1, 1],
                [1, 0, 0],
                [1, 0, 1],
                [1, 1, 0],
                [1, 1, 1],
            ]
        )
        d_points = np.dot(diagonals, supercell_matrix)

        mins = np.min(d_points, axis=0)
        maxes = np.max(d_points, axis=0) + 1

        ar = (
            np.arange(mins[0], maxes[0])[:, None]
            * np.array([1, 0, 0])[None, :]
        )
        br = (
            np.arange(mins[1], maxes[1])[:, None]
            * np.array([0, 1, 0])[None, :]
        )
        cr = (
            np.arange(mins[2], maxes[2])[:, None]
            * np.array([0, 0, 1])[None, :]
        )

        all_points = ar[:, None, None] + br[None, :, None] + cr[None, None, :]
        all_points = all_points.reshape((-1, 3))

        frac_points = np.dot(all_points, np.linalg.inv(supercell_matrix))

        tvects = frac_points[
            np.all(frac_points < 1 - 1e-10, axis=1)
            & np.all(frac_points >= -1e-10, axis=1)
        ]
        assert len(tvects) == round(abs(np.linalg.det(supercell_matrix)))
        return tvects

    def make_supercell_matrix(self, scaling_matrix):
        """
        Adapted from Pymatgen.

        Makes a supercell. Allowing to have sites outside the unit cell.

        Args:
            scaling_matrix: A scaling matrix for transforming the lattice
            vectors. Has to be all integers. Several options are possible:
            a. A full 3x3 scaling matrix defining the linear combination
             the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
             1]] generates a new structure with lattice vectors a' =
             2a + b, b' = 3b, c' = c where a, b, and c are the lattice
             vectors of the original structure.
            b. An sequence of three scaling factors. E.g., [2, 1, 1]
             specifies that the supercell should have dimensions 2a x b x
             c.
            c. A number, which simply scales all lattice vectors by the
             same factor.

        Returns:
            Supercell structure. Note that a Structure is always returned,
            even if the input structure is a subclass of Structure. This is
            to avoid different arguments signatures from causing problems. If
            you prefer a subclass to return its own type, you need to override
            this method in the subclass.
        """
        scale_matrix = np.array(scaling_matrix, np.int16)
        if scale_matrix.shape != (3, 3):
            scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
        new_lattice = Lattice(np.dot(scale_matrix, self.lattice_mat))

        f_lat = self.lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.cart_coords(f_lat)

        new_sites = []
        new_elements = []
        for site, el in zip(self.cart_coords, self.elements):
            for v in c_lat:
                new_elements.append(el)
                tmp = site + v
                new_sites.append(tmp)
        return Atoms(
            lattice_mat=new_lattice.lattice(),
            elements=new_elements,
            coords=new_sites,
            cartesian=True,
        )

    def make_supercell(self, dim=[2, 2, 2]):
        """Make supercell of dimension dim."""
        dim = np.array(dim)
        if dim.shape == (3, 3):
            dim = np.array([int(np.linalg.norm(v)) for v in dim])
        coords = self.frac_coords
        all_symbs = self.elements  # [i.symbol for i in s.species]
        nat = len(coords)

        new_nat = nat * dim[0] * dim[1] * dim[2]
        new_coords = np.zeros((new_nat, 3))
        new_symbs = []  # np.chararray((new_nat))
        props = []  # self.props

        ct = 0
        for i in range(nat):
            for j in range(dim[0]):
                for k in range(dim[1]):
                    for m in range(dim[2]):
                        props.append(self.props[i])
                        new_coords[ct][0] = (coords[i][0] + j) / float(dim[0])
                        new_coords[ct][1] = (coords[i][1] + k) / float(dim[1])
                        new_coords[ct][2] = (coords[i][2] + m) / float(dim[2])
                        new_symbs.append(all_symbs[i])
                        ct = ct + 1

        nat = new_nat

        nat = len(coords)  # int(s.composition.num_atoms)
        lat = np.zeros((3, 3))
        box = self.lattice_mat
        lat[0][0] = dim[0] * box[0][0]
        lat[0][1] = dim[0] * box[0][1]
        lat[0][2] = dim[0] * box[0][2]
        lat[1][0] = dim[1] * box[1][0]
        lat[1][1] = dim[1] * box[1][1]
        lat[1][2] = dim[1] * box[1][2]
        lat[2][0] = dim[2] * box[2][0]
        lat[2][1] = dim[2] * box[2][1]
        lat[2][2] = dim[2] * box[2][2]
        super_cell = Atoms(
            lattice_mat=lat,
            coords=new_coords,
            elements=new_symbs,
            props=props,
            cartesian=False,
        )
        return super_cell

    def get_lll_reduced_structure(self):
        """Get LLL algorithm based reduced structure."""
        reduced_latt = self.lattice.get_lll_reduced_lattice()
        if reduced_latt != self.lattice:
            return Atoms(
                lattice_mat=reduced_latt._lat,
                elements=self.elements,
                coords=self.frac_coords,
                cartesian=False,
            )
        else:
            return Atoms(
                lattice_mat=self.lattice_mat,
                elements=self.elements,
                coords=self.frac_coords,
                cartesian=False,
            )

    def __repr__(self):
        """Get representation during print statement."""
        return self.get_string()

    def get_string(self, cart=True, sort_order="X"):
        """
        Convert Atoms to string.

        Optional arguments below.

        Args:
          cart:True/False for cartesian/fractional coords.

          sort_order: sort by chemical properties of
                    elements. Default electronegativity.
        """
        system = str(self.composition.formula)
        header = (
            str(system)
            + str("\n1.0\n")
            + str(self.lattice_mat[0][0])
            + " "
            + str(self.lattice_mat[0][1])
            + " "
            + str(self.lattice_mat[0][2])
            + "\n"
            + str(self.lattice_mat[1][0])
            + " "
            + str(self.lattice_mat[1][1])
            + " "
            + str(self.lattice_mat[1][2])
            + "\n"
            + str(self.lattice_mat[2][0])
            + " "
            + str(self.lattice_mat[2][1])
            + " "
            + str(self.lattice_mat[2][2])
            + "\n"
        )
        if sort_order is None:
            order = np.argsort(self.elements)
        else:
            order = np.argsort(
                [Specie(i).element_property(sort_order) for i in self.elements]
            )
        if cart:
            coords = self.cart_coords
        else:
            coords = self.frac_coords
        coords_ordered = np.array(coords)[order]
        elements_ordered = np.array(self.elements)[order]
        props_ordered = np.array(self.props)[order]
        # check_selective_dynamics = False # TODO
        counts = get_counts(elements_ordered)
        cart_frac = ""
        if cart:
            cart_frac = "\nCartesian\n"
        else:
            cart_frac = "\nDirect\n"

        if "T" in "".join(map(str, self.props[0])):
            middle = (
                " ".join(map(str, counts.keys()))
                + "\n"
                + " ".join(map(str, counts.values()))
                + "\nSelective dynamics\n"
                + cart_frac
            )
        else:
            middle = (
                " ".join(map(str, counts.keys()))
                + "\n"
                + " ".join(map(str, counts.values()))
                + cart_frac
            )
        rest = ""
        if coords_ordered.ndim == 1:
            coords_ordered = np.array([coords])
        for ii, i in enumerate(coords_ordered):
            if self.show_props:
                rest = (
                    rest
                    + " ".join(map(str, i))
                    + " "
                    + str(props_ordered[ii])
                    + "\n"
                )
            else:
                rest = rest + " ".join(map(str, i)) + "\n"

        result = header + middle + rest
        return result


class VacuumPadding(object):
    """Adds vaccum padding to make 2D structure or making molecules."""

    def __init__(self, atoms, vacuum=20.0):
        """
        Initialize object.

        Args:
            atoms: Atoms object
            vacuum:  vacuum padding
        """
        self.atoms = atoms
        self.vacuum = vacuum

    def get_effective_2d_slab(self):
        """Add 2D vacuum to a system."""
        z_coord = []
        for i in self.atoms.frac_coords:
            tmp = i[2]
            if i[2] >= 0.5:
                tmp = i[2] - 1
            elif i[2] < -0.5:
                tmp = i[2] + 1
            z_coord.append(tmp)

        zmaxp = max(np.array(z_coord) * self.atoms.lattice_mat[2][2])
        zminp = min(np.array(z_coord) * self.atoms.lattice_mat[2][2])
        thickness = abs(zmaxp - zminp)
        padding = self.vacuum + thickness

        # lattice_mat = self.atoms.lattice_mat
        # lattice_mat[2][2] = padding
        # elements = self.atoms.elements
        # atoms = atoms.center(vacuum=0.0)
        new_lat = self.atoms.lattice_mat
        a1 = new_lat[0]
        a2 = new_lat[1]
        a3 = new_lat[2]
        new_lat = np.array(
            [
                a1,
                a2,
                np.cross(a1, a2)
                * np.dot(a3, np.cross(a1, a2))
                / np.linalg.norm(np.cross(a1, a2)) ** 2,
            ]
        )

        a1 = new_lat[0]
        a2 = new_lat[1]
        a3 = new_lat[2]

        latest_lat = np.array(
            [
                (np.linalg.norm(a1), 0, 0),
                (
                    np.dot(a1, a2) / np.linalg.norm(a1),
                    np.sqrt(
                        np.linalg.norm(a2) ** 2
                        - (np.dot(a1, a2) / np.linalg.norm(a1)) ** 2
                    ),
                    0,
                ),
                (0, 0, np.linalg.norm(a3)),
            ]
        )
        # latest_lat[2][2] = padding
        M = np.linalg.solve(new_lat, latest_lat)
        new_cart_coords = self.atoms.cart_coords
        new_coords = np.dot(new_cart_coords, M)
        new_atoms = Atoms(
            lattice_mat=latest_lat,
            elements=self.atoms.elements,
            coords=new_coords,
            cartesian=True,
        )
        frac_coords = new_atoms.frac_coords
        frac_coords[:] = frac_coords[:] % 1
        new_atoms = Atoms(
            lattice_mat=latest_lat,
            elements=self.atoms.elements,
            coords=frac_coords,
            cartesian=False,
        )
        new_lat = new_atoms.lattice_mat
        new_cart_coords = new_atoms.cart_coords
        elements = new_atoms.elements
        new_lat[2][2] = padding  # new_lat[2][2]+ self.vacuum
        with_vacuum_atoms = Atoms(
            lattice_mat=new_lat,
            elements=elements,
            coords=new_cart_coords,
            cartesian=True,
        )
        frac = np.array(with_vacuum_atoms.frac_coords)
        frac[:, 2] = frac[:, 2] - np.mean(frac[:, 2]) + 0.5
        frac[:, 2] = frac[:, 2] - np.mean(frac[:, 2]) + 0.5
        with_vacuum_atoms = Atoms(
            lattice_mat=new_lat,
            elements=elements,
            coords=frac,
            cartesian=False,
        )
        return with_vacuum_atoms

    def get_effective_molecule(self):
        """Add vacuum around a system."""
        x_coord = []
        y_coord = []
        z_coord = []
        for i in self.atoms.frac_coords:
            tmp = i[0]
            if i[0] >= 0.5:
                tmp = i[0] - 1
            elif i[0] < -0.5:
                tmp = i[0] + 1
            x_coord.append(tmp)
            tmp = i[1]
            if i[1] >= 0.5:
                tmp = i[1] - 1
            elif i[1] < -0.5:
                tmp = i[1] + 1
            y_coord.append(tmp)
            tmp = i[2]
            if i[2] >= 0.5:
                tmp = i[2] - 1
            elif i[2] < -0.5:
                tmp = i[2] + 1
            z_coord.append(tmp)

        xmaxp = max(np.array(x_coord) * self.atoms.lattice_mat[0][0])
        xminp = min(np.array(x_coord) * self.atoms.lattice_mat[0][0])
        ymaxp = max(np.array(y_coord) * self.atoms.lattice_mat[1][1])
        yminp = min(np.array(y_coord) * self.atoms.lattice_mat[1][1])
        zmaxp = max(np.array(z_coord) * self.atoms.lattice_mat[2][2])
        zminp = min(np.array(z_coord) * self.atoms.lattice_mat[2][2])
        thickness_x = abs(xmaxp - xminp)
        thickness_y = abs(ymaxp - yminp)
        thickness_z = abs(zmaxp - zminp)

        lattice_mat = np.zeros((3, 3))  # self.atoms.lattice_mat
        lattice_mat[0][0] = self.vacuum + thickness_x
        lattice_mat[1][1] = self.vacuum + thickness_y
        lattice_mat[2][2] = self.vacuum + thickness_z
        elements = self.atoms.elements
        atoms = Atoms(
            lattice_mat=lattice_mat,
            coords=self.atoms.cart_coords,
            elements=elements,
            cartesian=True,
        )
        frac = np.array(atoms.frac_coords)
        frac[:, 0] = frac[:, 0] - np.mean(frac[:, 0]) + 0.5
        frac[:, 1] = frac[:, 1] - np.mean(frac[:, 1]) + 0.5
        frac[:, 2] = frac[:, 2] - np.mean(frac[:, 2]) + 0.5
        with_vacuum_atoms = Atoms(
            lattice_mat=lattice_mat,
            elements=elements,
            coords=frac,
            cartesian=False,
        )
        return with_vacuum_atoms


def fix_pbc(atoms):
    """Use for making Atoms with vacuum."""
    new_f_coords = []
    for i in atoms.frac_coords:
        if i[2] > 0.5:
            i[2] = i[2] - 1
        if i[2] < -0.5:
            i[2] = i[2] + 1
        new_f_coords.append(i)
    return Atoms(
        lattice_mat=atoms.lattice_mat,
        elements=atoms.elements,
        coords=new_f_coords,
        cartesian=False,
    )


def add_atoms(top, bottom, distance=[0, 0, 1], apply_strain=False):
    """
    Add top and bottom Atoms with a distance array.

    Bottom Atoms lattice-matrix is chosen as final lattice.
    """
    top = top.center_around_origin([0, 0, 0])
    bottom = bottom.center_around_origin(distance)
    strain_x = (
        top.lattice_mat[0][0] - bottom.lattice_mat[0][0]
    ) / bottom.lattice_mat[0][0]
    strain_y = (
        top.lattice_mat[1][1] - bottom.lattice_mat[1][1]
    ) / bottom.lattice_mat[1][1]
    if apply_strain:
        top.apply_strain([strain_x, strain_y, 0])
    #  print("strain_x,strain_y", strain_x, strain_y)
    elements = []
    coords = []
    lattice_mat = bottom.lattice_mat
    for i, j in zip(bottom.elements, bottom.frac_coords):
        elements.append(i)
        coords.append(j)
    top_cart_coords = lattice_coords_transformer(
        new_lattice_mat=top.lattice_mat,
        old_lattice_mat=bottom.lattice_mat,
        cart_coords=top.cart_coords,
    )
    top_frac_coords = bottom.lattice.frac_coords(top_cart_coords)
    for i, j in zip(top.elements, top_frac_coords):
        elements.append(i)
        coords.append(j)

    order = np.argsort(np.array(elements))
    elements = np.array(elements)[order]
    coords = np.array(coords)[order]
    determnt = np.linalg.det(np.array(lattice_mat))
    if determnt < 0.0:
        lattice_mat = -1 * np.array(lattice_mat)
    determnt = np.linalg.det(np.array(lattice_mat))
    if determnt < 0.0:
        print("Serious issue, check lattice vectors.")
        print("Many software follow right hand basis rule only.")
    combined = Atoms(
        lattice_mat=lattice_mat,
        coords=coords,
        elements=elements,
        cartesian=False,
    ).center_around_origin()
    return combined


def get_supercell_dims(atoms, enforce_c_size=10, extend=1):
    """Get supercell dimensions."""
    a = atoms.lattice.lat_lengths()[0]
    b = atoms.lattice.lat_lengths()[1]
    c = atoms.lattice.lat_lengths()[2]
    dim1 = int(float(enforce_c_size) / float(a)) + extend
    dim2 = int(float(enforce_c_size) / float(b)) + extend
    dim3 = int(float(enforce_c_size) / float(c)) + extend
    return [dim1, dim2, dim3]


def pmg_to_atoms(pmg=""):
    """Convert pymatgen structure to Atoms."""
    return Atoms(
        lattice_mat=pmg.lattice.matrix,
        elements=[i.symbol for i in pmg.species],
        coords=pmg.frac_coords,
        cartesian=False,
    )


def ase_to_atoms(ase_atoms=""):
    """Convert ase structure to Atoms."""
    return Atoms(
        lattice_mat=ase_atoms.get_cell(),
        elements=ase_atoms.get_chemical_symbols(),
        coords=ase_atoms.get_positions(),
        #         pbc=True,
    )


def crop_square(atoms=None, csize=10):
    """Crop a sqaur portion from a surface/2D system."""
    sz = csize / 2
    # just to make sure we have enough material to crop from
    enforce_c_size = sz * 3
    dims = get_supercell_dims(atoms, enforce_c_size=enforce_c_size)
    b = atoms.make_supercell_matrix(dims).center_around_origin()
    lat_mat = [
        [enforce_c_size, 0, 0],
        [0, enforce_c_size, 0],
        [0, 0, b.lattice_mat[2][2]],
    ]
    # M = np.linalg.solve(b.lattice_mat, lat_mat)
    # tol = 3

    els = []
    coords = []
    for i, j in zip(b.cart_coords, b.elements):
        if i[0] <= sz and i[0] >= -sz and i[1] <= sz and i[1] >= -sz:
            els.append(j)
            coords.append(i)
    coords = np.array(coords)
    # new_mat = (
    #    [max(coords[:, 0]) - min(coords[:, 0]) + tol, 0, 0],
    #    [0, max(coords[:, 1]) - min(coords[:, 1]) + tol, 0],
    #    [0, 0, b.lattice_mat[2][2]],
    # )
    new_atoms = Atoms(
        lattice_mat=lat_mat, elements=els, coords=coords, cartesian=True
    ).center_around_origin([0.5, 0.5, 0.5])
    return new_atoms


def build_xanes_poscar(
    atoms=None,
    selected_element="Si",
    prefix="-",
    extend=1,
    enforce_c_size=12,
    dir=".",
    filename_with_prefix=False,
):
    """Generate POSCAR file for XANES, note the element ordering."""
    from jarvis.core.utils import rand_select
    from jarvis.analysis.structure.spacegroup import Spacegroup3D

    dims = get_supercell_dims(
        atoms, enforce_c_size=enforce_c_size, extend=extend
    )
    atoms = atoms.make_supercell_matrix(dims)
    spath = os.path.join(dir, "POSCAR-supercell.vasp")
    atoms.write_poscar(spath)
    spg = Spacegroup3D(atoms)
    wyckoffs = spg._dataset["wyckoffs"]
    atoms.props = wyckoffs
    props = rand_select(atoms.props)
    tmp_atoms = atoms
    for i, j in props.items():
        if tmp_atoms.elements[j] == selected_element:
            atoms = tmp_atoms
            defect_strt = atoms.remove_site_by_index(j)
            coords = tmp_atoms.frac_coords[j]
            added_strt = defect_strt.add_site(element="XANES", coords=coords)
            if filename_with_prefix:
                filename = (
                    "POSCAR"
                    + prefix
                    + tmp_atoms.elements[j]
                    + "_"
                    + str(j)
                    + "_"
                    + tmp_atoms.props[j]
                    + ".vasp"
                )
            else:
                filename = "POSCAR"

            filename = os.path.join(dir, filename)
            added_strt.props = ["" for i in range(added_strt.num_atoms)]
            added_strt.write_poscar(filename)
            f = open(filename, "r")
            filedata = f.read()
            f.close()

            newdata = filedata.replace("XANES", tmp_atoms.elements[j])
            f = open(filename, "w")
            f.write(newdata)
            f.close()


# ['Mn ', 'Mn ', 'Ru ', 'U ']
#
# def clear_elements(atoms=None):
# return info
#    info={}
#    info['lattice_mat']=atoms['lattice_mat']
#    info['coords']=atoms['coords']
#    info['props']=atoms['props']
#    info['cartesian']=atoms['cartesian']
#    elements=[i.strip() for i in atoms['elements']]
#    info['elements']=elements
#    return info

"""
if __name__ == "__main__":
    x=Atoms.from_cif('1000000.cif')
    f=open('1000000.cif','r')
    lines=f.read()
    f.close()
    x=Atoms.from_cif(from_string=lines)
    print (x)
    print (x.num_atoms)
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Si.write_xyz("atoms.xyz")
    from jarvis.io.vasp.inputs import Poscar

    Si = Atoms.from_poscar("/users/knc6/POSCAR")
    Si.write_cif()
    a = Atoms.from_cif("atoms.cif")
    print(a)
    fn = "/cluster/users/knc6/justback/desc_library/cod/cif/1000052.cif"
    # fn="/cluster/users/knc6/justback/desc_library/cod/cif/1000443.cif"
    a = Atoms.from_cif(filename=fn)
    print (a)
    Si.write_poscar()
    print (Si.composition.reduced_formula)
    #print (Si.get_string())
    print (Si.get_primitive_atoms)
    print (Si.raw_distance_matrix)
    import sys
    sys.exit()
    #print (Si.props)
    #print (Si.make_supercell().props)
    d = Si.to_dict()
    print (d)
    a=Atoms.from_dict(d)
    print (a)
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Si.props = ["a", "a"]
    # spg = Spacegroup3D(Si)
    #polar=Si.check_polar
    #print ('polar',polar)
    # Si = spg.conventional_standard_structure
    # print ('center',Si.center())
    # print ('propos',Si.props)
    #print("Supercell\n", Si.make_supercell([2, 2, 2]))
    # print (Si.make_supercell().props)
    #print(Si.make_supercell([2, 2, 2]).remove_site_by_index())
    print ('Si',Si)
    # print ('reduced',Si.get_lll_reduced_structure())
    # print ('pf',Si.packing_fraction,Si.make_supercell())
    pmg = Si.pymatgen_converter()
    pmg.make_supercell([2, 2, 2])
    #print (pmg)
    # print (Si.get_center_of_mass())
    # print (Si.get_string())
    a=Atoms.from_cif('ll.cif')
    print(a)
    from pymatgen.core.structure import Structure
    s=Structure.from_file('ll.cif')
    from pymatgen.io.vasp.inputs import Poscar
    p=Poscar(s)
    print (p)
"""
