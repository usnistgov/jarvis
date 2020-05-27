"""
This module provides classes to specify atomic structure
"""
from collections import Counter
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie
from jarvis.core.lattice import Lattice
import matplotlib.pyplot as plt
from collections import OrderedDict
import pprint
import math
from numpy.linalg import norm, solve
from jarvis.core.utils import get_counts
amu_gm = 1.66054e-24
ang_cm = 1e-8



class Atoms(object):
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
        
        Create atomic structure with lattice, coordinates, atom type and other information
        
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
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements,cartesian=True)
        >>> round(Si.density,2)
        2.33
        >>> Si.spacegroup()
        'C2/m (12)'
        >>> Si.pymatgen_converter()!={}
        True
        
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
        if self.cartesian == True:
            self.cart_coords = self.coords
            self.frac_coords = np.array(self.lattice.frac_coords(self.coords))
            # print ('TRUE')
        else:
            self.frac_coords = self.coords
            self.cart_coords = np.array(self.lattice.cart_coords(self.coords))
            # print ('FALSE')

    @property
    def check_polar(self):
        """
            Check if the surface structure is polar
            by comparing atom types at top and bottom.
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
            print("polar")
            polar = True
        if up == dn:
            print("Non-polar")
            polar = False
        return polar

    def apply_strain(self, strain):
        """
        Apply a strain(e.g. 0.01) to the lattice.
        """
        s = (1 + np.array(strain)) * np.eye(3)
        self.lattice_mat = np.dot(self.lattice_mat.T, s).T

    def to_dict(self):
        """
        Dictionary representation of the atoms object
        """
        d = OrderedDict()
        d["lattice_mat"] = self.lattice_mat.tolist()
        d["coords"] = np.array(self.coords).tolist()
        d["elements"] = self.elements
        d["abc"] = self.lattice.lat_lengths()
        d["angles"] = self.lattice.lat_angles()
        d["cartesian"] = self.cartesian
        d["props"] = self.props
        return d

    @classmethod
    def from_dict(self, d={}):
        """
        Form atoms object from the dictionary
        """
        return Atoms(
            lattice_mat=d["lattice_mat"],
            elements=d["elements"],
            props=d["props"],
            coords=d["coords"],
            cartesian=d["cartesian"],
        )

    def remove_site_by_index(self, site=0):
        """
        Remove an atom by its index number
        """
        new_els = []
        new_coords = []
        new_props = []
        for ii, i in enumerate(self.frac_coords):
            if ii != site:
                # print(self.elements, len(self.elements), len(self.frac_coords))
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
    
    @property
    def get_primitive_atoms(self):
        from jarvis.analysis.structure.spacegroup import Spacegroup3D
        return Spacegroup3D(self).primitive_atoms

    @property
    def raw_distance_matrix(self):
        coords = np.array(self.cart_coords)
        z = (coords[:, None, :] - coords[None, :, :]) ** 2
        return np.sum(z, axis=-1) ** 0.5

    def center(self, axis=2, vacuum=18.0, about=None):
        """
        Center structure with vacuum padding os size:vacuum in a direction:axis
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
            cosphi = np.dot(cell[i], dirs[i]) / np.sqrt(np.dot(cell[i], cell[i]))
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
            lattice_mat=cell, elements=self.elements, coords=new_coords, cartesian=True
        )
        return atoms

    @property
    def volume(self):
        """
        Get volume of the atoms object
        """
        m = self.lattice_mat
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol

    @property
    def composition(self):
        """
        Get composition of the atoms object
        """
        comp = {}
        for i in self.elements:
            comp[i] = comp.setdefault(i, 0) + 1
        return Composition(OrderedDict(comp))

    @property
    def density(self):
        """
        Get density in g/cm3 of the atoms object
        """
        den = float(self.composition.weight * amu_gm) / (
            float(self.volume) * (ang_cm) ** 3
        )
        return den

    @property
    def atomic_numbers(self):
        """
        Get list of atomic numbers of atoms in the atoms object
        """
        numbers = []
        for i in self.elements:
            numbers.append(Specie(i).Z)
        return numbers

    @property
    def num_atoms(self):
        """
        Get number of atoms
        """
        return len(self.coords)

    def get_center_of_mass(self):
        """
        Get center of mass of the atoms object
        """
        # atomic_mass
        m = []
        for i in self.elements:
            m.append(Specie(i).atomic_mass)
        m = np.array(m)
        com = np.dot(m, self.cart_coords) / m.sum()
        # com = np.linalg.solve(self.lattice_mat.T, com)
        return com

    def get_origin(self):
        """
        Get center of mass of the atoms object
        """
        # atomic_mass
        return self.frac_coords.mean(axis=0)

    def center_around_origin(self, new_origin=[0.0, 0.0, 0.5]):
        lat = self.lattice_mat
        typ_sp = self.elements
        natoms = self.num_atoms
        abc = self.lattice.lat_lengths()
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
        struct = Atoms(lattice_mat=lat, elements=typ_sp, coords=coords, cartesian=False)
        return struct

    def pymatgen_converter(self):
        """
        Get pymatgen representation of the atoms object
        """
        try:
            from pymatgen.core.structure import Structure

            return Structure(
                self.lattice_mat,
                self.elements,
                self.frac_coords,
                coords_are_cartesian=False,
            )
        except:
            pass

    def spacegroup(self, symprec=1e-3):
        """
        Get spacegroup of the atoms object
        """

        import spglib

        sg = spglib.get_spacegroup(
            (self.lattice_mat, self.frac_coords, self.atomic_numbers), symprec=symprec
        )
        return sg

    @property
    def packing_fraction(self):
        """
        Get packing fraction of the atoms object
        """
        total_rad = 0
        for i in self.elements:
            total_rad = total_rad + Specie(i).atomic_rad ** 3
        pf = np.array([4 * np.pi * total_rad / (3 * self.volume)])
        return round(pf[0], 5)

    def lattice_points_in_supercell(self, supercell_matrix):
        """
        Adapted from Pymatgen
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

        ar = np.arange(mins[0], maxes[0])[:, None] * np.array([1, 0, 0])[None, :]
        br = np.arange(mins[1], maxes[1])[:, None] * np.array([0, 1, 0])[None, :]
        cr = np.arange(mins[2], maxes[2])[:, None] * np.array([0, 0, 1])[None, :]

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
        Adapted from Pymatgen
        Makes a supercell. Allowing to have sites outside the unit cell

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
        """
        Make supercell of dimension dim
        """
        dim = np.array(dim)
        if dim.shape == (3, 3):
            dim = np.array([int(np.linalg.norm(v)) for v in dim])
        coords = self.frac_coords
        all_symbs = self.elements  # [i.symbol for i in s.species]
        nat = len(coords)

        new_nat = nat * dim[0] * dim[1] * dim[2]
        # print ('new_nat,dim',new_nat,dim)
        new_coords = np.zeros((new_nat, 3))
        new_symbs = []  # np.chararray((new_nat))
        props = []  # self.props

        count = 0
        for i in range(nat):
            for j in range(dim[0]):
                for k in range(dim[1]):
                    for l in range(dim[2]):
                        props.append(self.props[i])
                        new_coords[count][0] = (coords[i][0] + j) / float(dim[0])
                        new_coords[count][1] = (coords[i][1] + k) / float(dim[1])
                        new_coords[count][2] = (coords[i][2] + l) / float(dim[2])
                        new_symbs.append(all_symbs[i])
                        count = count + 1

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

    def get_string(self):
        """
        Get string representation of the atoms object
        """
        system = str(self.composition.reduced_formula)
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
        order = np.argsort(self.elements)
        coords = self.frac_coords

        coords_ordered = np.array(coords)[order]
        elements_ordered = np.array(self.elements)[order]
        props_ordered = np.array(self.props)[order]
        counts = get_counts(elements_ordered)

        middle = (
            " ".join(map(str, counts.keys()))
            + "\n"
            + " ".join(map(str, counts.values()))
            + "\ndirect\n"
        )
        rest = ""
        for ii, i in enumerate(coords_ordered):

            if self.show_props == True:
                rest = (
                    rest + " ".join(map(str, i)) + " " + str(props_ordered[ii]) + "\n"
                )
            else:
                rest = rest + " ".join(map(str, i)) + "\n"

        result = header + middle + rest

        return result

    def get_lll_reduced_structure(self):
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
        system = str(self.composition.reduced_formula)
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
        order = np.argsort(self.elements)
        coords = self.frac_coords
        coords_ordered = np.array(coords)#[order]
        elements_ordered = np.array(self.elements)#[order]
        props_ordered = np.array(self.props)#[order]
        check_selective_dynamics = False
        counts = get_counts(elements_ordered)
        if "T" in "".join(map(str, self.props[0])):
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

        for ii, i in enumerate(coords_ordered):
            if self.show_props == True:
                rest = (
                    rest + " ".join(map(str, i)) + " " + str(props_ordered[ii]) + "\n"
                )
            else:
                rest = rest + " ".join(map(str, i)) + "\n"

        result = header + middle + rest
        return result




class VacuumPadding(object):
    """
    Adds vaccum padding to make 2D structure or making molecules
    """

    def __init__(self, atoms, vacuum=20.0):
        self.atoms = atoms
        self.vacuum = vacuum

    def get_effective_2d_slab(self):
        """
        Adds 2D vacuum to a system
        """
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

        lattice_mat = self.atoms.lattice_mat
        # lattice_mat[2][2] = padding
        # elements = self.atoms.elements
        # atoms = Atoms(lattice_mat = lattice_mat, coords = self.atoms.cart_coords, elements = elements, cartesian = True)
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
                / norm(np.cross(a1, a2)) ** 2,
            ]
        )

        a1 = new_lat[0]
        a2 = new_lat[1]
        a3 = new_lat[2]
        # print("a1,a2,a3", new_lat)

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
            lattice_mat=new_lat, elements=elements, coords=frac, cartesian=False
        )
        return with_vacuum_atoms

    def get_effective_molecule(self):
        """
        Adds vacuum around a system
        """
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
            lattice_mat=lattice_mat, elements=elements, coords=frac, cartesian=False
        )
        return with_vacuum_atoms




"""
if __name__ == "__main__":
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
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
# if __name__=='__main__':
#    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
#    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
#    elements = ["Si", "Si"]
#    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
#    comp=Si.composition
#    #print (comp,Si.density)
#    print (Si.atomic_numbers)
#   print (Si.pymatgen_converter().composition.weight,Si.composition.weight,Si.density)
"""
