"""
This module provides classes to specify atomic structure
"""
# from jarvis.analysis.structure.spacegroup import Spacegroup3D
from collections import Counter
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie
from jarvis.core.lattice import Lattice
import matplotlib.pyplot as plt
from collections import OrderedDict
import pprint
import math

plt.switch_backend("agg")
amu_gm = 1.66054e-24
ang_cm = 1e-8


class Atoms(object):
    def __init__(
        self, lattice_mat=None, coords=None, elements=None, props=None, cartesian=False
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
        self.lattice = Lattice(lattice_mat)
        self.coords = coords
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
            for site,element in zip(self.frac_coords, self.elements):
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



    def to_dict(self):
        d = OrderedDict()
        d["lattice_mat"] = self.lattice_mat.tolist()
        d["coords"] = self.coords.tolist()
        d["elements"] = self.elements
        d["abc"] = self.lattice.lat_lengths()
        d["angles"] = self.lattice.lat_angles()
        d["cartesian"] = self.cartesian
        d["props"] = self.props
        return d

    def from_dict(self, d={}):
        return Atoms(
            lattice_mat=d["lattice_mat"],
            elements=d["elements"],
            props=d["props"],
            coords=d["coords"],
            cartesian=d["cartesian"],
        )

    def remove_site_by_index(self, site=0):
        # print ('elements',self.elements.pop(site))
        # print ('coords',self.frac_coords.tolist().pop(site))
        # print ('props',self.props.pop(site))
        new_els = []
        new_coords = []
        new_props = []
        for ii, i in enumerate(self.frac_coords):
            if ii != site:
                #print(self.elements, len(self.elements), len(self.frac_coords))
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

    def center(self, axis=2, vacuum=18.0, about=None):
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
        m = self.lattice_mat
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol

    @property
    def composition(self):
        comp = {}
        for i in self.elements:
            comp[i] = comp.setdefault(i, 0) + 1
        return Composition(comp)

    @property
    def density(self):
        den = float(self.composition.weight * amu_gm) / (
            float(self.volume) * (ang_cm) ** 3
        )
        return den

    @property
    def atomic_numbers(self):
        numbers = []
        for i in self.elements:
            numbers.append(Specie(i).Z)
        return numbers

    @property
    def num_atoms(self):
        return len(self.coords)

    def get_center_of_mass(self):
        # atomic_mass
        m = []
        for i in self.elements:
            m.append(Specie(i).atomic_mass)
        m = np.array(m)
        com = np.dot(m, self.cart_coords) / m.sum()
        # com = np.linalg.solve(self.lattice_mat.T, com)
        return com

    def pymatgen_converter(self):
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
        # try:
        import spglib

        sg = spglib.get_spacegroup(
            (self.lattice_mat, self.frac_coords, self.atomic_numbers), symprec=symprec
        )
        return sg

    # except:
    #    pass
    @property
    def packing_fraction(self):
        total_rad = 0
        for i in self.elements:
            total_rad = total_rad + Specie(i).atomic_rad ** 3
        pf = np.array([4 * np.pi * total_rad / (3 * self.volume)])
        return round(pf[0], 5)

    # def write(self,filename='POSCAR;,format='poscar'):

    def make_supercell(self, dim=[2, 2, 2]):
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
        header = (
            str("\nSystem\n1.0\n")
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
        middle = (
            " ".join(map(str, Counter(self.elements).keys()))
            + "\n"
            + " ".join(map(str, Counter(self.elements).values()))
            + "\ndirect\n"
        )
        rest = ""
        for i in self.frac_coords:
            rest = rest + " ".join(map(str, i)) + "\n"
        result = header + middle + rest
        return result

    # def __repr__(self,indent=4):
    #     return pprint.pformat(self.to_dict(), indent=indent)

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
    #def __str__(self):
    #     return "(%s(%r))" %(self.__class__,self.__dict__)
    def __repr__(self):
        header = (
            str("\nSystem\n1.0\n")
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
        middle = (
            " ".join(map(str, Counter(self.elements).keys()))
            + "\n"
            + " ".join(map(str, Counter(self.elements).values()))
            + "\ndirect\n"
        )
        rest = ""
        # print ('repr',self.frac_coords, self.frac_coords.shape)
        for ii, i in enumerate(self.frac_coords):
            rest = rest + " ".join(map(str, i)) + " " + str(self.props[ii]) + "\n"
        result = header + middle + rest
        return result


if __name__ == "__main__":
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
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
