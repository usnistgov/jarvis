"""Modules for making crystallographic plane surfaces."""
from jarvis.core.atoms import Atoms
from jarvis.core.utils import ext_gcd
import numpy as np
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from numpy.linalg import norm
from numpy import gcd
from collections import OrderedDict


def wulff_normals(miller_indices=[], surface_energies=[]):
    """Obtain Wulff Normals.

    Args:
         miller_indices : Miller indices

         surface_energies : corresponding surface energies

    Returns: Surface normals
    """
    max_s = min(surface_energies)
    normals = []
    for i, j in zip(miller_indices, surface_energies):
        normal = j * np.linalg.norm(i) / float(max_s)
        normals.append([normal, i])

    return normals


class Surface(object):
    """Get surface object of arbitrary atoms object and miller index."""

    def __init__(
        self,
        atoms=None,
        indices=[0, 0, 1],
        layers=3,
        thickness=None,
        vacuum=18.0,
        tol=1e-10,
        from_conventional_structure=True,
        use_thickness_c=True,
    ):
        """Initialize the class.

        Args:
             atoms: jarvis.core.Atoms object

             indices: Miller indices

             layers: Number of surface layers

             thickness: Provide thickness instead of layers

             vacuum: vacuum padding

             tol: tolerance during dot product

             from_conventional_structure: whether to use the conv. atoms
        """
        self.indices = np.array(indices)
        self.from_conventional_structure = from_conventional_structure
        if self.from_conventional_structure:
            self.atoms = Spacegroup3D(atoms).conventional_standard_structure
        else:
            self.atoms = atoms
        self.tol = tol
        self.vacuum = vacuum
        self.layers = layers
        self.thickness = thickness
        self.use_thickness_c = use_thickness_c
        # Note thickness overwrites layers

    def to_dict(self):
        """Convert to a dictionary."""
        d = OrderedDict()
        d["atoms"] = self.atoms.to_dict()
        d["indices"] = self.indices
        d["tol"] = self.tol
        d["vacuum"] = self.vacuum
        d["layers"] = self.layers
        d["from_conventional_structure"] = self.from_conventional_structure
        return d

    @classmethod
    def from_dict(self, d={}):
        """Construct class from a dictionary."""
        return Surface(
            atoms=Atoms.from_dict(d["atoms"]),
            indices=d["indices"],
            tol=d["tol"],
            vacuum=d["vacuum"],
            layers=d["layers"],
            from_conventional_structure=d["from_conventional_structure"],
        )

    def make_surface(self):
        """Generate specified surface. Modified from ase package."""
        atoms = self.atoms
        h_index, k_index, l_index = self.indices
        h0, k0, l0 = self.indices == 0
        if h0 and k0 or h0 and l0 or k0 and l0:  # if two indices are zero
            if not h0:
                c1, c2, c3 = [(0, 1, 0), (0, 0, 1), (1, 0, 0)]
            if not k0:
                c1, c2, c3 = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]
            if not l0:
                c1, c2, c3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        else:
            p, q = ext_gcd(k_index, l_index)
            a1, a2, a3 = self.atoms.lattice_mat  # .lat_lengths()

            # constants describing the dot product of basis c1 and c2:
            # dot(c1,c2) = k1+i*k2, i in Z
            k1 = np.dot(
                p * (k_index * a1 - h_index * a2)
                + q * (l_index * a1 - h_index * a3),
                l_index * a2 - k_index * a3,
            )
            k2 = np.dot(
                l_index * (k_index * a1 - h_index * a2)
                - k_index * (l_index * a1 - h_index * a3),
                l_index * a2 - k_index * a3,
            )

            if abs(k2) > self.tol:
                i = -int(round(k1 / k2))
                p, q = p + i * l_index, q - i * k_index

            a, b = ext_gcd(p * k_index + q * l_index, h_index)

            c1 = (p * k_index + q * l_index, -p * h_index, -q * h_index)
            c2 = np.array((0, l_index, -k_index)) // abs(gcd(l_index, k_index))
            c3 = (b, a * p, a * q)
        lattice = atoms.lattice_mat  # .lat_lengths()
        basis = np.array([c1, c2, c3])
        scaled = np.linalg.solve(basis.T, np.array(atoms.frac_coords).T).T
        scaled -= np.floor(scaled + self.tol)
        new_coords = scaled
        tmp_cell = np.dot(basis, lattice)
        M = np.linalg.solve(lattice, tmp_cell)
        cart_coords = np.dot(scaled, lattice)
        new_coords = np.dot(cart_coords, M)

        new_atoms = Atoms(
            lattice_mat=tmp_cell,
            coords=new_coords,
            elements=atoms.elements,
            cartesian=True,
        )
        if self.thickness is not None and (self.thickness) > 0:
            if not self.use_thickness_c:
                new_lat = new_atoms.lattice_mat  # lat_lengths()
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

                a3 = new_lat[2]
                self.layers = int(self.thickness / np.linalg.norm(a3)) + 1
            else:
                self.layers = int(self.thickness / new_atoms.lattice.c) + 1
            # print("self.layers", self.layers)
            # dims=get_supercell_dims(new_atoms,enforce_c_size=self.thickness)
            # print ('dims=',dims,self.layers)
            # surf_atoms = new_atoms.make_supercell_matrix([1, 1, dims[2]])
            # print('self.layers',self.layers,self.thickness,new_atoms.lattice.c)
        surf_atoms = new_atoms.make_supercell_matrix([1, 1, self.layers])
        # print("supercell_cart_coords", surf_atoms.frac_coords)

        new_lat = surf_atoms.lattice_mat  # lat_lengths()
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

        M = np.linalg.solve(new_lat, latest_lat)

        new_cart_coords = surf_atoms.cart_coords  # np.dot(scaled,lattice)

        new_coords = np.dot(new_cart_coords, M)

        new_atoms = Atoms(
            lattice_mat=latest_lat,
            elements=surf_atoms.elements,
            coords=new_coords,
            cartesian=True,
        ).center_around_origin()

        frac_coords = new_atoms.frac_coords

        # frac_coords[:] = frac_coords[:] % 1
        new_atoms = Atoms(
            lattice_mat=latest_lat,
            elements=surf_atoms.elements,
            coords=frac_coords,
            cartesian=False,
        )
        new_lat = new_atoms.lattice_mat
        new_cart_coords = new_atoms.cart_coords
        elements = new_atoms.elements
        new_lat[2][2] = new_lat[2][2] + self.vacuum
        with_vacuum_atoms = Atoms(
            lattice_mat=new_lat,
            elements=elements,
            coords=new_cart_coords,
            cartesian=True,
        )
        # new_atoms.center()
        # print (with_vacuum_atoms)
        return with_vacuum_atoms


"""
if __name__ == "__main__":
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Surface(atoms=Si, indices=[1, 1, 1]).make_surface()
    su = [
        0.8582640971273426,
        0.9334963319196496,
        0.9360461382184894,
        0.9419095687284446,
        0.9802042233627004,
        0.9875446840480956,
        1.0120634294466684,
        1.0126231880823566,
        1.0241538763302507,
        1.0315901848682645,
        1.0318271257831195,
        1.0331286888257398,
        1.0344297141291043,
        1.0388709097092674,
        1.040277640596931,
        1.042494119906149,
        1.04453679643896,
        1.0450598648770613,
        1.045076130339553,
        1.0469310544190567,
        1.0491015867538047,
        1.0495494553198788,
        1.0534717916897114,
        1.0535201391639715,
        1.054233162444997,
        1.0579157863887743,
        1.0595676718662346,
        1.0601381085497692,
        1.109580394178689,
    ]

    ml = [
        [0, 0, 1],
        [2, 0, 3],
        [2, 0, 1],
        [1, 0, 1],
        [3, 0, 2],
        [1, 0, 3],
        [3, 1, 1],
        [3, 0, 1],
        [3, 1, 3],
        [3, -1, 1],
        [3, 1, 0],
        [3, 2, 1],
        [3, 3, 1],
        [1, 0, 0],
        [2, 2, 1],
        [3, -1, 3],
        [3, -1, 2],
        [3, 3, 2],
        [3, 2, 2],
        [2, -1, 3],
        [3, 2, 0],
        [3, 2, 3],
        [1, 1, 1],
        [1, 0, 2],
        [3, 1, 2],
        [2, -1, 2],
        [3, -1, 0],
        [2, 2, 3],
        [1, 1, 0],
    ]
    nm = wulff_normals(miller_indices=ml, surface_energies=su)
    print(nm)
    from jarvis.core.lattice import Lattice
    lat = Lattice([[4.05, 0, 0], [0, 4.05, 0], [0, 0, 4.05]])
    pmg_wulff = WulffShape(lat, ml, su)
    print(pmg_wulff.facets)
"""
