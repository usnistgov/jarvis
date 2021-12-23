"""Modules for handling crystallographic Spacegroup related operations."""
from functools import reduce
from jarvis.core.lattice import Lattice
from jarvis.core.atoms import Atoms
import spglib
from jarvis.core.specie import Specie
import numpy as np
from numpy import sin, cos
import itertools

# from fractions import gcd
from jarvis.core.utils import gcd

# from numpy import gcd
# from math import gcd
import os


def unique_rows_2(a):
    """Remove duplicate rows."""
    order = np.lexsort(a.T)
    a = a[order]
    diff = np.diff(a, axis=0)
    ui = np.ones(len(a), "bool")
    ui[1:] = (diff != 0).any(axis=1)
    return a[ui]


def symmetrically_distinct_miller_indices(max_index=3, cvn_atoms=None):
    """Get unique miller indices for max_index."""
    # Need to work on this
    r1 = list(range(1, max_index + 1))
    r2 = list(range(-max_index, 1))
    r2.reverse()
    r3 = r1 + r2
    # r.reverse()
    r = r3
    # print ('sorted',r,sorted(r))
    conv_hkl_list = [
        miller
        for miller in itertools.product(r, r, r)
        if any([i != 0 for i in miller])
    ]
    spg = Spacegroup3D(cvn_atoms)._dataset
    rot = spg["rotations"]
    done = []
    unique_millers = []
    # print (conv_hkl_list)
    # print (sorted(conv_hkl_list, key=lambda x: x[0], reverse=True))
    for i in conv_hkl_list:
        d = abs(reduce(gcd, i))
        miller = tuple([int(k / d) for k in i])
        for j in rot:
            prod = list(np.dot(miller, j))
            if prod not in done:
                unique_millers.append(i)
                done.append(prod)
    uniq = unique_rows_2(np.array(unique_millers))
    return uniq


wyckoff_file = str(os.path.join(os.path.dirname(__file__), "Wyckoff.csv"))


def parse_wyckoff_csv(wyckoff_file):
    """Parse Wyckoff.csv from spglib.

    There are 530 data sets. For one example:
    9:C 1 2 1:::::::
    ::4:c:1:(x,y,z):(-x,y,-z)::
    ::2:b:2:(0,y,1/2):::
    ::2:a:2:(0,y,0):::
    """
    rowdata = []
    points = []
    hP_nums = [433, 436, 444, 450, 452, 458, 460]
    for i, line in enumerate(wyckoff_file):
        if line.strip() == "end of data":
            break
        rowdata.append(line.strip().split(":"))

        # 2:P -1  ::::::: <-- store line number if first element is number
        if rowdata[-1][0].isdigit():
            points.append(i)
    points.append(i)

    wyckoff = []
    for i in range(len(points) - 1):  # 0 to 529
        symbol = rowdata[points[i]][1]  # e.g. "C 1 2 1"
        if i + 1 in hP_nums:
            symbol = symbol.replace("R", "H", 1)
        wyckoff.append({"symbol": symbol.strip()})

    # When the number of positions is larger than 4,
    # the positions are written in the next line.
    # So those positions are connected.
    for i in range(len(points) - 1):
        count = 0
        wyckoff[i]["wyckoff"] = []
        for j in range(points[i] + 1, points[i + 1]):
            # Hook if the third element is a number (multiplicity), e.g.,
            #
            # 232:P 2/b 2/m 2/b:::::::  <- ignored
            # ::8:r:1:(x,y,z):(-x,y,-z):(x,-y+1/2,-z):(-x,-y+1/2,z)
            # :::::(-x,-y,-z):(x,-y,z):(-x,y+1/2,z):(x,y+1/2,-z)  <- ignored
            # ::4:q:..m:(x,0,z):(-x,0,-z):(x,1/2,-z):(-x,1/2,z)
            # ::4:p:..2:(0,y,1/2):(0,-y+1/2,1/2):(0,-y,1/2):(0,y+1/2,1/2)
            # ::4:o:..2:(1/2,y,0):(1/2,-y+1/2,0):(1/2,-y,0):(1/2,y+1/2,0)
            # ...
            if rowdata[j][2].isdigit():
                pos = []
                w = {
                    "letter": rowdata[j][3].strip(),
                    "multiplicity": int(rowdata[j][2]),
                    "site_symmetry": rowdata[j][4].strip(),
                    "positions": pos,
                }
                wyckoff[i]["wyckoff"].append(w)

                for k in range(4):
                    if rowdata[j][k + 5]:  # check if '(x,y,z)' or ''
                        count += 1
                        pos.append(rowdata[j][k + 5])
            else:
                for k in range(4):
                    if rowdata[j][k + 5]:
                        count += 1
                        pos.append(rowdata[j][k + 5])

        # assertion
        # for w in wyckoff[i]['wyckoff']:
        #    n_pos = len(w['positions'])
        #    n_pos *= len(lattice_symbols[wyckoff[i]['symbol'][0]])
        #    assert n_pos == w['multiplicity']

    return wyckoff


def read_wyckoff_csv(filename):
    """Read wyckoff_csv file."""
    with open(filename) as wyckoff_file:
        return parse_wyckoff_csv(wyckoff_file)


def get_wyckoff_position_operators(hall_number):
    """Get all Wyckoff operations for Hall number."""
    wyckoff = read_wyckoff_csv(wyckoff_file)
    operations = wyckoff[hall_number - 1]
    return operations


class Spacegroup3D(object):
    """
    Provide spacegroup related data for Atoms object.

    Currently uses spglib to derive spacegroup
    related information for 3D materials mainly
    """

    def __init__(self, atoms=[], dataset={}, symprec=1e-2, angle_tolerance=5):
        """
        Following information are needed for Spacegroup3D.

        If dataset is not provided, the default dataset is used.
        Args:
            atoms: jarvis.core.Atoms

            dataset: spacegroup dataset

            symprec: symmetry precision

            angle_tolerance: angle tolerance
        """
        self._dataset = dataset
        self._atoms = atoms
        self._symprec = symprec
        self._angle_tolerance = angle_tolerance
        if self._dataset == {}:
            spg = self.spacegroup_data()
            self._dataset = spg._dataset

    def spacegroup_data(self):
        """Provide spacegroup data from spglib."""
        phonopy_atoms = (
            self._atoms.lattice_mat,
            self._atoms.frac_coords,
            self._atoms.atomic_numbers,
        )
        dataset = spglib.get_symmetry_dataset(
            phonopy_atoms,
            symprec=self._symprec,
            angle_tolerance=self._angle_tolerance,
        )
        """
        keys = ('number',
        'hall_number',
        'international',
        'hall',
        'choice',
        'transformation_matrix',
        'origin_shift',
        'rotations',
        'translations',
        'wyckoffs',
        'site_symmetry_symbols',
        'equivalent_atoms',
        'mapping_to_primitive',
        'std_lattice',
        'std_types',
        'std_positions',
        'std_rotation_matrix',
        'std_mapping_to_primitive',
        'pointgroup')
        """
        return Spacegroup3D(
            atoms=self._atoms,
            symprec=self._symprec,
            angle_tolerance=self._angle_tolerance,
            dataset=dataset,
        )

    @property
    def space_group_symbol(self):
        """Get spacegroup symbol."""
        # spg = self.spacegroup_data()
        return self._dataset["international"]

    @property
    def space_group_number(self):
        """Get spacegroup number."""
        # spg = self.spacegroup_data()
        return self._dataset["number"]

    @property
    def primitive_atoms(self):
        """Get primitive atoms."""
        phonopy_atoms = (
            self._atoms.lattice_mat,
            self._atoms.frac_coords,
            self._atoms.atomic_numbers,
        )
        lattice, scaled_positions, numbers = spglib.find_primitive(
            phonopy_atoms, symprec=self._symprec
        )
        elements = self._atoms.elements
        el_dict = {}
        for i in elements:
            el_dict.setdefault(Specie(i).Z, i)
        prim_elements = [el_dict[i] for i in numbers]
        prim_atoms = Atoms(
            lattice_mat=lattice,
            elements=prim_elements,
            coords=scaled_positions,
            cartesian=False,
        )
        return prim_atoms

    @property
    def refined_atoms(self):
        """Refine atoms based on spacegroup data."""
        phonopy_atoms = (
            self._atoms.lattice_mat,
            self._atoms.frac_coords,
            self._atoms.atomic_numbers,
        )
        lattice, scaled_positions, numbers = spglib.refine_cell(
            phonopy_atoms, self._symprec, self._angle_tolerance
        )
        elements = self._atoms.elements
        el_dict = {}
        for i in elements:
            el_dict.setdefault(Specie(i).Z, i)
        ref_elements = [el_dict[i] for i in numbers]
        ref_atoms = Atoms(
            lattice_mat=lattice,
            elements=ref_elements,
            coords=scaled_positions,
            cartesian=False,
        )
        return ref_atoms

    @property
    def crystal_system(self):
        """Get crystal system."""
        n = self._dataset["number"]

        def f(i, j):
            return i <= n <= j

        cs = {
            "triclinic": (1, 2),
            "monoclinic": (3, 15),
            "orthorhombic": (16, 74),
            "tetragonal": (75, 142),
            "trigonal": (143, 167),
            "hexagonal": (168, 194),
            "cubic": (195, 230),
        }

        crystal_sytem = None

        for k, v in cs.items():
            if f(*v):
                crystal_sytem = k
                break
        return crystal_sytem

    @property
    def lattice_system(self):
        """Get lattice system."""
        n = self._dataset["number"]
        system = self.crystal_system
        if n in [146, 148, 155, 160, 161, 166, 167]:
            return "rhombohedral"
        elif system == "trigonal":
            return "hexagonal"
        else:
            return system

    @property
    def point_group_symbol(self):
        """Get pointgroup."""
        return self._dataset["pointgroup"]

    @property
    def conventional_standard_structure(
        self, tol=1e-5, international_monoclinic=True
    ):
        """
        Give a conventional cell according to certain conventions.

        The conventionss are defined in Setyawan, W., & Curtarolo,
        S. (2010). High-throughput electronic band structure calculations:
        Challenges and tools. Computational Materials Science,
        49(2), 299-312. doi:10.1016/j.commatsci.2010.05.010
        They basically enforce as much as possible
        norm(a1)<norm(a2)<norm(a3)
        Returns:
            The structure in a conventional standardized cell
        """
        struct = self.refined_atoms
        latt = struct.lattice
        latt_type = self.lattice_system
        sorted_lengths = sorted(latt.abc)
        sorted_dic = sorted(
            [
                {"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i}
                for i in [0, 1, 2]
            ],
            key=lambda k: k["length"],
        )

        if latt_type in ("orthorhombic", "cubic"):
            # you want to keep the c axis where it is
            # to keep the C- settings
            transf = np.zeros(shape=(3, 3))
            if self.space_group_symbol.startswith("C"):
                transf[2] = [0, 0, 1]
                a, b = sorted(latt.abc[:2])
                sorted_dic = sorted(
                    [
                        {
                            "vec": latt.matrix[i],
                            "length": latt.abc[i],
                            "orig_index": i,
                        }
                        for i in [0, 1]
                    ],
                    key=lambda k: k["length"],
                )
                for i in range(2):
                    transf[i][sorted_dic[i]["orig_index"]] = 1
                c = latt.abc[2]
            elif self.space_group_symbol.startswith(
                "A"
            ):  # change to C-centering to match Setyawan/Curtarolo convention
                transf[2] = [1, 0, 0]
                a, b = sorted(latt.abc[1:])
                sorted_dic = sorted(
                    [
                        {
                            "vec": latt.matrix[i],
                            "length": latt.abc[i],
                            "orig_index": i,
                        }
                        for i in [1, 2]
                    ],
                    key=lambda k: k["length"],
                )
                for i in range(2):
                    transf[i][sorted_dic[i]["orig_index"]] = 1
                c = latt.abc[0]
            else:
                for i in range(len(sorted_dic)):
                    transf[i][sorted_dic[i]["orig_index"]] = 1
                a, b, c = sorted_lengths
            latt = Lattice.orthorhombic(a, b, c)

        elif latt_type == "tetragonal":
            # find the "a" vectors
            # it is basically the vector repeated two times
            transf = np.zeros(shape=(3, 3))
            a, b, c = sorted_lengths
            for d in range(len(sorted_dic)):
                transf[d][sorted_dic[d]["orig_index"]] = 1

            if abs(b - c) < tol and abs(a - c) > tol:
                a, c = c, a
                transf = np.dot([[0, 0, 1], [0, 1, 0], [1, 0, 0]], transf)
            latt = Lattice.tetragonal(a, c)
        elif latt_type in ("hexagonal", "rhombohedral"):
            # for the conventional cell representation,
            # we allways show the rhombohedral lattices as hexagonal

            # check first if we have the refined structure shows a rhombohedral
            # cell
            # if so, make a supercell
            a, b, c = latt.abc
            if np.all(np.abs([a - b, c - b, a - c]) < 0.001):
                struct.make_supercell(((1, -1, 0), (0, 1, -1), (1, 1, 1)))
                a, b, c = sorted(struct.lattice.abc)

            if abs(b - c) < 0.001:
                a, c = c, a
            new_matrix = [
                [a / 2, -a * np.sqrt(3) / 2, 0],
                [a / 2, a * np.sqrt(3) / 2, 0],
                [0, 0, c],
            ]
            latt = Lattice(new_matrix)
            transf = np.eye(3, 3)

        elif latt_type == "monoclinic":
            # You want to keep the c axis where it is to keep the C- settings

            if self.space_group_symbol.startswith("C"):
                transf = np.zeros(shape=(3, 3))
                transf[2] = [0, 0, 1]
                sorted_dic = sorted(
                    [
                        {
                            "vec": latt.matrix[i],
                            "length": latt.abc[i],
                            "orig_index": i,
                        }
                        for i in [0, 1]
                    ],
                    key=lambda k: k["length"],
                )
                a = sorted_dic[0]["length"]
                b = sorted_dic[1]["length"]
                c = latt.abc[2]
                new_matrix = None
                for t in itertools.permutations(list(range(2)), 2):
                    m = latt.matrix
                    latt2 = Lattice([m[t[0]], m[t[1]], m[2]])
                    lengths = latt2.abc
                    angles = latt2.angles
                    if angles[0] > 90:
                        # if the angle is > 90 we invert a and b to get
                        # an angle < 90
                        a, b, c, alpha, beta, gamma = Lattice(
                            [-m[t[0]], -m[t[1]], m[2]]
                        ).parameters
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][2] = 1
                        alpha = np.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                        continue

                    elif angles[0] < 90:
                        transf = np.zeros(shape=(3, 3))
                        # print ('464-470')
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][2] = 1
                        a, b, c = lengths
                        alpha = np.pi * angles[0] / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]

                if new_matrix is None:
                    # print ('479-482')
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [[a, 0, 0], [0, b, 0], [0, 0, c]]
                    transf = np.zeros(shape=(3, 3))
                    for c in range(len(sorted_dic)):
                        transf[c][sorted_dic[c]["orig_index"]] = 1
            # if not C-setting
            else:
                # try all permutations of the axis
                # keep the ones with the non-90 angle=alpha
                # and b<c
                new_matrix = None
                for t in itertools.permutations(list(range(3)), 3):
                    m = latt.matrix
                    a, b, c, alpha, beta, gamma = Lattice(
                        [m[t[0]], m[t[1]], m[t[2]]]
                    ).parameters
                    if alpha > 90 and b < c:
                        a, b, c, alpha, beta, gamma = Lattice(
                            [-m[t[0]], -m[t[1]], m[t[2]]]
                        ).parameters
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = -1
                        transf[1][t[1]] = -1
                        transf[2][t[2]] = 1
                        alpha = np.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                        continue
                    elif alpha < 90 and b < c:
                        # print ('510-515')
                        transf = np.zeros(shape=(3, 3))
                        transf[0][t[0]] = 1
                        transf[1][t[1]] = 1
                        transf[2][t[2]] = 1
                        alpha = np.pi * alpha / 180
                        new_matrix = [
                            [a, 0, 0],
                            [0, b, 0],
                            [0, c * cos(alpha), c * sin(alpha)],
                        ]
                if new_matrix is None:
                    # print ('523-530')
                    # this if is to treat the case
                    # where alpha==90 (but we still have a monoclinic sg
                    new_matrix = [
                        [sorted_lengths[0], 0, 0],
                        [0, sorted_lengths[1], 0],
                        [0, 0, sorted_lengths[2]],
                    ]
                    transf = np.zeros(shape=(3, 3))
                    for c in range(len(sorted_dic)):
                        transf[c][sorted_dic[c]["orig_index"]] = 1

            if international_monoclinic:
                # The above code makes alpha the non-right angle.
                # The following will convert to proper international convention
                # that beta is the non-right angle.
                op = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
                transf = np.dot(op, transf)
                new_matrix = np.dot(op, new_matrix)
                beta = Lattice(new_matrix).beta
                if beta < 90:
                    op = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                    transf = np.dot(op, transf)
                    new_matrix = np.dot(op, new_matrix)

            latt = Lattice(new_matrix)

        elif latt_type == "triclinic":
            # we use a LLL Minkowski-like reduction for the triclinic cells
            struct = struct.get_lll_reduced_structure()

            a, b, c = latt.abc  # lengths
            alpha, beta, gamma = [np.pi * i / 180 for i in latt.angles]
            new_matrix = None
            test_matrix = [
                [a, 0, 0],
                [b * cos(gamma), b * sin(gamma), 0.0],
                [
                    c * cos(beta),
                    c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    c
                    * np.sqrt(
                        sin(gamma) ** 2
                        - cos(alpha) ** 2
                        - cos(beta) ** 2
                        + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]

            def is_all_acute_or_obtuse(m):
                recp_angles = np.array(Lattice(m).reciprocal_lattice().angles)
                return np.all(recp_angles <= 90) or np.all(recp_angles > 90)

            if is_all_acute_or_obtuse(test_matrix):
                transf = np.eye(3)
                new_matrix = test_matrix

            test_matrix = [
                [-a, 0, 0],
                [b * cos(gamma), b * sin(gamma), 0.0],
                [
                    -c * cos(beta),
                    -c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    -c
                    * np.sqrt(
                        sin(gamma) ** 2
                        - cos(alpha) ** 2
                        - cos(beta) ** 2
                        + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
                new_matrix = test_matrix

            test_matrix = [
                [-a, 0, 0],
                [-b * cos(gamma), -b * sin(gamma), 0.0],
                [
                    c * cos(beta),
                    c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    c
                    * np.sqrt(
                        sin(gamma) ** 2
                        - cos(alpha) ** 2
                        - cos(beta) ** 2
                        + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]

            if is_all_acute_or_obtuse(test_matrix):
                transf = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
                new_matrix = test_matrix

            test_matrix = [
                [a, 0, 0],
                [-b * cos(gamma), -b * sin(gamma), 0.0],
                [
                    -c * cos(beta),
                    -c * (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma),
                    -c
                    * np.sqrt(
                        sin(gamma) ** 2
                        - cos(alpha) ** 2
                        - cos(beta) ** 2
                        + 2 * cos(alpha) * cos(beta) * cos(gamma)
                    )
                    / sin(gamma),
                ],
            ]
            if is_all_acute_or_obtuse(test_matrix):
                transf = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
                new_matrix = test_matrix

            latt = Lattice(new_matrix)

        new_coords = np.dot(transf, np.transpose(struct.frac_coords)).T
        new_struct = Atoms(
            lattice_mat=latt.matrix,
            elements=struct.elements,
            coords=new_coords,
            cartesian=False,
        )
        return new_struct


def parse_xyz_string(xyz_string):
    """
    Convert xyz info to translation and rotation vectors.

    Adapted from pymatgen.
    Args:
        xyz_string: string of the form 'x, y, z', '-x, -y, z',
            '-2y+1/2, 3x+1/2, z-y+1/2', etc.
    Returns:
        translation and rotation vectors.
    """
    from jarvis.core.utils import parse_xyz_string

    return parse_xyz_string(xyz_string)


def operate_affine(cart_coord=[], affine_matrix=[]):
    """Operate affine method."""
    affine_point = np.array([cart_coord[0], cart_coord[1], cart_coord[2], 1])
    return np.dot(np.array(affine_matrix), affine_point)[0:3]


def get_new_coord_for_xyz_sym(frac_coord=[], xyz_string=""):
    """Obtain new coord from xyz string."""
    from jarvis.core.utils import get_new_coord_for_xyz_sym

    return get_new_coord_for_xyz_sym(
        frac_coord=frac_coord, xyz_string=xyz_string
    )


def check_duplicate_coords(coords=[], coord=[]):
    """Check if a coordinate exists."""
    from jarvis.core.utils import check_duplicate_coords

    return check_duplicate_coords(coords=coords, coord=coord)


"""
if __name__ == "__main__":
    x = get_wyckoff_position_operators(488)
    print (x)
    import sys
    #sys.exit()
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(atoms=Si)  # .spacegroup_data()
    #print(spg.space_group_symbol)
    #print(spg.space_group_number)
    #primt = spg.primitive_atoms
    #print("primt", primt)
    #print("cryst_sys", spg.crystal_system)
    #print("point group", spg.point_group_symbol)
    #print("cvn", spg.conventional_standard_structure)
    cvn = spg.conventional_standard_structure
    ml=symmetrically_distinct_miller_indices(max_index=3, cvn_atoms=cvn)
    print ('miller=',ml,len(ml))
#    ml_pmg=get_symmetrically_distinct_miller_indices(cvn.pymatgen_converter(),max_index=3)
#    print ('millerpmg=',ml_pmg,len(ml_pmg))
    # pmg=Si.pymatgen_converter()
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # spg=SpacegroupAnalyzer(pmg)
    # print (spg.get_space_group_symbol(),spg.get_space_group_number())
    # print (pmg.get_primitive_structure())
"""
