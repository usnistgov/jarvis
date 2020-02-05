from jarvis.core.lattice import Lattice
from jarvis.core.atoms import Atoms
import spglib
from jarvis.core.specie import Specie
import itertools
import numpy as np
from numpy import sin, cos
import itertools

def symmetrically_distinct_miller_indices(max_index=3,cvn_atoms=None):
    r = list(range(-max_index, max_index + 1))
    r.reverse()
    conv_hkl_list = [miller for miller in itertools.product(r, r, r) if any([i != 0 for i in miller])]
    return conv_hkl_list

class Spacegroup3D(object):
    def __init__(self, atoms=[], dataset={}, symprec=1e-2, angle_tolerance=5):
        self._dataset = dataset
        self._atoms = atoms
        self._symprec = symprec
        self._angle_tolerance = angle_tolerance
        if self._dataset == {}:
            spg = self.spacegroup_data()
            self._dataset = spg._dataset

    def spacegroup_data(self):
        phonopy_atoms = (
            self._atoms.lattice_mat,
            self._atoms.frac_coords,
            self._atoms.atomic_numbers,
        )
        dataset = spglib.get_symmetry_dataset(
            phonopy_atoms, symprec=self._symprec, angle_tolerance=self._angle_tolerance
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
        # spg = self.spacegroup_data()
        return self._dataset["international"]

    @property
    def space_group_number(self):
        # spg = self.spacegroup_data()
        return self._dataset["number"]

    @property
    def primitive_atoms(self):
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
        return self._dataset["pointgroup"]

    @property
    def conventional_standard_structure(self, tol=1e-5, international_monoclinic=True):
        """
        Gives a structure with a conventional cell according to certain
        standards. The standards are defined in Setyawan, W., & Curtarolo,
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
                        {"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i}
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
                        {"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i}
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
                        {"vec": latt.matrix[i], "length": latt.abc[i], "orig_index": i}
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
                        # print ('[-m[t[0]], -m[t[1]], m[2]]',[-m[t[0]], -m[t[1]], m[2]])
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
                recp_angles = np.array(Lattice(m).reciprocal_lattice.angles)
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


if __name__ == "__main__":
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    spg = Spacegroup3D(atoms=Si)  # .spacegroup_data()
    print(spg.space_group_symbol)
    print(spg.space_group_number)
    primt = spg.primitive_atoms
    print("primt", primt)
    print("cryst_sys", spg.crystal_system)
    print("point group", spg.point_group_symbol)
    print("cvn", spg.conventional_standard_structure)
    ml=symmetrically_distinct_miller_indices(max_index=1)
    print ('miller=',ml)
    # pmg=Si.pymatgen_converter()
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # spg=SpacegroupAnalyzer(pmg)
    # print (spg.get_space_group_symbol(),spg.get_space_group_number())
    # print (pmg.get_primitive_structure())
