"""Module for k-points used n various calculations."""

from jarvis.core.lattice import Lattice
import numpy as np
from numpy import linalg as LA
from collections import OrderedDict
import pprint
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from numpy import cos, sin
from math import tan, pi
from math import ceil
import math


def generate_kgrid(grid=[5, 5, 5]):
    """Generate k-mesh of size grid."""
    t = []
    for i in range(grid[0]):
        for j in range(grid[1]):
            for k in range(grid[2]):
                t.append(
                    [
                        float(i) / (float(grid[0])),
                        float(j) / (float(grid[1])),
                        float(k) / (float(grid[2])),
                    ]
                )
    return t


def generate_kpath(kpath=[[0, 0, 0], [0, 0.5, 0.5]], num_k=5):
    """Generate k-path with distance num_k k-points between them."""
    K = []
    for i in range(len(kpath) - 1):
        dk = np.array(kpath[i + 1]) - np.array(kpath[i])
        for j in range(num_k):
            K.append(np.array(kpath[i]) + dk * (float(j) / float(num_k)))
    K.append(kpath[-1])
    return K


class Kpoints3D(object):
    """Handle k-points python object."""

    def __init__(
        self,
        kpoints=[[1, 1, 1]],
        labels=[],
        kpoints_weights=[],
        kpoint_mode="automatic",
        header="Gamma",
    ):
        """Several types of k-points in the mesh or high-symmetry BZ."""
        self._kpoints = kpoints
        self._labels = labels
        self._kpoint_mode = kpoint_mode
        self._header = header
        self._kp_weights = kpoints_weights

    def automatic_length_mesh(self, lattice_mat=[], length=20, header="Gamma"):
        """Length based automatic k-points."""
        inv_lat = Lattice(lattice_mat=lattice_mat).inv_lattice()
        b1 = LA.norm(np.array(inv_lat[0]))
        b2 = LA.norm(np.array(inv_lat[1]))
        b3 = LA.norm(np.array(inv_lat[2]))
        n1 = int(max(1, length * b1 + 0.5))
        n2 = int(max(1, length * b2 + 0.5))
        n3 = int(max(1, length * b3 + 0.5))
        return Kpoints3D(
            kpoints=[[n1, n2, n3]], header=header, kpoint_mode="automatic"
        )

    @property
    def kpts(self):
        """Return k-points arrays."""
        return self._kpoints

    def kpoints_per_atom(self, atoms=None, kppa=1000):
        """Return Kpoints object for kpoints per atom for a cell."""
        if math.fabs((math.floor(kppa ** (1 / 3) + 0.5)) ** 3 - kppa) < 1:
            kppa += kppa * 0.01
        # latt = atoms.lattice_mat
        lengths = atoms.lattice.lat_lengths()
        ngrid = kppa / atoms.num_atoms
        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1 / 3)
        num_div = [int(math.floor(max(mult / lg, 1))) for lg in lengths]
        kpts = Kpoints3D(kpoints=[num_div])
        return kpts

    @property
    def labels(self):
        """Return k-points labels, used for high BZ points."""
        return self._labels

    def write_file(self, filename=""):
        """Write k-point object to a files."""
        if self._kpoint_mode == "automatic":
            f = open(filename, "w")
            f.write("Automatic kpoint scheme\n")
            f.write("0\n")
            line = str(self._header) + "\n"
            f.write(line)
            kp = self._kpoints
            line = (
                str(kp[0][0])
                + str(" ")
                + str(kp[0][1])
                + str(" ")
                + str(kp[0][2])
                + str("\n")
            )
            f.write(line)
            f.close()

        if self._kpoint_mode == "linemode":
            kp = self._kpoints
            lbls = self._labels
            f = open(filename, "w")
            line = str(self._header) + "\n"
            f.write(line)
            line = str(len(kp)) + "\n"
            f.write(line)
            f.write("Reciprocal\n")
            if self._kp_weights == []:
                self._kp_weights = np.ones(len(lbls))
            for i, ii in enumerate(lbls):
                f.write(
                    "%6f %6f %6f %d %s \n"
                    % (kp[i][0], kp[i][1], kp[i][2], self._kp_weights[i], ii)
                )
            f.close()

    def to_dict(self):
        """Provide dictionary representation."""
        d = OrderedDict()
        d["kpoints"] = list(np.array(self._kpoints).tolist())
        d["labels"] = list(self._labels)
        d["kpoint_mode"] = self._kpoint_mode
        d["header"] = self._header
        d["kpoints_weights"] = list(self._kp_weights)
        return d

    @classmethod
    def from_dict(self, d={}):
        """Build class from a dictionary representation."""
        return Kpoints3D(
            kpoints=d["kpoints"],
            labels=d["labels"],
            kpoints_weights=d["kpoints_weights"],
            kpoint_mode=d["kpoint_mode"],
            header=d["header"],
        )

    def high_symm_path(self, atoms):
        """Get high symmetry k-points for given Atoms."""
        spg = Spacegroup3D(atoms=atoms)
        lat_sys = spg.lattice_system
        spg_symb = spg.space_group_symbol
        kp = None
        if lat_sys == "cubic":
            if "P" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().cubic()
            elif "F" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().fcc()
            elif "I" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().bcc()
            else:
                print("kpath space group  is not implemeted ", spg_symb)

        elif lat_sys == "tetragonal":
            if "P" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().tet()
            elif "I" in spg_symb:
                cvn = spg.conventional_standard_structure
                a = cvn.lattice.a
                c = cvn.lattice.c
                if c < a:
                    kp = HighSymmetryKpoint3DFactory().bctet1(c, a)
                else:
                    kp = HighSymmetryKpoint3DFactory().bctet2(c, a)
            else:
                print("kpath space group is not implemeted ", spg_symb)

        elif lat_sys == "orthorhombic":
            cvn = spg.conventional_standard_structure
            a = cvn.lattice.a
            c = cvn.lattice.c
            b = cvn.lattice.b

            if "P" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().orc()

            elif "F" in spg_symb:
                if 1 / a ** 2 > 1 / b ** 2 + 1 / c ** 2:
                    kp = HighSymmetryKpoint3DFactory().orcf1(a, b, c)
                elif 1 / a ** 2 < 1 / b ** 2 + 1 / c ** 2:
                    kp = HighSymmetryKpoint3DFactory().orcf2(a, b, c)
                else:
                    kp = HighSymmetryKpoint3DFactory().orcf3(a, b, c)

            elif "I" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().orci(a, b, c)

            elif "C" in spg_symb or "A" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().orcc(a, b, c)
            else:
                print("kpath space group is not implemeted ", spg_symb)

        elif lat_sys == "hexagonal":
            kp = HighSymmetryKpoint3DFactory().hex()

        elif lat_sys == "rhombohedral":
            prim = spg.primitive_atoms
            alpha = prim.lattice.angles[0]
            if alpha < 90:
                kp = HighSymmetryKpoint3DFactory().rhl1(alpha * pi / 180)
            else:
                kp = HighSymmetryKpoint3DFactory().rhl2(alpha * pi / 180)

        elif lat_sys == "monoclinic":
            cvn = spg.conventional_standard_structure
            a, b, c = cvn.lattice.abc
            alpha = cvn.lattice.angles[0]

            if "P" in spg_symb:
                kp = HighSymmetryKpoint3DFactory().mcl(b, c, alpha * pi / 180)

            elif "C" in spg_symb:
                prim = spg.primitive_atoms.lattice.reciprocal_lattice()

                kgamma = prim.angles[2]
                if kgamma > 90:
                    kp = HighSymmetryKpoint3DFactory().mclc1(
                        a, b, c, alpha * pi / 180
                    )
                if kgamma == 90:
                    kp = HighSymmetryKpoint3DFactory().mclc2(
                        a, b, c, alpha * pi / 180
                    )
                if kgamma < 90:
                    if (
                        b * cos(alpha * pi / 180) / c
                        + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2
                        < 1
                    ):
                        kp = HighSymmetryKpoint3DFactory().mclc3(
                            a, b, c, alpha * pi / 180
                        )
                    if (
                        b * cos(alpha * pi / 180) / c
                        + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2
                        == 1
                    ):
                        kp = HighSymmetryKpoint3DFactory().mclc4(
                            a, b, c, alpha * pi / 180
                        )
                    if (
                        b * cos(alpha * pi / 180) / c
                        + b ** 2 * sin(alpha * pi / 180) ** 2 / a ** 2
                        > 1
                    ):
                        kp = HighSymmetryKpoint3DFactory().mclc5(
                            a, b, c, alpha * pi / 180
                        )
            else:
                print("kpath space group is not implemeted ", spg_symb)

        elif lat_sys == "triclinic":
            prim = spg.primitive_atoms.lattice.reciprocal_lattice()
            kalpha = prim.angles[0]
            kbeta = prim.angles[1]
            kgamma = prim.angles[2]
            if kalpha > 90 and kbeta > 90 and kgamma > 90:
                kp = HighSymmetryKpoint3DFactory().tria()
            if kalpha < 90 and kbeta < 90 and kgamma < 90:
                kp = HighSymmetryKpoint3DFactory().trib()
            if kalpha > 90 and kbeta > 90 and kgamma == 90:
                kp = HighSymmetryKpoint3DFactory().tria()
            if kalpha < 90 and kbeta < 90 and kgamma == 90:
                kp = HighSymmetryKpoint3DFactory().trib()

        else:
            print("kpath space group is not implemeted ", spg_symb)
        # print("kp",spg_symb)
        return kp

    def high_kpath(self, atoms):
        """Get high symmetry path as a dictionary."""
        return self.high_symm_path(atoms).to_dict()

    def interpolated_points(
        self, atoms, line_density=20, coords_are_cartesian=False
    ):
        """Provide bandstructure k-points, controlled by the line_density."""
        list_k_points = []
        sym_point_labels = []
        spg = Spacegroup3D(atoms=atoms)
        prim = spg.primitive_atoms
        self.kpath = self.high_kpath(atoms)
        self._prim_rec = prim.lattice.reciprocal_lattice()
        # print ('self._prim_rec',self._prim_rec)
        for b in self.kpath["path"]:
            for i in range(1, len(b)):
                start = np.array(self.kpath["kpoints"][b[i - 1]])
                end = np.array(self.kpath["kpoints"][b[i]])
                distance = np.linalg.norm(
                    self._prim_rec.cart_coords(start)
                    - self._prim_rec.cart_coords(end)
                )
                nb = int(ceil(distance * line_density))
                if nb == 0:
                    continue
                # print("nb", nb, distance, line_density)
                sym_point_labels.extend([b[i - 1]] + [""] * (nb - 1) + [b[i]])
                list_k_points.extend(
                    [
                        self._prim_rec.cart_coords(start)
                        + float(i)
                        / float(nb)
                        * (
                            self._prim_rec.cart_coords(end)
                            - self._prim_rec.cart_coords(start)
                        )
                        for i in range(0, nb + 1)
                    ]
                )
        if coords_are_cartesian:
            return list_k_points, sym_point_labels
        else:
            frac_k_points = [
                self._prim_rec.frac_coords(k) for k in list_k_points
            ]
            return frac_k_points, sym_point_labels

    def kpath(
        self,
        atoms,
        line_density=20,
        weights=[],
        unique_kp_only=False,
        coords_are_cartesian=False,
    ):
        """Get k-path for bandstructure calculations."""
        k_points, labels = self.interpolated_points(
            atoms,
            line_density=line_density,
            coords_are_cartesian=coords_are_cartesian,
        )
        if unique_kp_only:
            uniqueValues, indicesList = np.unique(
                k_points, axis=0, return_index=True
            )
            k_points = np.array(k_points)[indicesList.astype(int)]
            labels = np.array(labels)[indicesList.astype(int)]
        return Kpoints3D(
            kpoints=k_points,
            header="Non SCF run along symmetry lines",
            labels=labels,
            kpoint_mode="linemode",
        )

    def __repr__(self, indent=4):
        """Representation for print statements."""
        return pprint.pformat(self.to_dict(), indent=indent)


class HighSymmetryKpoint3DFactory(object):
    """High-symmetry k-points for different crystal-systems."""

    def __init__(self, kpoints=[], path=[], name=None):
        """Require kpoints and path."""
        self._kpoints = kpoints
        self._path = path
        self.name = name

    def cubic(self):
        """Cubic HighSymmKPath, return: Dict."""
        self.name = "CUB"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "X": np.array([0.0, 0.5, 0.0]),
            "R": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.5, 0.0]),
        }
        path = [["\\Gamma", "X", "M", "\\Gamma", "R", "X"], ["M", "R"]]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def fcc(self):
        """Fcc HighSymmKPath, return: Dict."""
        self.name = "FCC"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "K": np.array([3.0 / 8.0, 3.0 / 8.0, 3.0 / 4.0]),
            "L": np.array([0.5, 0.5, 0.5]),
            "U": np.array([5.0 / 8.0, 1.0 / 4.0, 5.0 / 8.0]),
            "W": np.array([0.5, 1.0 / 4.0, 3.0 / 4.0]),
            "X": np.array([0.5, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "X", "W", "K", "\\Gamma", "L", "U", "W", "L", "K"],
            ["U", "X"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def bcc(self):
        """Bcc HighSymmKPath, return: Dict."""
        self.name = "BCC"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "H": np.array([0.5, -0.5, 0.5]),
            "P": np.array([0.25, 0.25, 0.25]),
            "N": np.array([0.0, 0.0, 0.5]),
        }
        path = [["\\Gamma", "H", "N", "\\Gamma", "P", "H"], ["P", "N"]]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def to_dict(self):
        """Get dictionary representation."""
        d = OrderedDict()
        d["kpoints"] = self._kpoints
        d["path"] = self._path
        return d

    def tet(self):
        """Tetragonal HighSymmKPath, return: Dict."""
        self.name = "TET"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.5, 0.0]),
            "R": np.array([0.0, 0.5, 0.5]),
            "X": np.array([0.0, 0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "X", "M", "\\Gamma", "Z", "R", "A", "Z"],
            ["X", "R"],
            ["M", "A"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def bctet1(self, c, a):
        """BCT1 HighSymmKPath, return: Dict."""
        self.name = "BCT1"
        eta = (1 + c ** 2 / a ** 2) / 4.0
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "M": np.array([-0.5, 0.5, 0.5]),
            "N": np.array([0.0, 0.5, 0.0]),
            "P": np.array([0.25, 0.25, 0.25]),
            "X": np.array([0.0, 0.0, 0.5]),
            "Z": np.array([eta, eta, -eta]),
            "Z_1": np.array([-eta, 1 - eta, eta]),
        }
        path = [
            ["\\Gamma", "X", "M", "\\Gamma", "Z", "P", "N", "Z_1", "M"],
            ["X", "P"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def bctet2(self, c, a):
        """BCT2 HighSymmKPath, return: Dict."""
        self.name = "BCT2"
        eta = (1 + a ** 2 / c ** 2) / 4.0
        zeta = a ** 2 / (2 * c ** 2)
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "N": np.array([0.0, 0.5, 0.0]),
            "P": np.array([0.25, 0.25, 0.25]),
            "\\Sigma": np.array([-eta, eta, eta]),
            "\\Sigma_1": np.array([eta, 1 - eta, -eta]),
            "X": np.array([0.0, 0.0, 0.5]),
            "Y": np.array([-zeta, zeta, 0.5]),
            "Y_1": np.array([0.5, 0.5, -zeta]),
            "Z": np.array([0.5, 0.5, -0.5]),
        }
        path = [
            [
                "\\Gamma",
                "X",
                "Y",
                "\\Sigma",
                "\\Gamma",
                "Z",
                "\\Sigma_1",
                "N",
                "P",
                "Y_1",
                "Z",
            ],
            ["X", "P"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def orc(self):
        """Orthorhombic HighSymmKPath, return: Dict."""
        self.name = "ORC"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "R": np.array([0.5, 0.5, 0.5]),
            "S": np.array([0.5, 0.5, 0.0]),
            "T": np.array([0.0, 0.5, 0.5]),
            "U": np.array([0.5, 0.0, 0.5]),
            "X": np.array([0.5, 0.0, 0.0]),
            "Y": np.array([0.0, 0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "X", "S", "Y", "\\Gamma", "Z", "U", "R", "T", "Z"],
            ["Y", "T"],
            ["U", "X"],
            ["S", "R"],
        ]

        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def orcf1(self, a, b, c):
        """Orthorhombic f1 HighSymmKPath, return: Dict."""
        self.name = "ORCF1"
        zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4

        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5 + zeta, zeta]),
            "A_1": np.array([0.5, 0.5 - zeta, 1 - zeta]),
            "L": np.array([0.5, 0.5, 0.5]),
            "T": np.array([1, 0.5, 0.5]),
            "X": np.array([0.0, eta, eta]),
            "X_1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([0.5, 0.0, 0.5]),
            "Z": np.array([0.5, 0.5, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
            ["T", "X_1"],
            ["X", "A", "Z"],
            ["L", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def orcf2(self, a, b, c):
        """Orthorhombic f2 HighSymmKPath, return: Dict."""
        self.name = "ORCF2"
        phi = (1 + c ** 2 / b ** 2 - c ** 2 / a ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        delta = (1 + b ** 2 / a ** 2 - b ** 2 / c ** 2) / 4
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "C": np.array([0.5, 0.5 - eta, 1 - eta]),
            "C_1": np.array([0.5, 0.5 + eta, eta]),
            "D": np.array([0.5 - delta, 0.5, 1 - delta]),
            "D_1": np.array([0.5 + delta, 0.5, delta]),
            "L": np.array([0.5, 0.5, 0.5]),
            "H": np.array([1 - phi, 0.5 - phi, 0.5]),
            "H_1": np.array([phi, 0.5 + phi, 0.5]),
            "X": np.array([0.0, 0.5, 0.5]),
            "Y": np.array([0.5, 0.0, 0.5]),
            "Z": np.array([0.5, 0.5, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "C", "D", "X", "\\Gamma", "Z", "D_1", "H", "C"],
            ["C_1", "Z"],
            ["X", "H_1"],
            ["H", "Y"],
            ["L", "\\Gamma"],
        ]

        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def orcf3(self, a, b, c):
        """Orthorhombic f3 HighSymmKPath, return: Dict."""
        self.name = "ORCF3"
        zeta = (1 + a ** 2 / b ** 2 - a ** 2 / c ** 2) / 4
        eta = (1 + a ** 2 / b ** 2 + a ** 2 / c ** 2) / 4
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5 + zeta, zeta]),
            "A_1": np.array([0.5, 0.5 - zeta, 1 - zeta]),
            "L": np.array([0.5, 0.5, 0.5]),
            "T": np.array([1, 0.5, 0.5]),
            "X": np.array([0.0, eta, eta]),
            "X_1": np.array([1, 1 - eta, 1 - eta]),
            "Y": np.array([0.5, 0.0, 0.5]),
            "Z": np.array([0.5, 0.5, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "T", "Z", "\\Gamma", "X", "A_1", "Y"],
            ["X", "A", "Z"],
            ["L", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def orci(self, a, b, c):
        """Orthorhombic I HighSymmKPath, return: Dict."""
        self.name = "ORCI"
        zeta = (1 + a ** 2 / c ** 2) / 4
        eta = (1 + b ** 2 / c ** 2) / 4
        delta = (b ** 2 - a ** 2) / (4 * c ** 2)
        mu = (a ** 2 + b ** 2) / (4 * c ** 2)
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "L": np.array([-mu, mu, 0.5 - delta]),
            "L_1": np.array([mu, -mu, 0.5 + delta]),
            "L_2": np.array([0.5 - delta, 0.5 + delta, -mu]),
            "R": np.array([0.0, 0.5, 0.0]),
            "S": np.array([0.5, 0.0, 0.0]),
            "T": np.array([0.0, 0.0, 0.5]),
            "W": np.array([0.25, 0.25, 0.25]),
            "X": np.array([-zeta, zeta, zeta]),
            "X_1": np.array([zeta, 1 - zeta, -zeta]),
            "Y": np.array([eta, -eta, eta]),
            "Y_1": np.array([1 - eta, eta, -eta]),
            "Z": np.array([0.5, 0.5, -0.5]),
        }
        path = [
            [
                "\\Gamma",
                "X",
                "L",
                "T",
                "W",
                "R",
                "X_1",
                "Z",
                "\\Gamma",
                "Y",
                "S",
                "W",
            ],
            ["L_1", "Y"],
            ["Y_1", "Z"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def orcc(self, a, b, c):
        """Orthorhombic C HighSymmKPath, return: Dict."""
        self.name = "ORCC"
        zeta = (1 + a ** 2 / b ** 2) / 4
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([zeta, zeta, 0.5]),
            "A_1": np.array([-zeta, 1 - zeta, 0.5]),
            "R": np.array([0.0, 0.5, 0.5]),
            "S": np.array([0.0, 0.5, 0.0]),
            "T": np.array([-0.5, 0.5, 0.5]),
            "X": np.array([zeta, zeta, 0.0]),
            "X_1": np.array([-zeta, 1 - zeta, 0.0]),
            "Y": np.array([-0.5, 0.5, 0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            [
                "\\Gamma",
                "X",
                "S",
                "R",
                "A",
                "Z",
                "\\Gamma",
                "Y",
                "X_1",
                "A_1",
                "T",
                "Y",
            ],
            ["Z", "T"],
        ]

        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def hex(self):
        """Hexagonal HighSymmKPath, return: Dict."""
        self.name = "HEX"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.0, 0.0, 0.5]),
            "H": np.array([1.0 / 3.0, 1.0 / 3.0, 0.5]),
            "K": np.array([1.0 / 3.0, 1.0 / 3.0, 0.0]),
            "L": np.array([0.5, 0.0, 0.5]),
            "M": np.array([0.5, 0.0, 0.0]),
        }
        path = [
            ["\\Gamma", "M", "K", "\\Gamma", "A", "L", "H", "A"],
            ["L", "M"],
            ["K", "H"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def rhl1(self, alpha):
        """Rhombohedral 1 HighSymmKPath, return: Dict."""
        self.name = "RHL1"
        eta = (1 + 4 * cos(alpha)) / (2 + 4 * cos(alpha))
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "B": np.array([eta, 0.5, 1.0 - eta]),
            "B_1": np.array([1.0 / 2.0, 1.0 - eta, eta - 1.0]),
            "F": np.array([0.5, 0.5, 0.0]),
            "L": np.array([0.5, 0.0, 0.0]),
            "L_1": np.array([0.0, 0.0, -0.5]),
            "P": np.array([eta, nu, nu]),
            "P_1": np.array([1.0 - nu, 1.0 - nu, 1.0 - eta]),
            "P_2": np.array([nu, nu, eta - 1.0]),
            "Q": np.array([1.0 - nu, nu, 0.0]),
            "X": np.array([nu, 0.0, -nu]),
            "Z": np.array([0.5, 0.5, 0.5]),
        }
        path = [
            ["\\Gamma", "L", "B_1"],
            ["B", "Z", "\\Gamma", "X"],
            ["Q", "F", "P_1", "Z"],
            ["L", "P"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def rhl2(self, alpha):
        """Rhombohedral 2 HighSymmKPath, return: Dict."""
        self.name = "RHL2"
        eta = 1 / (2 * tan(alpha / 2.0) ** 2)
        nu = 3.0 / 4.0 - eta / 2.0
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([0.5, -0.5, 0.0]),
            "L": np.array([0.5, 0.0, 0.0]),
            "P": np.array([1 - nu, -nu, 1 - nu]),
            "P_1": np.array([nu, nu - 1.0, nu - 1.0]),
            "Q": np.array([eta, eta, eta]),
            "Q_1": np.array([1.0 - eta, -eta, -eta]),
            "Z": np.array([0.5, -0.5, 0.5]),
        }
        path = [
            ["\\Gamma", "P", "Z", "Q", "\\Gamma", "F", "P_1", "Q_1", "L", "Z"]
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def mcl(self, b, c, beta):
        """Monoclinic 1 HighSymmKPath, return: Dict."""
        self.name = "MCL"
        eta = (1 - b * cos(beta) / c) / (2 * sin(beta) ** 2)
        nu = 0.5 - eta * c * cos(beta) / b
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "A": np.array([0.5, 0.5, 0.0]),
            "C": np.array([0.0, 0.5, 0.5]),
            "D": np.array([0.5, 0.0, 0.5]),
            "D_1": np.array([0.5, 0.5, -0.5]),
            "E": np.array([0.5, 0.5, 0.5]),
            "H": np.array([0.0, eta, 1.0 - nu]),
            "H_1": np.array([0.0, 1.0 - eta, nu]),
            "H_2": np.array([0.0, eta, -nu]),
            "M": np.array([0.5, eta, 1.0 - nu]),
            "M_1": np.array([0.5, 1 - eta, nu]),
            "M_2": np.array([0.5, 1 - eta, nu]),
            "X": np.array([0.0, 0.5, 0.0]),
            "Y": np.array([0.0, 0.0, 0.5]),
            "Y_1": np.array([0.0, 0.0, -0.5]),
            "Z": np.array([0.5, 0.0, 0.0]),
        }
        path = [
            ["\\Gamma", "Y", "H", "C", "E", "M_1", "A", "X", "H_1"],
            ["M", "D", "Z"],
            ["Y", "D"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def mclc1(self, a, b, c, alpha):
        """Monoclinic C1 HighSymmKPath, return: Dict."""
        self.name = "MCLC1"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F_1": np.array([zeta, zeta, eta]),
            "F_2": np.array([-zeta, -zeta, 1 - eta]),
            # 'F_3': np.array([1 - zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 0.5]),
            "I_1": np.array([1 - phi, phi - 1, 0.5]),
            "L": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "X": np.array([1 - psi, psi - 1, 0.0]),
            "X_1": np.array([psi, 1 - psi, 0.0]),
            "X_2": np.array([psi - 1, -psi, 0.0]),
            "Y": np.array([0.5, 0.5, 0.0]),
            "Y_1": np.array([-0.5, -0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "L", "I"],
            ["I_1", "Z", "F_1"],
            ["Y", "X_1"],
            ["X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def mclc2(self, a, b, c, alpha):
        """Monoclinic C2 HighSymmKPath, return: Dict."""
        self.name = "MCLC2"
        zeta = (2 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        psi = 0.75 - a ** 2 / (4 * b ** 2 * sin(alpha) ** 2)
        phi = psi + (0.75 - psi) * b * cos(alpha) / c
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "F": np.array([1 - zeta, 1 - zeta, 1 - eta]),
            "F_1": np.array([zeta, zeta, eta]),
            "F_2": np.array([-zeta, -zeta, 1 - eta]),
            "F_3": np.array([1 - zeta, -zeta, 1 - eta]),
            "I": np.array([phi, 1 - phi, 0.5]),
            "I_1": np.array([1 - phi, phi - 1, 0.5]),
            "L": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "X": np.array([1 - psi, psi - 1, 0.0]),
            "X_1": np.array([psi, 1 - psi, 0.0]),
            "X_2": np.array([psi - 1, -psi, 0.0]),
            "Y": np.array([0.5, 0.5, 0.0]),
            "Y_1": np.array([-0.5, -0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "L", "I"],
            ["I_1", "Z", "F_1"],
            ["N", "\\Gamma", "M"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def mclc3(self, a, b, c, alpha):
        """Monoclinic C3 HighSymmKPath, return: Dict."""
        self.name = "MCLC3"
        mu = (1 + b ** 2 / a ** 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ** 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F_1": np.array([phi, phi - 1, psi]),
            "F_2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H_1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([0.5, -0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "X": np.array([0.5, -0.5, 0.0]),
            "Y": np.array([mu, mu, delta]),
            "Y_1": np.array([1 - mu, -mu, -delta]),
            "Y_2": np.array([-mu, -mu, -delta]),
            "Y_3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "H", "Z", "I", "F_1"],
            ["H_1", "Y_1", "X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def mclc4(self, a, b, c, alpha):
        """Monoclinic C4 HighSymmKPath, return: Dict."""
        self.name = "MCLC4"
        mu = (1 + b ** 2 / a ** 2) / 4.0
        delta = b * c * cos(alpha) / (2 * a ** 2)
        zeta = mu - 0.25 + (1 - b * cos(alpha) / c) / (4 * sin(alpha) ** 2)
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        phi = 1 + zeta - 2 * mu
        psi = eta - 2 * delta
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([1 - phi, 1 - phi, 1 - psi]),
            "F_1": np.array([phi, phi - 1, psi]),
            "F_2": np.array([1 - phi, -phi, 1 - psi]),
            "H": np.array([zeta, zeta, eta]),
            "H_1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([0.5, -0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "X": np.array([0.5, -0.5, 0.0]),
            "Y": np.array([mu, mu, delta]),
            "Y_1": np.array([1 - mu, -mu, -delta]),
            "Y_2": np.array([-mu, -mu, -delta]),
            "Y_3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "H", "Z", "I"],
            ["H_1", "Y_1", "X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def mclc5(self, a, b, c, alpha):
        """Monoclinic C5 HighSymmKPath, return: Dict."""
        self.name = "MCLC5"
        zeta = (
            b ** 2 / a ** 2 + (1 - b * cos(alpha) / c) / sin(alpha) ** 2
        ) / 4
        eta = 0.5 + 2 * zeta * c * cos(alpha) / b
        mu = (
            eta / 2 + b ** 2 / (4 * a ** 2) - b * c * cos(alpha) / (2 * a ** 2)
        )
        nu = 2 * mu - zeta
        rho = 1 - zeta * a ** 2 / b ** 2
        omega = (
            (4 * nu - 1 - b ** 2 * sin(alpha) ** 2 / a ** 2)
            * c
            / (2 * b * cos(alpha))
        )
        delta = zeta * c * cos(alpha) / b + omega / 2 - 0.25
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "F": np.array([nu, nu, omega]),
            "F_1": np.array([1 - nu, 1 - nu, 1 - omega]),
            "F_2": np.array([nu, nu - 1, omega]),
            "H": np.array([zeta, zeta, eta]),
            "H_1": np.array([1 - zeta, -zeta, 1 - eta]),
            "H_2": np.array([-zeta, -zeta, 1 - eta]),
            "I": np.array([rho, 1 - rho, 0.5]),
            "I_1": np.array([1 - rho, rho - 1, 0.5]),
            "L": np.array([0.5, 0.5, 0.5]),
            "M": np.array([0.5, 0.0, 0.5]),
            "N": np.array([0.5, 0.0, 0.0]),
            "N_1": np.array([0.0, -0.5, 0.0]),
            "X": np.array([0.5, -0.5, 0.0]),
            "Y": np.array([mu, mu, delta]),
            "Y_1": np.array([1 - mu, -mu, -delta]),
            "Y_2": np.array([-mu, -mu, -delta]),
            "Y_3": np.array([mu, mu - 1, delta]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["\\Gamma", "Y", "F", "L", "I"],
            ["I_1", "Z", "H", "F_1"],
            ["H_1", "Y_1", "X", "\\Gamma", "N"],
            ["M", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def tria(self):
        """Trigonal a HighSymmKPath, return: Dict."""
        self.name = "TRI1a"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "L": np.array([0.5, 0.5, 0.0]),
            "M": np.array([0.0, 0.5, 0.5]),
            "N": np.array([0.5, 0.0, 0.5]),
            "R": np.array([0.5, 0.5, 0.5]),
            "X": np.array([0.5, 0.0, 0.0]),
            "Y": np.array([0.0, 0.5, 0.0]),
            "Z": np.array([0.0, 0.0, 0.5]),
        }
        path = [
            ["X", "\\Gamma", "Y"],
            ["L", "\\Gamma", "Z"],
            ["N", "\\Gamma", "M"],
            ["R", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)

    def trib(self):
        """Trigonal b HighSymmKPath, return: Dict."""
        self.name = "TRI1b"
        kpoints = {
            "\\Gamma": np.array([0.0, 0.0, 0.0]),
            "L": np.array([0.5, -0.5, 0.0]),
            "M": np.array([0.0, 0.0, 0.5]),
            "N": np.array([-0.5, -0.5, 0.5]),
            "R": np.array([0.0, -0.5, 0.5]),
            "X": np.array([0.0, -0.5, 0.0]),
            "Y": np.array([0.5, 0.0, 0.0]),
            "Z": np.array([-0.5, 0.0, 0.5]),
        }
        path = [
            ["X", "\\Gamma", "Y"],
            ["L", "\\Gamma", "Z"],
            ["N", "\\Gamma", "M"],
            ["R", "\\Gamma"],
        ]
        return HighSymmetryKpoint3DFactory(kpoints=kpoints, path=path)


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    lattice_mat = Si.lattice_mat
    kp = Kpoints3D().automatic_length_mesh(lattice_mat=lattice_mat)
    # print(kp.__repr__(0))
    kp.write_file("KPOINTS")
    sym = kp.high_symm_path(Si)._path
    # print (sym)
    x, y = kp.interpolated_points(Si)
    # for i,j in zip(x,y):
    #   print (i,j)
    Si = Poscar.from_file(
        "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/CONTCAR"
    ).atoms
    kp = Kpoints3D().kpath(atoms=Si)
    kp.write_file("KPOINTS")
"""
