"""Modules for handing crystallographic lattice-parameters."""

import numpy as np
from numpy import dot
from collections import OrderedDict


def abs_cap(val, max_abs_val=1):
    """
    Return the value with its absolute value capped at max_abs_val.

    Particularly useful in passing values to trignometric functions where
    numerical errors may result in an argument > 1 being passed in.

    Args:

        val (float): Input value.

        max_abs_val (float): The maximum absolute value for val. Defaults to 1.

    Returns:
        val if abs(val) < 1 else sign of val * max_abs_val.
    """
    return max(min(val, max_abs_val), -max_abs_val)


class Lattice(object):
    """Construct Lattice parameter object."""

    def __init__(self, lattice_mat=None, round_off=5):
        """
        Hold lattice parameter information.

        >>> box=[[10,0,0],[0,10,0],[0,0,10]]
        >>> lat=Lattice(box)
        >>> lat.lat_lengths()
        [10.0, 10.0, 10.0]
        >>> Lattice(box).lattice()[0][0]
        10.0
        >>> lat.lat_angles()
        [90.0, 90.0, 90.0]
        >>> round(lat.inv_lattice()[0][0],2)
        0.1
        >>> [round(i,2) for i in lat.lat_angles(radians=True)]
        [1.57, 1.57, 1.57]
        >>> frac_coords=[[0,0,0],[0.5,0.5,0.5]]
        >>> lat.cart_coords(frac_coords)[1][1]
        5.0
        >>> cart_coords=[[0,0,0],[5,5,5]]
        >>> lat.frac_coords(cart_coords)[1][1]
        0.5
        >>> box=[[1e-8,0,0],[0,1e-8,0],[0,0,1e-8]]
        >>> lat=Lattice(box)
        >>> [round(i,2) for i in lat.lat_angles()]
        [90.0, 90.0, 90.0]
        """
        tmp = np.array(lattice_mat, dtype=np.float64).reshape((3, 3))
        self._lat = np.around(tmp, decimals=round_off)
        self._inv_lat = None
        self._lll_matrix_mappings = {}

    def lat_lengths(self):
        """Return lattice vectors' lengths."""
        return [
            round(i, 6)
            for i in (np.sqrt(np.sum(self._lat ** 2, axis=1)).tolist())
        ]
        # return [round(np.linalg.norm(v), 6) for v in self._lat]

    @property
    def volume(self):
        """Return volume given a lattice object."""
        m = self._lat
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol

    @property
    def a(self):
        """Return a lattice vector length."""
        return self.lat_lengths()[0]

    @property
    def b(self):
        """Return b lattice vector length."""
        return self.lat_lengths()[1]

    @property
    def c(self):
        """Return c lattice vector length."""
        return self.lat_lengths()[2]

    @property
    def alpha(self):
        """Return alpha lattice vector angle."""
        return self.lat_angles()[0]

    @property
    def beta(self):
        """Return beta lattice vector angle."""
        return self.lat_angles()[1]

    @property
    def gamma(self):
        """Return gamma lattice vector angle."""
        return self.lat_angles()[2]

    @property
    def abc(self):
        """Return lattice vector lengths."""
        return self.lat_lengths()

    @property
    def angles(self):
        """Return lattice vector angles."""
        return self.lat_angles()

    def lat_angles(self, tol=1e-2, radians=False):
        """Return lattice vector angles in radians or degree."""
        lengths = self.lat_lengths()
        angles = []
        for i in range(3):
            j = i - 1
            k = i - 2
            leng2 = lengths[j] * lengths[k]
            if leng2 > tol:
                tmp = np.dot(self._lat[j], self._lat[k]) / leng2
                angle = round(180.0 * np.arccos(tmp) / np.pi, 4)
            else:
                angle = 90.0
            angles.append(angle)
        if radians:
            angles = [round(angleX * np.pi / 180.0, 4) for angleX in angles]
        return angles

    @property
    def parameters(self):
        """Return lattice vector angles in radians or degree."""
        return [
            self.a,
            self.b,
            self.c,
            self.angles[0],
            self.angles[1],
            self.angles[2],
        ]

    @staticmethod
    def from_parameters(a, b, c, alpha, beta, gamma):
        """Construct Lattice from lattice parameter information."""
        angles_r = np.radians([alpha, beta, gamma])
        cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
        sin_alpha, sin_beta, sin_gamma = np.sin(angles_r)
        tmp = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta)
        val = abs_cap(tmp)
        gamma_star = np.arccos(val)
        vector_a = [a * sin_beta, 0.0, a * cos_beta]
        vector_b = [
            -b * sin_alpha * np.cos(gamma_star),
            b * sin_alpha * np.sin(gamma_star),
            b * cos_alpha,
        ]
        vector_c = [0.0, 0.0, float(c)]
        return Lattice(lattice_mat=[vector_a, vector_b, vector_c])

    @staticmethod
    def cubic(a):
        """Construct cubic Lattice from lattice parameter information."""
        return Lattice.from_parameters(a, a, a, 90, 90, 90)

    @staticmethod
    def tetragonal(a, c):
        """Construct tetragonal Lattice from lattice parameter information."""
        return Lattice.from_parameters(a, a, c, 90, 90, 90)

    @staticmethod
    def orthorhombic(a, b, c):
        """Construct orthorhombic Lattice."""
        return Lattice.from_parameters(a, b, c, 90, 90, 90)

    @staticmethod
    def monoclinic(a, b, c, beta):
        """Construct monoclinic Lattice from lattice parameter information."""
        return Lattice.from_parameters(a, b, c, 90, beta, 90)

    @staticmethod
    def hexagonal(a, c):
        """Construct hexagonal Lattice from lattice parameter information."""
        return Lattice.from_parameters(a, a, c, 90, 90, 120)

    @staticmethod
    def rhombohedral(a, alpha):
        """Construct rhombohedral Lattice."""
        return Lattice.from_parameters(a, a, a, alpha, alpha, alpha)

    def to_dict(self):
        """Return lattice parameter information as a dictionary."""
        d = OrderedDict()
        d["matrix"] = self.matrix
        return d

    @classmethod
    def from_dict(self, d):
        """Construct Lattice from lattice matrix dictionary."""
        return Lattice(lattice_mat=d["matrix"])

    @property
    def matrix(self):
        """Return lattice matrix."""
        return self.lattice()

    @property
    def inv_matrix(self):
        """Return inverse lattice matrix."""
        return self.inv_lattice()

    def lattice(self):
        """Return lattice matrix."""
        return self._lat

    def inv_lattice(self):
        """Return inverse lattice matrix."""
        # if self._inv_lat == None:
        self._inv_lat = np.linalg.inv(self._lat)
        return self._inv_lat

    def cart_coords(self, frac_coords):
        """Return cartesian coords from fractional coords using Lattice."""
        return np.dot(np.array(frac_coords), self._lat)

    def frac_coords(self, cart_coords):
        """Return fractional coords from cartesian coords using Lattice."""
        return np.dot(np.array(cart_coords), self.inv_lattice())

    def reciprocal_lattice(self):
        """Return reciprocal Lattice."""
        return Lattice(2 * np.pi * np.linalg.inv(self._lat).T)

    def reciprocal_lattice_crystallographic(self):
        """Return reciprocal Lattice without 2 * pi."""
        return Lattice(self.reciprocal_lattice().matrix / (2 * np.pi))

    def get_points_in_sphere(self, frac_points, center, r):
        """
        Find all points within a sphere from the point.

        Takes into account periodic boundary conditions.
        This includes sites in other periodic images.
        Adapted from pymatgen.
        """
        recp_len = np.array(self.reciprocal_lattice().lat_lengths()) / (
            2 * np.pi
        )
        nmax = float(r) * recp_len + 0.01

        # Get the fractional coordinates of the center of the sphere
        pcoords = self.frac_coords(center)
        center = np.array(center)

        # Prepare the list of output atoms
        n = len(frac_points)
        fcoords = np.array(frac_points) % 1
        indices = np.arange(n)

        # Generate all possible images that could be within `r` of `center`
        mins = np.floor(pcoords - nmax)
        maxes = np.ceil(pcoords + nmax)
        arange = np.arange(start=mins[0], stop=maxes[0], dtype=np.int)
        brange = np.arange(start=mins[1], stop=maxes[1], dtype=np.int)
        crange = np.arange(start=mins[2], stop=maxes[2], dtype=np.int)
        arange = arange[:, None] * np.array([1, 0, 0], dtype=np.int)[None, :]
        brange = brange[:, None] * np.array([0, 1, 0], dtype=np.int)[None, :]
        crange = crange[:, None] * np.array([0, 0, 1], dtype=np.int)[None, :]
        tmp_cr = crange[None, None, :]
        images = arange[:, None, None] + brange[None, :, None] + tmp_cr

        # Generate the coordinates of all atoms within these images
        tmp_img = images[None, :, :, :, :]
        shifted_coords = fcoords[:, None, None, None, :] + tmp_img

        # Determine distance from `center`
        cart_coords = self.cart_coords(fcoords)
        cart_images = self.cart_coords(images)
        tmp_img = cart_images[None, :, :, :, :]
        coords = cart_coords[:, None, None, None, :] + tmp_img
        coords -= center[None, None, None, None, :]
        coords **= 2
        d_2 = np.sum(coords, axis=4)

        # Determine which points are within `r` of `center`
        within_r = np.where(d_2 <= r ** 2)
        return (
            shifted_coords[within_r],
            np.sqrt(d_2[within_r]),
            indices[within_r[0]],
            images[within_r[1:]],
        )

    def find_all_matches(self, other_lattice, ltol=1e-5, atol=1):
        """Find all lattice mappings, adapted from pymatgen."""
        lengths = other_lattice.lat_lengths()
        angles = other_lattice.lat_angles()
        # print ('angles',angles)
        alpha = angles[0]
        beta = angles[1]
        gamma = angles[2]
        frac, dist, _, _ = self.get_points_in_sphere(
            [[0, 0, 0]], [0, 0, 0], max(lengths) * (1 + ltol)
        )
        cart = self.cart_coords(frac)
        # this can't be broadcast because they're different lengths
        inds = [
            np.logical_and(
                dist / leng < 1 + ltol, dist / leng > 1 / (1 + ltol)
            )  # type: ignore
            for leng in lengths
        ]
        c_a, c_b, c_c = (cart[i] for i in inds)
        f_a, f_b, f_c = (frac[i] for i in inds)
        l_a, l_b, l_c = (
            np.sum(c ** 2, axis=-1) ** 0.5 for c in (c_a, c_b, c_c)
        )

        def get_angles(v1, v2, l1, l2):
            x = np.inner(v1, v2) / l1[:, None] / l2
            x[x > 1] = 1
            x[x < -1] = -1
            angles = np.arccos(x) * 180.0 / np.pi
            return angles

        alphab = np.abs(get_angles(c_b, c_c, l_b, l_c) - alpha) < atol
        betab = np.abs(get_angles(c_a, c_c, l_a, l_c) - beta) < atol
        gammab = np.abs(get_angles(c_a, c_b, l_a, l_b) - gamma) < atol

        for i, all_j in enumerate(gammab):
            inds = np.logical_and(
                all_j[:, None], np.logical_and(alphab, betab[i][None, :])
            )
            for j, k in np.argwhere(inds):
                scale_m = np.array(
                    (f_a[i], f_b[j], f_c[k]), dtype=np.int
                )  # type: ignore
                if abs(np.linalg.det(scale_m)) < 1e-8:
                    continue

                aligned_m = np.array((c_a[i], c_b[j], c_c[k]))

                rotation_m = np.linalg.solve(aligned_m, other_lattice._lat)
                yield Lattice(aligned_m), rotation_m, scale_m

    def find_matches(self, other_lattice, ltol=1e-5, atol=1):
        """Find matches with length and angle tolerances."""
        for x in self.find_all_matches(other_lattice, ltol, atol):
            return x
        return None

    def _calculate_lll(self, delta=0.75):
        """
        Perform a Lenstra-Lenstra-Lovasz lattice basis reduction.

        Obtain a c-reduced basis.
        This method returns a basis which is as "good" as
        possible, with "good" defined by orthongonality of the lattice vectors.
        This basis is used for all the periodic boundary condition calcs.
        Adapted from pymatgen.

        Args:

            delta (float): Reduction parameter.
                           Default of 0.75 is usually fine.

        Returns:
            Reduced lattice matrix, mapping to get to that lattice.
        """
        # Transpose the lattice matrix first so that basis vectors are columns.
        # Makes life easier.
        a = self._lat.copy().T

        b = np.zeros((3, 3))  # Vectors after the Gram-Schmidt process
        u = np.zeros((3, 3))  # Gram-Schmidt coeffieicnts
        m = np.zeros(3)  # These are the norm squared of each vec.

        b[:, 0] = a[:, 0]
        m[0] = dot(b[:, 0], b[:, 0])
        for i in range(1, 3):
            u[i, 0:i] = dot(a[:, i].T, b[:, 0:i]) / m[0:i]
            b[:, i] = a[:, i] - dot(b[:, 0:i], u[i, 0:i].T)
            m[i] = dot(b[:, i], b[:, i])

        k = 2

        mapping = np.identity(3, dtype=np.double)
        while k <= 3:
            # Size reduction.
            for i in range(k - 1, 0, -1):
                q = round(u[k - 1, i - 1])
                if q != 0:
                    # Reduce the k-th basis vector.
                    a[:, k - 1] = a[:, k - 1] - q * a[:, i - 1]
                    mapping[:, k - 1] = (
                        mapping[:, k - 1] - q * mapping[:, i - 1]
                    )
                    uu = list(u[i - 1, 0 : (i - 1)])
                    uu.append(1)
                    # Update the GS coefficients.
                    u[k - 1, 0:i] = u[k - 1, 0:i] - q * np.array(uu)

            # Check the Lovasz condition.
            if dot(b[:, k - 1], b[:, k - 1]) >= (
                delta - abs(u[k - 1, k - 2]) ** 2
            ) * dot(b[:, (k - 2)], b[:, (k - 2)]):
                # Increment k if the Lovasz condition holds.
                k += 1
            else:
                # If the Lovasz condition fails,
                # swap the k-th and (k-1)-th basis vector
                v = a[:, k - 1].copy()
                a[:, k - 1] = a[:, k - 2].copy()
                a[:, k - 2] = v

                v_m = mapping[:, k - 1].copy()
                mapping[:, k - 1] = mapping[:, k - 2].copy()
                mapping[:, k - 2] = v_m

                # Update the Gram-Schmidt coefficients
                for s in range(k - 1, k + 1):
                    u[s - 1, 0 : (s - 1)] = (
                        dot(a[:, s - 1].T, b[:, 0 : (s - 1)]) / m[0 : (s - 1)]
                    )
                    b[:, s - 1] = a[:, s - 1] - dot(
                        b[:, 0 : (s - 1)], u[s - 1, 0 : (s - 1)].T
                    )
                    m[s - 1] = dot(b[:, s - 1], b[:, s - 1])

                if k > 2:
                    k -= 1
                else:
                    # We have to do p/q, so do lstsq(q.T, p.T).T instead.
                    p = dot(a[:, k:3].T, b[:, (k - 2) : k])
                    q = np.diag(m[(k - 2) : k])
                    result = np.linalg.lstsq(q.T, p.T, rcond=None)[
                        0
                    ].T  # type: ignore
                    u[k:3, (k - 2) : k] = result

        return a.T, mapping.T

    def get_lll_reduced_lattice(self, delta=0.75):
        """
        Get LLL reduced lattice.

        Adpted from pymatgen.
        Args:

           delta: Delta parameter.

        Returns:
             LLL reduced Lattice.
        """
        if delta not in self._lll_matrix_mappings:
            self._lll_matrix_mappings[delta] = self._calculate_lll()
        return Lattice(self._lll_matrix_mappings[delta][0])


def lattice_coords_transformer(
    old_lattice_mat=[], new_lattice_mat=[], cart_coords=[]
):
    """Transform coords to a new lattice."""
    M = np.linalg.solve(new_lattice_mat, old_lattice_mat)
    #  Maintains the z-distances - kfg
    M[2, 2] = 1.0
    new_cart_coords = np.dot(cart_coords, M)
    return new_cart_coords


"""
if __name__ == "__main__":

    box = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    lat = Lattice(box)
    print("lll", lat._calculate_lll())
    print("lll_educed", lat.get_lll_reduced_lattice()._lat)
    frac_coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
    print(lat.cart_coords(frac_coords)[1][1])

    cart_coords = [[0, 0, 0], [5, 5, 5]]
    print(lat.frac_coords(cart_coords)[1][1])
"""
