"""Module to predict X-ray diffraction."""
import numpy as np
import json
import collections
import os
from jarvis.core.specie import Specie

# from jarvis.core.spectrum import Spectrum
import itertools

# from jarvis.analysis.structure.spacegroup import Spacegroup3D,


class XRD(object):
    """Constryct an XRD class."""

    def __init__(
        self,
        wavelength=1.54184,
        thetas=[0, 180],
        two_theta_array=[],
        dhkl_array=[],
        intensity_array=[],
        scaling_factor=100,
        two_theta_tol=1e-5,
        intensity_tol=0.5,
        max_index=5,
    ):
        """Initialize class with incident wavelength."""
        self.wavelength = wavelength
        self.min2theta = np.radians(thetas[0])
        self.max2theta = np.radians(thetas[1])
        self.thetas = thetas
        self.two_theta_array = two_theta_array
        self.dhkl_array = dhkl_array
        self.intensity_array = intensity_array
        self.two_theta_tol = two_theta_tol
        self.intensity_tol = intensity_tol
        self.max_index = max_index
        self.scaling_factor = scaling_factor

    def simulate(self, atoms=None):
        """
        Simulate XRD pattern.

        Forked from https://github.com/qzhu2017/XRD.
        """
        # atoms=Spacegroup3D(atoms).conventional_standard_structure
        rec_matrix = atoms.lattice.reciprocal_lattice_crystallographic().matrix

        min_r = self.wavelength / np.sin(self.max2theta / 2) / 2
        r1 = list(range(1, self.max_index + 1))
        r2 = list(range(-self.max_index, 1))
        r2.reverse()
        r3 = r1 + r2
        r = r3
        hkl_index = [
            miller
            for miller in itertools.product(r, r, r)
            if any([i != 0 for i in miller])
        ]
        # hkl_index=symmetrically_distinct_miller_indices(cvn_atoms=atoms,max_index=self.max_index)
        hkl_max = np.array([self.max_index, self.max_index, self.max_index])
        for index in hkl_index:
            d = np.linalg.norm(np.dot(index, rec_matrix))
            multiple = int(np.ceil(1 / d / min_r))
            index *= multiple
            for i in range(len(hkl_max)):
                if hkl_max[i] < index[i]:
                    hkl_max[i] = index[i]

        h1, k1, l1 = hkl_max
        h = np.arange(-h1, h1 + 1)
        k = np.arange(-k1, k1 + 1)
        L = np.arange(-l1, l1 + 1)

        hkl = np.array((np.meshgrid(h, k, L))).transpose()
        hkl_list = np.reshape(hkl, [len(h) * len(k) * len(L), 3])
        hkl_list = hkl_list[np.where(hkl_list.any(axis=1))[0]]
        d_hkl = 1 / np.linalg.norm(np.dot(hkl_list, rec_matrix), axis=1)

        shortlist = d_hkl > (min_r)
        d_hkl = d_hkl[shortlist]
        hkl_list = hkl_list[shortlist]
        sintheta = self.wavelength / 2 / d_hkl
        self.theta = np.arcsin(sintheta)
        self.hkl_list = np.array(hkl_list)
        self.d_hkl = d_hkl
        fp = os.path.join(
            os.path.dirname(__file__), "atomic_scattering_params.json"
        )
        with open(fp, "r") as f:
            ATOMIC_SCATTERING_PARAMS = json.load(f)
        d0 = (1 / 2 / self.d_hkl) ** 2
        coeffs = []
        zs = []
        for elem in atoms.elements:
            if elem == "D":
                elem = "H"
            c = ATOMIC_SCATTERING_PARAMS[elem]
            z = Specie(elem).Z
            coeffs.append(c)
            zs.append(z)
        coeffs = np.array(coeffs)
        self.peaks = {}
        two_thetas = []

        for hkl, s2, theta, d_hkl in zip(
            self.hkl_list, d0, self.theta, self.d_hkl
        ):

            # Calculate the scattering factor sf
            g_dot_r = np.dot(atoms.frac_coords, np.transpose([hkl])).T[0]
            sf = zs - 41.78214 * s2 * np.sum(
                coeffs[:, :, 0] * np.exp(-coeffs[:, :, 1] * s2), axis=1
            )

            # Calculate the structure factor f
            f = np.sum(sf * np.exp(2j * np.pi * g_dot_r))

            # Calculate the lorentz polarization factor lf
            lf = (1 + np.cos(2 * theta) ** 2) / (
                np.sin(theta) ** 2 * np.cos(theta)
            )

            po = 1
            # Calculate the intensity I
            intensity_tmp = (f * f.conjugate()).real

            # calculate 2*theta
            two_theta = np.degrees(2 * theta)

            # Find where the scattered angles are equal
            ind = np.where(
                np.abs(np.subtract(two_thetas, two_theta)) < self.two_theta_tol
            )

            # Append intensity, hkl plane, and thetas to lists
            if len(ind[0]) > 0:
                self.peaks[two_thetas[ind[0][0]]][0] += intensity_tmp * lf * po
                self.peaks[two_thetas[ind[0][0]]][1].append(tuple(hkl))
            else:
                self.peaks[two_theta] = [
                    intensity_tmp * lf * po,
                    [tuple(hkl)],
                    d_hkl,
                ]
                two_thetas.append(two_theta)

        # max_intensity = max([v[0] for v in self.peaks.values()])
        x = []
        y = []
        d_hkls = []
        hkl_families = []
        for k in sorted(self.peaks.keys()):
            v = self.peaks[k]
            x.append(k)
            y.append(v[0])
            d_hkls.append(v[2])
            family = self.get_unique_families(v[1])
            hkl_families.append(family)
        x = np.array(x)
        y = np.array(y)
        d_hkls = np.array(d_hkls)
        const = np.sum(y, axis=0)
        y = y * self.scaling_factor / const
        screen = y > self.intensity_tol
        self.intensity_array = y[screen]
        self.two_theta_array = x[screen]
        self.dhkl_array = d_hkls[screen]
        return (
            x[screen].tolist(),
            d_hkls[screen].tolist(),
            y[screen].tolist(),
            # hkl_families[screen],
        )

    def get_unique_families(self, hkls):
        """
        Return unique families of Miller indices.

        Families must be permutations of each other.
        Args:
            hkls ([h, k, l]): List of Miller indices.
        Returns:
            {hkl: multiplicity}: A dict with unique hkl and multiplicity.
        """

        def is_perm(hkl1, hkl2):
            h1 = np.abs(hkl1)
            h2 = np.abs(hkl2)
            return all([i == j for i, j in zip(sorted(h1), sorted(h2))])

        unique = collections.defaultdict(list)
        for hkl1 in hkls:
            found = False
            for hkl2 in unique.keys():
                if is_perm(hkl1, hkl2):
                    found = True
                    unique[hkl2].append(hkl1)
                    break
            if not found:
                unique[hkl1].append(hkl1)

        pretty_unique = {}
        for k, v in unique.items():
            pretty_unique[sorted(v)[-1]] = len(v)

        return pretty_unique


"""
if __name__ == "__main__":
    # h=create_index()
    # print (h,len(h))
    from jarvis.core.atoms import Atoms

    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
    elements = ["Si", "Si"]
    atoms = Atoms(lattice_mat=box, coords=coords, elements=elements)
    a, b, c, d = XRD().simulate(atoms=atoms)
    # print("theta,d_hkls,intens", a, b, c)
    print("a=", a)
    print("b=", b)
    print("c=", c)
"""
