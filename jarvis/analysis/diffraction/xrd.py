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

def baseline_als(y, lam, p, niter=10):
    """ALS baseline correction to remove broad background trends."""
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z

def recast_array(x_original, y_original, x_new, tol=0.1):
    """Recast original spectrum onto a new grid, accumulating close values."""
    x_original = np.array(x_original)
    y_new = np.zeros_like(x_new, dtype=np.float64)

    # Accumulate intensities for sharpness preservation
    for x_val, y_val in zip(x_original, y_original):
        closest_index = np.abs(x_new - x_val).argmin()
        y_new[closest_index] += y_val

    # Remove noise below tolerance level
    y_new[y_new < tol] = 0
    return x_new, y_new

def sharpen_peaks(y, sigma=0.5):
    """Sharpen peaks using a narrow Gaussian filter."""
    # Use a very small sigma to reduce peak broadening
    y_sharp = gaussian_filter1d(y, sigma=sigma, mode='constant')
    return y_sharp

def processed(x, y, x_range=[0, 90], intvl=0.1, sigma=.05,recast=True,tol=0.1,background_subs=True):
    """Process the spectrum: background removal and peak sharpening."""
    y = np.array(y,dtype='float')
    if background_subs:

      # 1. Baseline correction
      background = baseline_als(y, lam=10000, p=0.01)
      y_corrected = y - background
    else:
      y_corrected = y

    # 2. Normalize the corrected spectrum
    y_corrected = y_corrected / np.max(y_corrected)

    # 3. Generate new x-axis values
    x_new = np.arange(x_range[0], x_range[1], intvl)

    # 4. Recast the spectrum onto the new grid
    if recast:
       x_new, y_corrected = recast_array(x, y_corrected, x_new,tol=tol)



    # 5. Sharpen the peaks using Gaussian filtering
    y_sharp = sharpen_peaks(y_corrected, sigma=sigma)

    # 6. Final normalization
    if np.max(y_sharp) > 0:
        y_sharp = y_sharp / np.max(y_sharp)

    return x_new, y_sharp

def smooth_xrd(atoms=None,thetas=[0, 90],intvl=0.5):
    a, b, c = XRD(thetas=thetas).simulate(atoms=atoms)
    a = np.array(a)
    c = np.array(c)
    c=c/np.max(c)
    a, c = recast_array(
        x_original=a,
        y_original=c,
        x_new=np.arange(thetas[0], thetas[1], intvl),
    )
    c=c/np.max(c)
    #c_str = "\n".join(["{0:.3f}".format(x) for x in c])
    c_str = "\n".join(["{0:.2f}".format(x) for x in c])

    return c_str,c

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
