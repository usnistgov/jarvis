"""Module to simulate STEM images using convoltuin approximation."""
# Adapted from https://github.com/jacobjma/fourier-scale-calibration
import numpy as np
from scipy.interpolate import interp1d
from numbers import Number
from jarvis.core.utils import gaussian
from jarvis.core.utils import lorentzian2 as lorentzian


class STEMConv(object):
    """Module to simulate STEM images using convoltuin approximation."""

    def __init__(
        self,
        atoms=None,
        output_size=[50, 50],
        power_factor=1.4,
        gaussian_width=0.5,
        lorentzian_width=0.5,
        intensity_ratio=0.5,
        nbins=100,
        tol=0.5,
    ):
        """
        Intitialize the class.
        """
        self.atoms = atoms
        self.output_size = output_size
        self.power_factor = power_factor
        self.gaussian_width = gaussian_width
        self.lorentzian_width = lorentzian_width
        self.intensity_ratio = intensity_ratio
        self.nbins = nbins
        self.tol = tol

    @staticmethod
    def superpose_deltas(positions, array):
        """Superpose deltas."""
        z = 0
        shape = array.shape[-2:]
        rounded = np.floor(positions).astype(np.int32)
        rows, cols = rounded[:, 0], rounded[:, 1]

        array[z, rows, cols] += (1 - (positions[:, 0] - rows)) * (
            1 - (positions[:, 1] - cols)
        )
        array[z, (rows + 1) % shape[0], cols] += (positions[:, 0] - rows) * (
            1 - (positions[:, 1] - cols)
        )
        array[z, rows, (cols + 1) % shape[1]] += (
            1 - (positions[:, 0] - rows)
        ) * (positions[:, 1] - cols)
        array[z, (rows + 1) % shape[0], (cols + 1) % shape[1]] += (
            rows - positions[:, 0]
        ) * (cols - positions[:, 1])

    def simulate_surface(self):
        """Simulate a STEM image."""

        extent = np.diag(self.atoms.lattice_mat)[:2]

        shape = 1
        sampling = (
            extent[0] / self.output_size[0],
            extent[1] / self.output_size[1],
        )

        margin = int(np.ceil(5 / min(sampling)))  # int like 20
        shape_w_margin = (
            self.output_size[0] + 2 * margin,
            self.output_size[1] + 2 * margin,
        )
        # Set up a grid
        x = np.fft.fftfreq(shape_w_margin[0]) * shape_w_margin[1] * sampling[0]
        y = np.fft.fftfreq(shape_w_margin[1]) * shape_w_margin[1] * sampling[1]
        r = np.sqrt(x[:, None] ** 2 + y[None] ** 2)

        # proble profile

        x = np.linspace(0, 4 * self.lorentzian_width, self.nbins)
        profile = gaussian(
            x, self.gaussian_width
        ) + self.intensity_ratio * lorentzian(x, self.lorentzian_width)

        profile /= profile.max()
        f = interp1d(x, profile, fill_value=0, bounds_error=False)
        intensity = f(r)
        positions = self.atoms.cart_coords[:, :2] / sampling - self.tol
        # Check if atoms are within the specified range
        inside = (
            (positions[:, 0] > -margin)
            & (positions[:, 1] > -margin)
            & (positions[:, 0] < self.output_size[0] + margin)
            & (positions[:, 1] < self.output_size[1] + margin)
        )
        positions = positions[inside] + margin
        numbers = np.array(self.atoms.atomic_numbers)[inside]

        array = np.zeros((1,) + shape_w_margin)  # adding extra 1
        for number in np.unique(np.array(self.atoms.atomic_numbers)):
            temp = np.zeros((1,) + shape_w_margin)
            self.superpose_deltas(positions[numbers == number], temp)

            array += temp * number ** self.power_factor

        array = np.fft.ifft2(np.fft.fft2(array) * np.fft.fft2(intensity)).real
        array = array[0, margin:-margin, margin:-margin]
        return array


"""
if __name__ == "__main__":
    from jarvis.core.atoms import crop_squre
    from jarvis.db.figshare import data, get_jid_data
    import matplotlib.pyplot as plt
    from jarvis.core.atoms import Atoms, ase_to_atoms, get_supercell_dims

    plt.switch_backend("agg")

    a = Atoms.from_dict(get_jid_data("JVASP-667")["atoms"])
    c = crop_square(a)
    # c = a.make_supercell_matrix([2, 2, 1])
    p = STEMConv(atoms=c).simulate_surface()
    plt.imshow(p, interpolation="gaussian", cmap="plasma")
    plt.savefig("stem_example.png")
    plt.close()
"""
