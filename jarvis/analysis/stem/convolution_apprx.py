"""Module to simulate STEM images using convoltuin approximation."""
# Adapted from https://github.com/jacobjma/fourier-scale-calibration
import numpy as np
from scipy.interpolate import interp1d

# from numbers import Number
from jarvis.core.utils import gaussian
from jarvis.core.utils import lorentzian2 as lorentzian
from jarvis.core.atoms import Atoms  # , get_supercell_dims, crop_square
from typing import List


class STEMConv(object):
    """Module to simulate STEM images using convoltuin approximation."""

    def __init__(
        self,
        atoms=None,
        output_size=[50, 50],
        power_factor=1.7,
        gaussian_width=0.5,
        lorentzian_width=0.5,
        intensity_ratio=0.5,
        nbins=100,
        tol=0.5,
        crop=False,
    ):
        """Intitialize the class."""
        self.atoms = atoms
        self.output_size = output_size
        self.power_factor = power_factor
        self.gaussian_width = gaussian_width
        self.lorentzian_width = lorentzian_width
        self.intensity_ratio = intensity_ratio
        self.nbins = nbins
        self.tol = tol
        self.crop = crop

    def superpose_deltas(self, positions, array):
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
        return array

    def simulate_surface(
        self,
        atoms: Atoms,
        px_scale: float = 0.2,
        eps: float = 0.6,
        rot: float = 0,
        shift: List = [0, 0],
    ):
        """Simulate a STEM image.

        atoms: jarvis.core.Atoms material slab
        px_scale: pixel size in angstroms/px
        eps: tolerance factor (angstroms)
        for rendering atoms outside the field of view
        rot: rotation about the image center (degrees)
        shift: rigid translation of field of view [dx, dy] (angstroms)

        """
        shift = np.squeeze(shift)
        output_px = np.squeeze(self.output_size)  # px

        # field of view size in angstroms
        view_size = px_scale * (output_px - 1)

        # construct a supercell grid big enough to fill the field of view
        cell_extent = atoms.lattice.abc[0:2]  # np.diag(atoms.lattice_mat)[:2]

        cells = ((view_size // cell_extent) + 1).astype(int)
        # print ('cells',cells)
        atoms = atoms.make_supercell_matrix((3 * cells[0], 3 * cells[1], 1))

        # Set up real-space grid (in angstroms)
        # construct the probe array with the output target size
        # fftshift, pad, un-fftshift
        x = np.fft.fftfreq(output_px[0]) * output_px[0] * px_scale
        y = np.fft.fftfreq(output_px[1]) * output_px[1] * px_scale
        r = np.sqrt(x[:, None] ** 2 + y[None] ** 2)

        # construct the probe profile centered
        # at (0,0) on the periodic spatial grid
        x = np.linspace(0, 4 * self.lorentzian_width, self.nbins)
        profile = gaussian(
            x, self.gaussian_width
        ) + self.intensity_ratio * lorentzian(x, self.lorentzian_width)
        profile /= profile.max()
        f = interp1d(x, profile, fill_value=0, bounds_error=False)
        intensity = f(r)

        # shift the probe profile to the center
        # apply zero-padding, and shift back to the origin
        margin = int(np.ceil(5 / px_scale))  # int like 20
        intensity = np.fft.fftshift(intensity)
        intensity = np.pad(intensity, (margin, margin))
        intensity = np.fft.fftshift(intensity)

        # project atomic coordinates onto the image
        # center them as well
        centroid = np.mean(atoms.cart_coords[:, :2], axis=0)

        # center atom positions around (0,0)
        pos = atoms.cart_coords[:, :2] - centroid

        # apply field of view rotation
        # (actually rotate the lattice coordinates)
        if rot != 0:
            rot = np.radians(rot)
            R = np.array(
                [[np.cos(rot), -np.sin(rot)], [np.sin(rot), np.cos(rot)]]
            )
            pos = pos @ R

        # shift to center of image
        pos += view_size / 2

        # apply rigid translation of atoms wrt image field of view
        pos += shift

        # select only atoms in field of view
        in_view = (
            (pos[:, 0] > -eps)
            & (pos[:, 0] < view_size[0] + eps)
            & (pos[:, 1] > -eps)
            & (pos[:, 1] < view_size[1] + eps)
        )
        # pos = pos[in_view]
        numbers = np.array(atoms.atomic_numbers)
        # numbers = numbers[in_view]

        atom_px = pos / px_scale  # AA / (AA/px) -> px

        # atom_px = atom_px + margin

        render = in_view
        # render = (
        #     (pos[:, 0] > 0)
        #     & (pos[:, 0] < view_size[0])
        #     & (pos[:, 1] > 0)
        #     & (pos[:, 1] < view_size[1])
        # )

        numbers_render = numbers[render]
        # # shift atomic positions to offset zero padding
        atom_px_render = atom_px[render] + margin

        # initialize arrays with zero padding
        array = np.zeros((1,) + intensity.shape)  # adding extra 1
        mask = np.zeros((1,) + intensity.shape)
        # print(f"intensity: {array.shape}")
        for number in np.unique(np.array(atoms.atomic_numbers)):

            temp = np.zeros((1,) + intensity.shape)
            temp = self.superpose_deltas(
                atom_px_render[numbers_render == number], temp
            )
            array += temp * number ** self.power_factor
            temp = np.where(temp > 0, number, temp)
            mask += temp[0]

        # FFT convolution of beam profile and atom position delta functions
        array = np.fft.ifft2(np.fft.fft2(array) * np.fft.fft2(intensity)).real

        # crop the FFT padding and fix atom coordinates relative to
        # the image field of view
        sel = slice(margin, -margin)
        array = array[0, sel, sel]
        mask = mask[0, sel, sel]
        # atom_px = atom_px - margin

        atom_px = pos[in_view] / px_scale
        numbers = numbers[in_view]

        return array, mask, atom_px, numbers


"""
if __name__ == "__main__":
    from jarvis.core.atoms import crop_squre
    from jarvis.db.figshare import data, get_jid_data
    import matplotlib.pyplot as plt
    from jarvis.core.atoms import Atoms, ase_to_atoms, get_supercell_dims

    plt.switch_backend("agg")

    a = Atoms.from_dict(get_jid_data("JVASP-667")["atoms"])
    stem = STEMConv(output_size=(256, 256))
    p, mask, atom_x, nb = stem.simulate_surface(
        a, px_scale=0.15, eps=0.6, rot=2, shift=[-0.2, 0.3]
    )
    plt.imshow(p, origin="lower", cmap="plasma")
    plt.scatter(atom_x[:, 1], atom_x[:, 0], marker="o",
    facecolors="none", edgecolors="r"
    )
    plt.savefig("stem_example.png")
    plt.close()
"""
