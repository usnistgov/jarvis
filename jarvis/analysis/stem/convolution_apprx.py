"""Module to simulate STEM images using convoltuin approximation."""
# Adapted from https://github.com/jacobjma/fourier-scale-calibration
import numpy as np
from scipy.interpolate import interp1d
from numbers import Number

gaussian = lambda x, sigma: np.exp(-(x ** 2) / (2 * sigma ** 2))
lorentzian = (
    lambda x, gamma: gamma
    / 2
    / (np.pi * (x ** 2 + (gamma / 2) ** 2))
    / (2 / (np.pi * gamma))
)


def probe_profile(
    gaussian_width=1, lorentzian_width=1, intensity_ratio=1, n=100
):
    """Probe STEM profile."""
    gaussian_width = gaussian_width / 2.355

    x = np.linspace(0, 4 * lorentzian_width, n)
    profile = gaussian(x, gaussian_width) + intensity_ratio * lorentzian(
        x, lorentzian_width
    )

    profile /= profile.max()
    f = interp1d(x, profile, fill_value=0, bounds_error=False)
    return f


def superpose_deltas(positions, z, array):
    """Probe STEM profile."""
    shape = array.shape[-2:]
    rounded = np.floor(positions).astype(np.int32)
    rows, cols = rounded[:, 0], rounded[:, 1]

    array[z, rows, cols] += (1 - (positions[:, 0] - rows)) * (
        1 - (positions[:, 1] - cols)
    )
    array[z, (rows + 1) % shape[0], cols] += (positions[:, 0] - rows) * (
        1 - (positions[:, 1] - cols)
    )
    array[z, rows, (cols + 1) % shape[1]] += (1 - (positions[:, 0] - rows)) * (
        positions[:, 1] - cols
    )
    array[z, (rows + 1) % shape[0], (cols + 1) % shape[1]] += (
        rows - positions[:, 0]
    ) * (cols - positions[:, 1])


def simulate_surface(
    atoms, probe_profile, sampling=None, shape=None, power_law=1.4
):
    """
    Simulate a STEM image of a surface using the convolution approximation.

    Args:

      atoms : Atoms object
      sampling : The shape of the output image.
      probe_profile : Function for calculating the probe profile.
      power_law : The assumed Z-contrast powerlaw

    Returns:
           ndarray
    """
    extent = np.diag(atoms.lattice_mat)[:2]

    if shape is None:
        if isinstance(sampling, Number):
            sampling = (sampling,) * 2
        shape = (
            int(np.ceil(extent[0] / sampling[0])),
            int(np.ceil(extent[1] / sampling[1])),
        )
    else:
        if isinstance(shape, Number):
            shape = (shape,) * 2

        sampling = (extent[0] / shape[0], extent[1] / shape[1])

    margin = int(np.ceil(5 / min(sampling)))
    shape_w_margin = (shape[0] + 2 * margin, shape[1] + 2 * margin)

    x = np.fft.fftfreq(shape_w_margin[0]) * shape_w_margin[1] * sampling[0]
    y = np.fft.fftfreq(shape_w_margin[1]) * shape_w_margin[1] * sampling[1]

    r = np.sqrt(x[:, None] ** 2 + y[None] ** 2)
    intensity = probe_profile(r)

    positions = atoms.cart_coords[:, :2] / sampling - 0.5

    inside = (
        (positions[:, 0] > -margin)
        & (positions[:, 1] > -margin)
        & (positions[:, 0] < shape[0] + margin)
        & (positions[:, 1] < shape[1] + margin)
    )

    positions = positions[inside] + margin
    numbers = np.array(atoms.atomic_numbers)[inside]

    array = np.zeros((1,) + shape_w_margin)
    for number in np.unique(np.array(atoms.atomic_numbers)):
        temp = np.zeros((1,) + shape_w_margin)
        superpose_deltas(positions[numbers == number], 0, temp)
        array += temp * number ** power_law

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
    c = crop_squre(a)
    c = a.make_supercell_matrix([2, 2, 1])
    p = simulate_surface(
        c,
        probe_profile=probe_profile(
            gaussian_width=1, lorentzian_width=1, intensity_ratio=1, n=100
        ),
        shape=[50, 50],
    )
    plt.imshow(p, interpolation="gaussian", cmap="plasma")
    plt.savefig("stem_example.png")
    plt.close()
"""
