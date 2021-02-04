"""Module to process spectrum like data."""

import numpy as np
from jarvis.core.utils import lorentzian
from scipy.signal import find_peaks_cwt


class Spectrum(object):
    """Module for spectrum like data, e.g. IR, Raman, DOS, epsilon."""

    def __init__(self, x=[], y=[], linewidth=5.0, resolution=0.1):
        """Initialize the class."""
        self.x = np.array(x)
        self.y = np.array(y)
        self.resolution = resolution
        self.linewidth = linewidth

    def rescale(self, mode="max", scaling_factor=1.0):
        """Rescale the spectrum."""
        if mode == "sum":
            const = np.sum(self.y, axis=0)
        if mode == "max":
            const = np.max(self.y, axis=0)
        return self.y * scaling_factor / const

    @property
    def num_modes(self):
        """Get number of modes."""
        return len(self.x)

    @property
    def min_x(self):
        """Get minimum mode frequency."""
        return min(self.x)

    @property
    def max_x(self):
        """Get maximum mode frequency."""
        return max(self.x)

    def get_peak_indices(self, window=np.arange(1, 10)):
        """Get peak indices for non-zero peaks."""
        return find_peaks_cwt(self.y, window)

    def smoothen_spiky_spectrum(self):
        """Smoothen peak for delta function like peaks."""
        lwidth_list = [float(self.linewidth) for i in range(self.num_modes)]
        spect_x = np.arange(
            self.min_x,
            self.max_x + self.resolution,
            self.resolution,
            dtype=np.float64,
        )
        spect_y = np.zeros_like(spect_x, dtype=np.float64)
        for i, j, k in zip(self.x, self.y, lwidth_list):
            spect_y += lorentzian(spect_x, j, i, k)
        return spect_x, spect_y

    def get_interpolated_values(self, new_dist=np.arange(0, 15, 0.05)):
        """Get interpolated grid on a fixed grid."""
        interp = np.interp(new_dist, self.x, self.y)
        return interp


"""
from jarvis.io.vasp.outputs import Vasprun, Outcar
from jarvis.analysis.phonon.ir import ir_intensity

import os

out = Outcar(
    os.path.join(os.path.dirname(__file__),
     "../tests/testfiles/io/vasp/OUTCAR.JVASP-39")
)
vrun = Vasprun(
    os.path.join(
        os.path.dirname(__file__),
        "../tests/testfiles/io/vasp/vasprun.xml.JVASP-39"
    )
)


phonon_eigenvectors = vrun.dfpt_data["phonon_eigenvectors"]
vrun_eigs = vrun.dfpt_data["phonon_eigenvalues"]
phonon_eigenvalues = out.phonon_eigenvalues
masses = vrun.dfpt_data["masses"]
born_charges = vrun.dfpt_data["born_charges"]
x, y = ir_intensity(
    phonon_eigenvectors=phonon_eigenvectors,
    phonon_eigenvalues=phonon_eigenvalues,
    masses=masses,
    born_charges=born_charges,
)
assert x[0] == 713.8676686817399
sp=Spectrum(x=x,y=y)
xx,yy=sp.smoothen_spiky_spectrum()
print (sp.get_peak_indices())
#print (xx.tolist())
print ()
#print (yy.tolist())
"""
