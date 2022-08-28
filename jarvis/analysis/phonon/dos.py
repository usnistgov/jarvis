"""Module to analyze phonon dos."""

import numpy as np

icm_to_eV = 1.23981e-4
hbar = 6.582119569e-16  # eV*s
kB = 8.617333262145e-5  # eV/K
e = 1.60217662e-19
Na = 6.0221409e23


class PhononDos(object):
    """Module to analyze phonon dos."""

    def __init__(self, phonon_dos=[], phonon_freq_cm=[]):
        """Initialize class."""
        self.phonon_dos = phonon_dos
        self.phonon_freq_cm = phonon_freq_cm

    def debye_temperature(self, atoms=None):
        """Get Debye temperature."""
        # http://dx.doi.org/10.1103/PhysRevB.89.024304
        # Eq. 10
        n = atoms.num_atoms
        omega = np.array(self.phonon_freq_cm) * icm_to_eV
        gomega = np.array(self.phonon_dos)
        integ = np.trapz(omega ** 2 * gomega, omega) / np.trapz(gomega, omega)
        prefact = 1 / kB
        # TODO: check if np.pi / 2 is required
        moment_debye = (
            n ** (-1 / 3)
            * (prefact)
            * np.sqrt(5 / 3 * integ)
            # np.pi / 2 * n ** (-1 / 3) * (prefact) * np.sqrt(5 / 3 * integ)
        )
        return moment_debye

    def heat_capacity(self, temperature=300):
        """Get heat capacity at a temperature."""
        omega = np.array(self.phonon_freq_cm) * icm_to_eV
        # http://www.columbia.edu/~jh2228/phonons_thermal_hone.pdf
        # Eq. 1
        dos = np.array(self.phonon_dos) / icm_to_eV
        x = (omega) / (kB * temperature)
        Cp = (
            kB
            * x[1:] ** 2
            * (np.exp(x[1:]) / (np.exp(x[1:]) - 1) ** 2)
            * dos[1:]
        )
        Cp = np.insert(Cp, 0, 0)
        return np.trapz(Cp, omega) * e * Na

    def vibrational_entropy(self, temperature=300):
        """Get heat vibrational entropy at a temperature."""
        omega = np.array(self.phonon_freq_cm) * icm_to_eV
        dos = np.array(self.phonon_dos) / icm_to_eV
        x = (omega) / (kB * temperature)
        n = 1 / (np.exp(x[1:]) - 1)
        S_vib = kB * ((n + 1) * np.log(n + 1) + n * np.log(n)) * dos[1:]
        S_vib = np.insert(S_vib, 0, S_vib[0])
        return S_vib
