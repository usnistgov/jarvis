"""Module to analyze phonon dos."""

import numpy as np
from phonopy.structure.atoms import isotope_data
from math import pi as pi

icm_to_eV = 1.23981e-4
icm_to_thz = 2.99792458e-2
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
        prefix = kB * x[1:] ** 2 * (np.exp(x[1:]) / (np.exp(x[1:]) - 1) ** 2)
        Cp = prefix * dos[1:]
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

    def phonon_isotope_scattering(self, atoms=None):
        """
        Get phonon-isotope scattering rate at natural isotopic abundance.
        Returns scattering rate in units of Hz.
        """
        omega = np.array(self.phonon_freq_cm)
        dos = np.array(self.phonon_dos)

        def isotopic_gamma(atoms):
            formula = atoms.composition.reduce()
            natoms = sum([v for v in formula[0].values()])
            ave_m = 0
            gamma = 0
            for k, v in formula[0].items():
                iso_list = isotope_data[k]
                ave_m_n = sum([iso[2] * iso[1] for iso in iso_list])
                g = [iso[2] * (iso[1] - ave_m_n) ** 2 for iso in iso_list]
                gamma_n = sum(g)
                ave_m += ave_m_n * (v / natoms)
                gamma += gamma_n * (v / natoms)
            return gamma / (ave_m ** 2)

        gamma = isotopic_gamma(atoms)
        atmV = (atoms.volume / atoms.num_atoms) * 1e-30
        omega = omega * icm_to_thz
        dos = dos / icm_to_thz / (atmV * atoms.num_atoms)
        tau = (pi / 6) * (atmV * gamma * omega ** 2) * dos
        return np.trapz(tau, omega) * 1e12


if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    from jarvis.db.figshare import get_jid_data

    dos_entry = get_jid_data(jid="JVASP-1459", dataset="edos_pdos")
    dft3d_entry = get_jid_data(jid="JVASP-1459", dataset="dft_3d")
    ph_dos = dos_entry["pdos_elast"]
    ph_freq = np.arange(0, 1000, 5)
    atoms = Atoms.from_dict(dft3d_entry["atoms"])

    ph = PhononDos(phonon_dos=ph_dos, phonon_freq_cm=ph_freq)
    debye_temp = ph.debye_temperature(atoms)
    iso_scatt = ph.phonon_isotope_scattering(atoms)

    print("Debye temperature:", debye_temp)
    print("Isotope scattering rate:", iso_scatt)
