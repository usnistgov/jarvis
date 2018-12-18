from __future__ import unicode_literals, print_function

# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp;
#               Kamal Choudhary, National Institute of Standards and Technology
# Distributed under the terms of the GNU License
# Initially forked and (extensively) adjusted from https://github.com/ldwillia/SL3ME

import os, math, cmath, pdb, json

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as constants

from scipy.interpolate import interp1d
from scipy.integrate import simps
from scipy.constants import physical_constants, speed_of_light
from math import pi
from monty.json import MSONable, MontyDecoder, MontyEncoder
from monty.io import zopen
from fnmatch import fnmatch
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.electronic_structure.core import Spin
from xml.etree.ElementTree import ParseError

"""
Module for calculating the SLME metric, including several classes for representing 
the optical properties of the material or an electromagnetic absorption spectrum.
"""

# Defining constants for tidy equations
c = constants.c  # speed of light, m/s
h = constants.h  # Planck's constant J*s (W)
h_e = constants.h / constants.e  # Planck's constant eV*s
k = constants.k  # Boltzmann's constant J/K
k_e = constants.k / constants.e  # Boltzmann's constant eV/K
e = constants.e  # Coulomb


class DielTensor(MSONable):
    """
    Class that represents the energy-dependent dielectric tensor of a solid
    state material.

    """

    def __init__(self, dielectric_data):
        """
        Initializes a DielTensor instance from the dielectric data.

        Args:
            dielectric_data (tuple): tuple of length N lists of dielectric data.

                - vasprun: The first tuple element contains the energies, the
                    second/third the real/imaginary part of the dielectric tensor.
                    Each dielectric entry should be a list of ``[xx, yy, zz, xy, xz,
                    yz ]`` dielectric tensor elements.
                - outcar: The first tuple element
        """
        self._energies, self._dielectric_tensor = self.parse_dielectric_data(
            dielectric_data
        )

    def parse_dielectric_data(self, dielectric_data):
        """
        Convert the set of - currently exclusively - vasprun formatted dielectric data to
        the corresponding 3x3 symmetric numpy matrices.

        Returns:
            (np.array), (np.array):  a Nx1 and Nx3 numpy array. The first array contains the
                energies corresponding to the dielectric tensor in the second.
        """

        # Check if the provided dielectric data is in the Vasprun tuple format (len == 3)
        if type(dielectric_data) is tuple and len(dielectric_data) == 3:
            energies = np.array(dielectric_data[0])

            dielectric_tensor = np.array(
                [to_matrix(*real_data) + 1j * to_matrix(*imag_data)
                 for real_data, imag_data in zip(dielectric_data[1],
                                                 dielectric_data[2])]
            )

            return energies, dielectric_tensor

        # Check if the provided dielectric data is in the Outcar tuple format (len == 2)
        elif type(dielectric_data) is tuple and len(dielectric_data) == 2:
            energies = dielectric_data[0]
            dielectric_tensor = dielectric_data[1]
            return energies, dielectric_tensor

        else:
            raise ImportError("Format of dielectric data not recognized.")

    def check_dielectric_data(self, dielectric_data):
        """
        Function that performs some tests on the dielectric data, to make sure
        input satisfies some constrains based on what we know about the dielectric
        tensor.

        Returns:
            None

        """
        pass  # TODO

    @property
    def energies(self):
        """
        Energy grid for which the dielectric tensor is defined in the original data.

        Returns:
            (numpy.array): (N,) shaped array with the energies of the grid in eV.

        """
        return self._energies

    @property
    def dielectric_tensor(self):
        """
        Dielectric tensor of the material, calculated for each energy in the energy grid.

        Returns:
            (numpy.array): (N, 3, 3) shaped array, where N corresponds to the number
                of energy
            grid points, and 3 to the different directions x,y,z.

        """
        return self._dielectric_tensor

    @property
    def dielectric_function(self):
        """
        The averaged dielectric function, derived from the tensor components by
        first diagonalizing the dielectric tensor for every energy and then
        averaging the diagonal elements.

        Returns:
            (np.array):
        """
        return np.array([np.mean(np.linalg.eigvals(tensor))
                         for tensor in self.dielectric_tensor])

    @property
    def absorption_coefficient(self):
        """
        Calculate the optical absorption coefficient from the dielectric data.
        For now the script only calculates the averaged absorption coefficient,
        i.e. by first averaging the diagonal elements and then using this
        dielectric function to calculate the absorption coefficient.

        Notes:
            The absorption coefficient is calculated as
            .. math:: \\alpha = \\frac{2 E}{ \hbar c} k(E)
            with $k(E)$ the imaginary part of the square root of the dielectric function

        Returns:
            (np.array): Energy (eV) dependent absorption coefficient in m^{-1}, where the
                energies correspond to self.energies.
        """

        energy = self.energies
        ext_coeff = np.array([cmath.sqrt(v).imag for v in self.dielectric_function])

        return 2.0 * energy * ext_coeff / (constants.hbar / constants.e * constants.c)

    def get_absorptivity(self, thickness, method="beer-lambert"):
        """
        Calculate the absorptivity for an absorber layer with a specified thickness
        and cell construction.

        Args:
            thickness (float): Thickness of the absorber layer.
            method (str): Method for calculating the absorptivity.

        Returns:

        """
        if method == "beer-lambert":
            return 1.0 - np.exp(-2.0 * self.absorption_coefficient * thickness)

    def plot(self, part="diel", variable_range=None, diel_range=None):
        """
        Plot the real and/or imaginary part of the dielectric function.

        Args:
            part (str): Which part of the dielectric function to plot, i.e. either "real",
            "imag" or "all".

        Returns:
            None

        """
        if part == "diel":
            f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

            ax1.plot(self.energies, self.dielectric_function.real)
            ax2.plot(self.energies, self.dielectric_function.imag)
            if variable_range:
                ax1.set(xlim=variable_range)
                ax2.set(xlim=variable_range)
            if diel_range:
                ax1.set(ylim=diel_range)
                ax2.set(ylim=diel_range)
            ax1.set(ylabel=r"$\varepsilon_1$")
            ax2.set(xlabel="Energy (eV)", ylabel=r"$\varepsilon_2$")
            f.subplots_adjust(hspace=0.1)
            plt.show()

        elif part == "real":

            plt.plot(self.energies, self.dielectric_function.real)
            plt.xlabel("Energy (eV)")
            if variable_range:
                plt.xlim(variable_range)
            if diel_range:
                plt.ylim(diel_range)
            plt.ylabel(r"$\varepsilon_1$")
            plt.show()

        elif part == "imag":

            plt.plot(self.energies, self.dielectric_function.imag)
            plt.xlabel("Energy (eV)")
            if variable_range:
                plt.xlim(variable_range)
            if diel_range:
                plt.ylim(diel_range)
            plt.ylabel(r"$\varepsilon_2$")
            plt.show()

        elif part == "abs_coeff":

            plt.plot(self.energies, self.absorption_coefficient)
            plt.xlabel("Energy (eV)")
            if variable_range:
                plt.xlim(variable_range)
            if diel_range:
                plt.ylim(diel_range)
            plt.ylabel(r"$\alpha(E)$")
            plt.yscale("log")
            plt.show()

    def as_dict(self):
        """
        Note: stores the real and imaginary part of the dielectric tensor
        separately, due to issues with JSON serializing complex numbers.

        Returns:
            dict:
        """
        d = dict()
        d["energies"] = self.energies
        d["real_diel"] = self.dielectric_tensor.real
        d["imag_diel"] = self.dielectric_tensor.imag
        return d

    @classmethod
    def from_dict(cls, d):
        energies = np.array(d["energies"]["data"])
        real_diel = np.array(d["real_diel"]["data"])
        imag_diel = np.array(d["imag_diel"]["data"])
        return cls((energies, real_diel + 1j*imag_diel))

    def to(self, filename):
        """
        Write the DielTensor to a JSON file.

        Args:
            filename (str): Path to the file in which the DielTensor should
                be written.

        Returns:
            None

        """
        with zopen(filename, "w") as f:
            f.write(self.to_json())

    @classmethod
    def from_file(cls, filename, format=None):
        """
        Initialize a DielTensor instance from a file.

        Args:
            filename (str): Path to file from which the dielectric data will be
                loaded. Can (so far) either be a vasprun.xml, OUTCAR or json file.
            fmt (str): Format of the file that contains the dielectric function
                data. Is optional, as the method can also figure out the format
                based on the filename.

        Returns:
            (DielTensor): Dielectric tensor object from the dielectric data.

        """
        # Vasprun format: dielectric data is length 3 tuple
        if format == "vasprun" or filename.endswith(".xml"):
            dielectric_data = Vasprun(filename, parse_potcar_file=False).dielectric
            return cls(dielectric_data)

        # OUTCAR format: dielectric data is length 2 tuple
        elif format == "outcar" or fnmatch(filename, "*OUTCAR*"):
            outcar = Outcar(filename)
            outcar.read_freq_dielectric()
            dielectric_data = (outcar.frequencies,
                               outcar.dielectric_tensor_function)
            return cls(dielectric_data)

        # JSON format
        if format == "json" or filename.endswith(".json"):
            with zopen(filename, "r") as f:
                return cls.from_dict(json.loads(f.read()))

        else:
            raise IOError("Format of file not recognized. Note: Currently "
                          "only vasprun.xml and OUTCAR files are supported.")


class EMRadSpectrum(MSONable):
    """
    Class that represents a electromagnetic radiation spectrum, e.g. the solar spectrum
    or the black body spectrum. The standard form we choose to express the spectrum is
    as the energy-dependent photon flux per square meter, where the energy is expressed
    in electronvolts (Units ~ m^{-2} s^{-1} eV^{-1}).

    """

    def __init__(self, energy, photon_flux):
        """
        Initialize the Radiation Spectrum object from the energy grid and the
        photon flux. The input is expected to be provided in the units described
        below.

        Args:
            energy (numpy.array): Energy grid for which the photon flux is given.
                Has to be given in electronvolt (eV).
            photon_flux (numpy.array): Number of photons per square meter per
                second per electronvolt. (~ m^{-2} s^{-1} eV^{-1})

        """
        self._energy = energy
        self._photon_flux = photon_flux

    @property
    def energy(self):
        """
        Energy grid for which the electromagnetic radiation spectrum is defined.

        Returns:
            (numpy.array)

        """
        return self._energy

    @property
    def photon_flux(self):
        """
        Electromagnetic radiation spectrum expressed in photon flux.

        Returns:
            (numpy.array)

        """
        return self._photon_flux

    def get_total_power_density(self):
        """
        Get the total power density in W m^{-2}.

        Returns:
            (float)

        """
        return simps(self.photon_flux * self.energy * e, self.energy)

    def get_interp_function(self, variable="energy", spectrum_units="flux"):
        """
        Obtain the 1D interpolation function using the scipy.interpolate.interp1d
        method. Linear interpolation is chosen, as it is more robust.

        Still have to implement the option to change the variable and type of
        spectrum.

        Args:
            variable (str):
            spectrum_units (str):

        Returns:

        """
        # TODO complete for other variables and spectrum choices
        return interp1d(self._energy, self._photon_flux, kind='linear',
                        fill_value=0.0, bounds_error=False)

    def to(self, filename):
        """
        Write the EMRadSpectrum to a JSON file.

        Args:
            filename (str): Path to the file in which the EMRadSpectrum should
                be written.

        Returns:
            None

        """
        with zopen(filename, "w") as f:
            f.write(self.to_json())

    @classmethod
    def from_file(cls, filename):
        """
        Load the EMRadSpectrum from a file.

        Args:
            filename (str): Path to the file from which the EMRadSpectrum should
                be loaded.

        Returns:
            (EMRadSpectrum)

        """
        if filename.endswith(".json"):
            with zopen(filename, "r") as f:
                return cls.from_dict(json.loads(f.read()))
        else:
            raise IOError("Filename does not have .json extension.")

    @classmethod
    def from_data(cls, data, variable="energy", spectrum_units="flux"):
        """

        Args:
            data (tuple): Tuple with length 2 that contains the data from which to
                construct the radiation spectrum. data[0] must contain a (N,) shaped
                numpy.array with the grid of the variable (e.g. energy, wavelength),
                data[1] should contain the spectral distribution.
            variable (str):
            spectrum_units (str):

        Returns:

        """
        if variable == "energy":
            energy = data[0]
            spectrum = data[1]

        elif variable == "wavelength":
            energy = np.flip(h_e * c / data[0])
            spectrum = np.flip(data[1]) * h_e * c / energy ** 2

        else:
            raise NotImplementedError

        if spectrum_units == "flux":
            photon_flux = spectrum
            return cls(energy, photon_flux)

        elif spectrum_units == "irradiance":
            photon_flux = spectrum / (e * energy)
            return cls(energy, photon_flux)

        else:
            raise NotImplementedError

    @classmethod
    def get_solar_spectrum(cls, spectrum="am1.5g"):
        """

        Args:
            spectrum:

        Returns:

        """

        data_file = os.path.join(os.path.dirname(__file__), spectrum + ".dat")

        wavelength, irradiance = np.loadtxt(
            data_file, usecols=[0, 1], unpack=True, skiprows=2
        )

        # Transfer units to m instead of nm
        wavelength *= 1e-9
        irradiance *= 1e9

        return cls.from_data((wavelength, irradiance), variable="wavelength",
                             spectrum_units="irradiance")

    @classmethod
    def get_blackbody(cls, temperature, grid, variable="energy",
                      spectrum_units="flux"):
        """
        Construct the blackbody spectrum of a specific temperature.

        Args:
            temperature (float): Temperature of the black body.
            grid (numpy.array): Grid of the spectral variable which the black
                body is dependent on.
            variable (str): Spectral variable of the distribution.
            spectrum_units (str): Units in which the spectrum should be expressed.
                Currently only support expressing the black body as a photon flux.

        Returns:
            (EMRadSpectrum)

        """
        if variable == "energy":
            energy = grid
        else:
            raise NotImplementedError

        # Define a exponential function that does not cause a range overflow
        def exponential(x):
            try:
                return math.exp(x)
            except OverflowError:
                return math.exp(700)  # ~= 1e304

        # Calculation of energy-dependent blackbody spectrum (~ W m^{-2})
        if spectrum_units == "flux":
            photon_flux = 2 * energy ** 2 / (h_e ** 3 * c ** 2) * (
                    1 / (np.array([exponential(energy / (k_e * temperature)) - 1
                                   for energy in energy]))
            )
        else:
            raise NotImplementedError

        return cls(energy, photon_flux)


class SolarCell(MSONable):
    """
    Class that represents a single p-n junction solar cell. Contains several modeling
    techniques for the calculation of the theoretical efficiency using metrics
    based on the optical properties and electronic structure (e.g. band gap)
    of the material being considered as the absorber layer.

    """

    def __init__(self, dieltensor, bandgaps):
        """
        Initialize an instance of the EfficiencyCalculator class.

        Args:
            dieltensor (DielTensor): Dielectric tensor of the absorber material.
            bandgaps (tuple): Tuple that contains the fundamental and direct
                allowed band gap of the absorber material, in that order.

        Returns:
            (SolarCell)

        """
        self._dieltensor = dieltensor
        self._bandgaps = bandgaps

    @property
    def dieltensor(self):
        return self._dieltensor

    @property
    def bandgaps(self):
        return self._bandgaps

    def slme(self, temperature=298.15, thickness=5e-7, interp_mesh=0.001,
             plot_iv_curve=False, cut_abs_below_bandgap=False):
        """
        Calculate the Spectroscopic Limited Maximum Efficiency.

        Args:
            temperature (float): Temperature of the solar cell.
            thickness (float): Thickness of the absorber layer.
            interp_mesh (float): Distance between two energy points in the grid
                used for the interpolation.
            plot_iv_curve (bool): Defaults to False. If set to true,

        Returns:
            (tuple) efficiency, v_oc, j_sc, j_0

        """
        # Set up the energy grid for the calculation
        energy = self.dieltensor.energies
        energy = np.linspace(
            np.min(energy) + interp_mesh, np.max(energy),
            np.ceil((np.max(energy) - np.min(energy)) / interp_mesh)
        )

        # Interpolation of the absorptivity to the new energy grid
        absorptivity = interp1d(
            self._dieltensor.energies,
            self._dieltensor.get_absorptivity(thickness, "beer-lambert"),
            kind='linear',
            fill_value=0,
            bounds_error=False
        )(energy)

        # If the user has requested the onset below the band gap to be removed
        if cut_abs_below_bandgap:
            # Set the absorptivity to zero for energies below the *direct* band gap
            cut_array = np.array([int(el) for el in energy > self.bandgaps[1]])
            absorptivity *= cut_array

        # Load energy-dependent total solar spectrum photon flux (~m^{-2}s^{-1}eV^{-1})
        solar_spectrum = \
            EMRadSpectrum.get_solar_spectrum("am1.5g").get_interp_function()(energy)

        # Calculation of energy-dependent blackbody spectrum (~m^{-2}s^{-1}eV^{-1})
        blackbody_spectrum = EMRadSpectrum.get_blackbody(temperature, energy).photon_flux

        # Numerically integrating photon flux over energy grid
        j_0_r = e * np.pi * simps(blackbody_spectrum * absorptivity, energy)

        # Calculate the fraction of radiative recombination
        delta = self._bandgaps[1] - self._bandgaps[0]
        fr = np.exp(-delta / (k_e * temperature))
        j_0 = j_0_r / fr

        # Numerically integrating irradiance over wavelength array ~ A/m**2
        j_sc = e * simps(solar_spectrum * absorptivity, energy)

        # Determine the open circuit voltage.
        v_oc = 0
        voltage_step = 0.001
        while j_sc - j_0 * (np.exp(e * v_oc / (k * temperature)) - 1.0) > 0:
            v_oc += voltage_step

        if plot_iv_curve:  # TODO Add some more details.
            voltage = np.linspace(0, v_oc, 2000)

            current = j_sc - j_0 * (np.exp(e * voltage / (k * temperature)) - 1.0)
            power = current * voltage

            plt.plot(voltage, current)
            plt.plot(voltage, power)
            plt.show()

        # Maximize the power density versus the voltage
        max_power = self.maximize_power(j_sc, j_0, temperature)

        # Calculation of integrated solar spectrum
        power_in = EMRadSpectrum.get_solar_spectrum().get_total_power_density()

        # Calculate the maximized efficiency
        efficiency = max_power / power_in

        return efficiency, v_oc, j_sc, j_0

    def plot_slme_vs_thickness(self, temperature=298.15, add_sq_limit=True,
                               cut_abs_below_bandgap=False):
        """
        Make a plot of the calculated SLME for a large range of thickness values,
        for a specific temperature.

        Args:
            temperature (float): Temperature of the solar cell. Defaults to 298.15 K.
            add_sq_limit (bool): Specifies whether the user would like to add a
                line representing the Shockley-Queisser limit that corresponds to
                the band gap of the absorber material at the temperature specified.

        Returns:
            None

        """
        thickness = 10 ** np.linspace(-9, -3, 40)
        efficiency = np.array([self.slme(thickness=d, temperature=temperature,
                                         cut_abs_below_bandgap=cut_abs_below_bandgap)[0]
                               for d in thickness])

        plt.plot(thickness, efficiency)
        if add_sq_limit:
            plt.plot(thickness, np.ones(thickness.shape)
                     * self.calculate_bandgap_sq(temperature=temperature)[0], 'k--')
            plt.legend(("SLME", "SQ"))
        else:
            plt.legend(("SLME",))
        plt.xlabel("Thicknesss (m)")
        plt.ylabel("Efficiency")
        plt.xscale("log")

        plt.show()

    def get_currents(self, temperature):
        pass

    def get_iv_curve(self, j_sc, j_0, temperature):
        pass

    @classmethod
    def from_file(cls, filename):
        """
        Loads a SolarCell instance from a vasprun.xml. # TODO extend

        Returns:
            (SolarCell)

        """
        try:
            vasprun = Vasprun(filename, parse_potcar_file=False)
        except ParseError:
            raise IOError("Error while parsing the input file. Currently the "
                          "SolarCell class can only be constructed from "
                          "the vasprun.xml file. If you have provided this "
                          "file, check if the run has completed.")

        diel_tensor = DielTensor(vasprun.dielectric)

        # Extract the information on the direct and indirect band gap
        bandstructure = vasprun.get_band_structure()
        bandgaps = (bandstructure.get_band_gap()["energy"],
                    bandstructure.get_direct_band_gap())

        return cls(diel_tensor, bandgaps)

    @staticmethod
    def maximize_power(j_sc, j_0, temperature):
        """
        Maximize the power density based on the short-circuit current j_sc, the
        recombination current j_0 and the temperature of the solar cell.

        Args:
            j_sc (float): Short-Circuit current density.
            j_0 (float): Recombination current density.
            temperature (float): Temperature of the solar cell.

        Returns:
            (float): The calculated maximum power.

        """

        # Calculate the current density J for a specified voltage V
        def current_density(voltage):
            j = j_sc - j_0 * (np.exp(e * voltage / (k * temperature)) - 1.0)
            return j

        # Calculate the corresponding power density P
        def power(voltage):
            p = current_density(voltage) * voltage
            return p

        # A somewhat primitive, but perfectly robust way of getting a reasonable
        # estimate for the maximum power.
        test_voltage = 0
        voltage_step = 0.001
        while power(test_voltage + voltage_step) > power(test_voltage):
            test_voltage += voltage_step

        return power(test_voltage)

    @staticmethod
    def calculate_slme_from_vasprun(filename, temperature=298.15, thickness=5e-7):
        return SolarCell.from_file(filename).slme(temperature, thickness)

    def calculate_bandgap_sq(self, temperature=298.15, fr=1.0, interp_mesh=0.001):
        """
        Calculate the Shockley-Queisser limit of the corresponding fundamental
        band gap.

        Args:
            temperature (float): Temperature of the solar cell. Defaults to 25 °C,
                 or 298.15 K.
            fr (float): Fraction of radiative recombination.
            interp_mesh (float): Distance between two energy points in the grid
                used for the interpolation.

        Returns:
            (float): Shockley-Queisser detailed balance limit of the band gap of
                the material.

        """
        return self.sq(self._bandgaps[0], temperature, fr, interp_mesh,
                       np.max(self.dieltensor.energies))

    @staticmethod
    def sq(bandgap, temperature=298.15, fr=1.0, interp_mesh=0.001, max_energy=20.0):
        """
        Calculate the Shockley-Queisser limit for a specified bandgap, temperature
        and fraction of radiative recombination.

        Args:
            bandgap (float): Band gap of the absorber material in eV.
            temperature (float): Temperature of the solar cell. Defaults to 25 °C,
                 or 298.15 K.
            fr (float): Fraction of radiative recombination.
            interp_mesh (float): Distance between two energy points in the grid
                used for the interpolation.
            max_energy (float): Maximum energy in the energy grid.

        Returns:
            (float): Shockley-Queisser detailed balance limit.

        """
        # Set up the energy grid for the calculation
        energy = np.linspace(
            interp_mesh, max_energy, np.ceil(max_energy) / interp_mesh
        )

        # Set up the absorption coefficient (Step function for SQ)
        absorptivity = np.array(
            [float(ener > bandgap) for ener in energy]
        )

        # Get total solar_spectrum
        solar_spectrum = \
            EMRadSpectrum.get_solar_spectrum("am1.5g").get_interp_function()(energy)

        # Calculation of energy-dependent blackbody spectrum, in units of W / m**2
        blackbody_spectrum = EMRadSpectrum.get_blackbody(temperature, energy).photon_flux

        # Numerically integrating irradiance over energy grid ~ A/m**2
        j_0_r = e * np.pi * simps(blackbody_spectrum * absorptivity, energy)
        j_0 = j_0_r / fr

        # Numerically integrating irradiance over wavelength array ~ A/m**2
        j_sc = e * simps(solar_spectrum * absorptivity, energy)

        # Maximize the power versus the voltage
        max_power = SolarCell.maximize_power(j_sc, j_0, temperature)

        # Calculation of integrated solar spectrum
        power_in = EMRadSpectrum.get_solar_spectrum().get_total_power_density()

        # Calculate the maximized efficiency
        efficiency = max_power / power_in

        return efficiency, j_sc, j_0


# Utility method
def to_matrix(xx, yy, zz, xy, yz, xz):
    """
    Convert a list of matrix components to a symmetric 3x3 matrix.
    Inputs should be in the order xx, yy, zz, xy, yz, xz.

    Args:
        xx (float): xx component of the matrix.
        yy (float): yy component of the matrix.
        zz (float): zz component of the matrix.
        xy (float): xy component of the matrix.
        yz (float): yz component of the matrix.
        xz (float): xz component of the matrix.

    Returns:
        (np.array): The matrix, as a 3x3 numpy array.

    """
    matrix = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
    return matrix


######################################
# Previous code, kept for comparison #
######################################

eV_to_recip_cm = 1.0 / (physical_constants['Planck constant in eV s'][0] * speed_of_light * 1e2)


def nelec_out(out=''):
    f = open(out, 'r')
    lines = f.read().splitlines()
    f.close()
    for i in lines:
        if 'NELECT =' in i:
            nelec = int(float(
                i.split()[2]))  # int(i.split('NELECT =')[1].split('total number of electrons')[0])
    return nelec


def get_dir_indir_gap(run='', ispin=1):
    """
    Get direct and indirect bandgaps
    Implemented for non-spin polarized case only right now
    """
    v = Vasprun(run)
    outcar = run.replace('vasprun.xml', 'OUTCAR')
    nbands = len(v.eigenvalues[Spin.up][1])
    nelec = nelec_out(outcar)  # ispin*len(v.eigenvalues.items()[0][1])
    nkpts = int(len(v.eigenvalues[Spin.up]))
    eigvals = np.zeros((ispin, nkpts, nbands))
    for kp, i in enumerate(v.eigenvalues[Spin.up]):
        for band, j in enumerate(i):
            eigvals[0, kp, band] = j[0]
    noso_homo = np.max(eigvals[0, :, nelec // 2 - 1])
    noso_lumo = np.min(eigvals[0, :, nelec // 2])
    noso_indir = np.min(eigvals[0, :, nelec // 2]) - np.max(eigvals[0, :, nelec // 2 - 1])
    noso_direct = np.min(eigvals[0, :, nelec // 2] - eigvals[0, :, nelec // 2 - 1])
    print('noso_direct,noso_indir', noso_direct, noso_indir)
    return noso_direct, noso_indir


def matrix_eigvals(matrix):
    """
    Calculate the eigenvalues of a matrix.
    Args:
        matrix (np.array): The matrix to diagonalise.
    Returns:
        (np.array): Array of the matrix eigenvalues.
    """
    eigvals, eigvecs = np.linalg.eig(matrix)
    return eigvals


def parse_dielectric_data(data):
    """
    Convert a set of 2D vasprun formatted dielectric data to
    the eigenvalues of each corresponding 3x3 symmetric numpy matrices.
    Args:
        data (list): length N list of dielectric data. Each entry should be
                     a list of ``[xx, yy, zz, xy, xz, yz ]`` dielectric
                     tensor elements.
    Returns:
        (np.array):  a Nx3 numpy array. Each row contains the eigenvalues
                     for the corresponding row in `data`.
    """
    return np.array([matrix_eigvals(to_matrix(*e)) for e in data])


def absorption_coefficient(dielectric):
    """
    Calculate the optical absorption coefficient from an input set of
    pymatgen vasprun dielectric constant data.
    Args:
        dielectric (list): A list containing the dielectric response function
                           in the pymatgen vasprun format.
                           | element 0: list of energies
                           | element 1: real dielectric tensors, in ``[xx, yy, zz, xy, xz, yz]`` format.
                           | element 2: imaginary dielectric tensors, in ``[xx, yy, zz, xy, xz, yz]`` format.
    Returns:
        (np.array): absorption coefficient using eV as frequency units (cm^-1).
    Notes:
        The absorption coefficient is calculated as
        .. math:: \\alpha = \\frac{2\sqrt{2} \pi}{\lambda} \sqrt{-\epsilon_1+\sqrt{\epsilon_1^2+\epsilon_2^2}}
    """
    energies_in_eV = np.array(dielectric[0])
    real_dielectric = parse_dielectric_data(dielectric[1])
    imag_dielectric = parse_dielectric_data(dielectric[2])
    epsilon_1 = np.mean(real_dielectric, axis=1)
    epsilon_2 = np.mean(imag_dielectric, axis=1)
    return energies_in_eV, (2.0 * np.sqrt(2.0) * pi * eV_to_recip_cm * energies_in_eV
                            * np.sqrt(-epsilon_1 + np.sqrt(epsilon_1 ** 2 + epsilon_2 ** 2)))


def optics(ru=''):
    dirgap, indirgap = get_dir_indir_gap(ru)

    run = Vasprun(ru, occu_tol=1e-2, )
    new_en, new_abs = absorption_coefficient(run.dielectric)
    return np.array(new_en, dtype=np.float64), np.array(new_abs, dtype=np.float64), dirgap, indirgap


def calculate_SQ(bandgap_ev, temperature=300, fr=1,
                 plot_current_voltage=False):
    """
    Args:
        bandgap_ev: bandga in electron-volt
        temperature: temperature in K
    Returns:

    """

    # Defining constants for tidy equations
    c = constants.c  # speed of light, m/s
    h = constants.h  # Planck's constant J*s (W)
    h_e = constants.h / constants.e  # Planck's constant eV*s
    k = constants.k  # Boltzmann's constant J/K
    k_e = constants.k / constants.e  # Boltzmann's constant eV/K
    e = constants.e  # Coulomb

    # Load the Air Mass 1.5 Global tilt solar spectrum
    solar_spectrum_data_file = str(os.path.join(os.path.dirname(__file__), "am1.5g.dat"))
    solar_spectra_wavelength, solar_spectra_irradiance = np.loadtxt(
        solar_spectrum_data_file, usecols=[0, 1], unpack=True, skiprows=2
    )

    solar_spectra_wavelength_meters = solar_spectra_wavelength * 1e-9

    # need to convert solar irradiance from Power/m**2(nm) into
    # photon#/s*m**2(nm) power is Watt, which is Joule / s
    # E = hc/wavelength
    # at each wavelength, Power * (wavelength(m)/(h(Js)*c(m/s))) = ph#/s
    solar_spectra_photon_flux = solar_spectra_irradiance * (
            solar_spectra_wavelength_meters / (h * c))

    ### Calculation of total solar power incoming
    power_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    # calculation of blackbody irradiance spectra
    # units of W/(m**3), different than solar_spectra_irradiance!!! (This
    # is intentional, it is for convenience)
    blackbody_irradiance = (2.0 * h * c ** 2 /
                            (solar_spectra_wavelength_meters ** 5)) \
                           * (1.0 / ((np.exp(h * c / (
            solar_spectra_wavelength_meters * k * temperature))) - 1.0))

    # I've removed a pi in the equation above - Marnik Bercx

    # now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (
            solar_spectra_wavelength_meters / (h * c))

    # absorbance interpolation onto each solar spectrum wavelength
    from numpy import interp

    # Get the bandgap in wavelength in meters
    bandgap_wavelength = h_e * c / bandgap_ev

    # Only take the part of the wavelength-dependent solar spectrum and
    # blackbody spectrum below the bandgap wavelength
    bandgap_index = np.searchsorted(solar_spectra_wavelength_meters, bandgap_wavelength)

    bandgap_irradiance = interp(
        np.array([bandgap_wavelength, ]), solar_spectra_wavelength_meters,
        solar_spectra_photon_flux)

    bandgap_blackbody = (2.0 * h * c ** 2 /
                         (bandgap_wavelength ** 5)) \
                        * (1.0 / ((np.exp(h * c / (
            bandgap_wavelength * k * temperature))) - 1.0)) * (
                                bandgap_wavelength / (h * c))

    integration_wavelength = np.concatenate(
        (solar_spectra_wavelength_meters[:bandgap_index],
         np.array([bandgap_wavelength, ])),
        axis=0
    )

    integration_solar_flux = np.concatenate(
        (solar_spectra_photon_flux[:bandgap_index],
         bandgap_irradiance),
        axis=0
    )

    integration_blackbody = np.concatenate(
        (blackbody_photon_flux[:bandgap_index], np.array([bandgap_blackbody])),
        axis=0
    )

    #  Numerically integrating irradiance over wavelength array
    # Note: elementary charge, not math e!  ## units of A/m**2   W/(V*m**2)
    J_0_r = e * np.pi * simps(integration_blackbody,
                              integration_wavelength)

    J_0 = J_0_r / fr

    #  Numerically integrating irradiance over wavelength array
    # elementary charge, not math e!  ### units of A/m**2   W/(V*m**2)
    J_sc = e * simps(integration_solar_flux * 1e9, integration_wavelength)

    #    J[i] = J_sc - J_0*(1 - exp( e*V[i]/(k*T) ) )
    #   #This formula from the original paper has a typo!!
    #    J[i] = J_sc - J_0*(exp( e*V[i]/(k*T) ) - 1)
    #   #Bercx chapter and papers have the correct formula (well,
    #   the correction on one paper)
    def J(V):
        J = J_sc - J_0 * (np.exp(e * V / (k * temperature)) - 1.0)
        return J

    def power(V):
        p = J(V) * V
        return p

    # A more primitive, but perfectly robust way of getting a reasonable
    # estimate for the maximum power.
    test_voltage = 0
    voltage_step = 0.001
    while power(test_voltage + voltage_step) > power(test_voltage):
        test_voltage += voltage_step

    max_power = power(test_voltage)

    # Calculate the maximized efficience
    efficiency = max_power / power_in

    # This segment isn't needed for functionality at all, but can display a
    # plot showing how the maximization of power by choosing the optimal
    # voltage value works
    if plot_current_voltage:
        V = np.linspace(0, 2, 200)
        plt.plot(V, J(V))
        plt.plot(V, power(V), linestyle='--')
        plt.show()
        print(max_power)

    return efficiency


def slme(material_energy_for_absorbance_data,
         material_absorbance_data,
         material_direct_allowed_gap,
         material_indirect_gap, thickness=50E-6, temperature=293.15,
         absorbance_in_inverse_centimeters=False,
         cut_off_absorbance_below_direct_allowed_gap=True,
         plot_current_voltage=False):
    """
    Calculate the
    IMPORTANT NOTES:
    1) Material calculated absorbance is assumed to be in m^-1, not cm^-1!
        (Most sources will provide absorbance in cm^-1, so be careful.)
    2) The default is to remove absorbance below the direct allowed gap.
        This is for dealing with broadening applied in DFT absorbance
        calculations. Probably not desired for experimental data.
    3) We can calculate at different temperatures if we want to, but 25 C /
        293.15 K is the standard temperature assumed if not specified
    4) If absorbance is in cm^-1, multiply values by 100 to match units
        assumed in code
    Args:
        material_energy_for_absorbance_data:
        material_absorbance_data:
        material_direct_allowed_gap:
        material_indirect_gap:
        thickness:
        temperature:
        absorbance_in_inverse_centimeters:
        cut_off_absorbance_below_direct_allowed_gap:
        plot_current_voltage:
    Returns:
        The calculated maximum efficiency.
    """

    # Defining constants for tidy equations
    c = constants.c  # speed of light, m/s
    h = constants.h  # Planck's constant J*s (W)
    h_e = constants.h / constants.e  # Planck's constant eV*s
    k = constants.k  # Boltzmann's constant J/K
    k_e = constants.k / constants.e  # Boltzmann's constant eV/K
    e = constants.e  # Coulomb

    # Make sure the absorption coefficient has the right units (m^{-1})
    if absorbance_in_inverse_centimeters:
        material_absorbance_data = material_absorbance_data * 100

    # Load the Air Mass 1.5 Global tilt solar spectrum
    solar_spectrum_data_file = str(os.path.join(os.path.dirname(__file__), "am1.5g.dat"))

    solar_spectra_wavelength, solar_spectra_irradiance = np.loadtxt(
        solar_spectrum_data_file, usecols=[0, 1], unpack=True, skiprows=2
    )

    solar_spectra_wavelength_meters = solar_spectra_wavelength * 1e-9

    delta = material_direct_allowed_gap - material_indirect_gap
    fr = np.exp(-delta / (k_e * temperature))

    # need to convert solar irradiance from Power/m**2(nm) into
    # photon#/s*m**2(nm) power is Watt, which is Joule / s
    # E = hc/wavelength
    # at each wavelength, Power * (wavelength(m)/(h(Js)*c(m/s))) = ph#/s
    solar_spectra_photon_flux = solar_spectra_irradiance * (
            solar_spectra_wavelength_meters / (h * c))

    # Calculation of total solar power incoming
    power_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    # calculation of blackbody irradiance spectra
    # units of W/(m**3), different than solar_spectra_irradiance!!! (This
    # is intentional, it is for convenience)
    blackbody_irradiance = (2.0 * h * c ** 2 /
                            (solar_spectra_wavelength_meters ** 5)) \
                           * (1.0 / ((np.exp(h * c / (
            solar_spectra_wavelength_meters * k * temperature))) - 1.0))

    # I've removed a pi in the equation above - Marnik Bercx

    # now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (
            solar_spectra_wavelength_meters / (h * c))

    # units of nm
    material_wavelength_for_absorbance_data = ((c * h_e) / (
            material_energy_for_absorbance_data + 0.00000001)) * 10 ** 9

    # absorbance interpolation onto each solar spectrum wavelength
    from scipy.interpolate import interp1d
    # creates cubic spline interpolating function, set up to use end values
    #  as the guesses if leaving the region where data exists
    material_absorbance_data_function = interp1d(
        material_wavelength_for_absorbance_data, material_absorbance_data,
        kind='cubic',
        fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
        bounds_error=False
    )

    material_interpolated_absorbance = np.zeros(
        len(solar_spectra_wavelength_meters))
    for i in range(0, len(solar_spectra_wavelength_meters)):
        ## Cutting off absorption data below the gap. This is done to deal
        # with VASPs broadening of the calculated absorption data

        if solar_spectra_wavelength[i] < 1e9 * ((c * h_e) /
                                                material_direct_allowed_gap) \
                or cut_off_absorbance_below_direct_allowed_gap == False:
            material_interpolated_absorbance[
                i] = material_absorbance_data_function(
                solar_spectra_wavelength[i])

    absorbed_by_wavelength = 1.0 - np.exp(-2.0 *
                                          material_interpolated_absorbance
                                          * thickness)

    #  Numerically integrating irradiance over wavelength array
    # Note: elementary charge, not math e!  ## units of A/m**2   W/(V*m**2)
    J_0_r = e * np.pi * simps(blackbody_photon_flux * absorbed_by_wavelength,
                              solar_spectra_wavelength_meters)

    J_0 = J_0_r / fr

    #  Numerically integrating irradiance over wavelength array
    # elementary charge, not math e!  ### units of A/m**2   W/(V*m**2)
    J_sc = e * simps(solar_spectra_photon_flux * absorbed_by_wavelength,
                     solar_spectra_wavelength)

    #    J[i] = J_sc - J_0*(1 - exp( e*V[i]/(k*T) ) )
    #   #This formula from the original paper has a typo!!
    #    J[i] = J_sc - J_0*(exp( e*V[i]/(k*T) ) - 1)
    #   #Bercx chapter and papers have the correct formula (well,
    #   the correction on one paper)
    def J(V):
        J = J_sc - J_0 * (np.exp(e * V / (k * temperature)) - 1.0)
        return J

    def power(V):
        p = J(V) * V
        return p

    # A more primitive, but perfectly robust way of getting a reasonable
    # estimate for the maximum power.
    test_voltage = 0
    voltage_step = 0.001
    while power(test_voltage + voltage_step) > power(test_voltage):
        test_voltage += voltage_step

    max_power = power(test_voltage)

    # Calculate the maximized efficience
    efficiency = max_power / power_in

    # This segment isn't needed for functionality at all, but can display a
    # plot showing how the maximization of power by choosing the optimal
    # voltage value works
    if plot_current_voltage:
        V = np.linspace(0, 2, 200)
        plt.plot(V, J(V))
        plt.plot(V, power(V), linestyle='--')
        plt.savefig('pp.png')
        plt.close()
        # print(power(V_Pmax))

    return efficiency


if __name__ == '__main__':
    path = str(os.path.join(os.path.dirname(__file__),
                            '../vasp/examples/SiOptb88/MAIN-MBJ-bulk@mp_149/vasprun.xml'))
    en, abz, dirgap, indirgap = optics(path)
    abz = abz * 100.0
    eff = slme(en, abz, indirgap, indirgap, plot_current_voltage=False)
    print('SLME', 100 * eff)
    eff = calculate_SQ(indirgap)
    print('SQ', 100 * eff)
