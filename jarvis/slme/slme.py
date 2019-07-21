from __future__ import unicode_literals, print_function

# Author: Marnik Bercx (University of Antwerp), Kamal Choudhary (NIST)
# Forked and adjusted from https://github.com/ldwillia/SL3ME

"""
Calculate spectroscopy limited maximum efficiency (SLME) given dielectric function data
"""

import os
import numpy as np
import matplotlib as plt
import scipy.constants as constants
from scipy.integrate import simps
import matplotlib.pyplot as plt

plt.switch_backend("agg")
import glob, os, math
from pymatgen.io.vasp.outputs import Vasprun
from numpy import loadtxt, arange, logspace
from math import pi, sqrt
from scipy.constants import physical_constants, speed_of_light
from pymatgen.electronic_structure.core import Spin


eV_to_recip_cm = 1.0 / (
    physical_constants["Planck constant in eV s"][0] * speed_of_light * 1e2
)


def nelec_out(out=""):
    f = open(out, "r")
    lines = f.read().splitlines()
    f.close()
    for i in lines:
        if "NELECT =" in i:
            nelec = int(
                float(i.split()[2])
            )  # int(i.split('NELECT =')[1].split('total number of electrons')[0])
    return nelec


def get_dir_indir_gap(run=""):
    """
 Get direct and indirect bandgaps for a vasprun.xml
 """

    v = Vasprun(run)
    bandstructure = v.get_band_structure()
    dir_gap = bandstructure.get_direct_band_gap()
    indir_gap = bandstructure.get_band_gap()["energy"]
    return dir_gap, indir_gap


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
    return (
        energies_in_eV,
        (
            2.0
            * np.sqrt(2.0)
            * pi
            * eV_to_recip_cm
            * energies_in_eV
            * np.sqrt(-epsilon_1 + np.sqrt(epsilon_1 ** 2 + epsilon_2 ** 2))
        ),
    )


def optics(ru=""):
    dirgap, indirgap = get_dir_indir_gap(ru)

    run = Vasprun(ru, occu_tol=1e-2)
    new_en, new_abs = absorption_coefficient(run.dielectric)
    return (
        np.array(new_en, dtype=np.float64),
        np.array(new_abs, dtype=np.float64),
        dirgap,
        indirgap,
    )


def calculate_SQ(bandgap_ev, temperature=300, fr=1, plot_current_voltage=False):
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
    solar_spectrum_data_file = str(
        os.path.join(os.path.dirname(__file__), "am1.5G.dat")
    )
    solar_spectra_wavelength, solar_spectra_irradiance = np.loadtxt(
        solar_spectrum_data_file, usecols=[0, 1], unpack=True, skiprows=2
    )

    solar_spectra_wavelength_meters = solar_spectra_wavelength * 1e-9

    # need to convert solar irradiance from Power/m**2(nm) into
    # photon#/s*m**2(nm) power is Watt, which is Joule / s
    # E = hc/wavelength
    # at each wavelength, Power * (wavelength(m)/(h(Js)*c(m/s))) = ph#/s
    solar_spectra_photon_flux = solar_spectra_irradiance * (
        solar_spectra_wavelength_meters / (h * c)
    )

    ### Calculation of total solar power incoming
    power_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    # calculation of blackbody irradiance spectra
    # units of W/(m**3), different than solar_spectra_irradiance!!! (This
    # is intentional, it is for convenience)
    blackbody_irradiance = (
        2.0 * h * c ** 2 / (solar_spectra_wavelength_meters ** 5)
    ) * (
        1.0
        / ((np.exp(h * c / (solar_spectra_wavelength_meters * k * temperature))) - 1.0)
    )

    # I've removed a pi in the equation above - Marnik Bercx

    # now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (
        solar_spectra_wavelength_meters / (h * c)
    )

    # absorbance interpolation onto each solar spectrum wavelength
    from numpy import interp

    # Get the bandgap in wavelength in meters
    bandgap_wavelength = h_e * c / bandgap_ev

    # Only take the part of the wavelength-dependent solar spectrum and
    # blackbody spectrum below the bandgap wavelength
    bandgap_index = np.searchsorted(solar_spectra_wavelength_meters, bandgap_wavelength)

    bandgap_irradiance = interp(
        np.array([bandgap_wavelength]),
        solar_spectra_wavelength_meters,
        solar_spectra_photon_flux,
    )

    bandgap_blackbody = (
        (2.0 * h * c ** 2 / (bandgap_wavelength ** 5))
        * (1.0 / ((np.exp(h * c / (bandgap_wavelength * k * temperature))) - 1.0))
        * (bandgap_wavelength / (h * c))
    )

    integration_wavelength = np.concatenate(
        (
            solar_spectra_wavelength_meters[:bandgap_index],
            np.array([bandgap_wavelength]),
        ),
        axis=0,
    )

    integration_solar_flux = np.concatenate(
        (solar_spectra_photon_flux[:bandgap_index], bandgap_irradiance), axis=0
    )

    integration_blackbody = np.concatenate(
        (blackbody_photon_flux[:bandgap_index], np.array([bandgap_blackbody])), axis=0
    )

    #  Numerically integrating irradiance over wavelength array
    # Note: elementary charge, not math e!  ## units of A/m**2   W/(V*m**2)
    J_0_r = e * np.pi * simps(integration_blackbody, integration_wavelength)

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
        plt.plot(V, power(V), linestyle="--")
        plt.show()
        print(max_power)

    return efficiency


def slme(
    material_energy_for_absorbance_data,
    material_absorbance_data,
    material_direct_allowed_gap,
    material_indirect_gap,
    thickness=50e-6,
    temperature=293.15,
    absorbance_in_inverse_centimeters=False,
    cut_off_absorbance_below_direct_allowed_gap=True,
    plot_current_voltage=False,
):
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
    solar_spectrum_data_file = str(
        os.path.join(os.path.dirname(__file__), "am1.5G.dat")
    )

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
        solar_spectra_wavelength_meters / (h * c)
    )

    # Calculation of total solar power incoming
    power_in = simps(solar_spectra_irradiance, solar_spectra_wavelength)

    # calculation of blackbody irradiance spectra
    # units of W/(m**3), different than solar_spectra_irradiance!!! (This
    # is intentional, it is for convenience)
    blackbody_irradiance = (
        2.0 * h * c ** 2 / (solar_spectra_wavelength_meters ** 5)
    ) * (
        1.0
        / ((np.exp(h * c / (solar_spectra_wavelength_meters * k * temperature))) - 1.0)
    )

    # I've removed a pi in the equation above - Marnik Bercx

    # now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
    blackbody_photon_flux = blackbody_irradiance * (
        solar_spectra_wavelength_meters / (h * c)
    )

    # units of nm
    material_wavelength_for_absorbance_data = (
        (c * h_e) / (material_energy_for_absorbance_data + 0.00000001)
    ) * 10 ** 9

    # absorbance interpolation onto each solar spectrum wavelength
    from scipy.interpolate import interp1d

    # creates cubic spline interpolating function, set up to use end values
    #  as the guesses if leaving the region where data exists
    material_absorbance_data_function = interp1d(
        material_wavelength_for_absorbance_data,
        material_absorbance_data,
        kind="cubic",
        fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
        bounds_error=False,
    )

    material_interpolated_absorbance = np.zeros(len(solar_spectra_wavelength_meters))
    for i in range(0, len(solar_spectra_wavelength_meters)):
        ## Cutting off absorption data below the gap. This is done to deal
        # with VASPs broadening of the calculated absorption data

        if (
            solar_spectra_wavelength[i]
            < 1e9 * ((c * h_e) / material_direct_allowed_gap)
            or cut_off_absorbance_below_direct_allowed_gap == False
        ):
            material_interpolated_absorbance[i] = material_absorbance_data_function(
                solar_spectra_wavelength[i]
            )

    absorbed_by_wavelength = 1.0 - np.exp(
        -2.0 * material_interpolated_absorbance * thickness
    )

    #  Numerically integrating irradiance over wavelength array
    # Note: elementary charge, not math e!  ## units of A/m**2   W/(V*m**2)
    J_0_r = (
        e
        * np.pi
        * simps(
            blackbody_photon_flux * absorbed_by_wavelength,
            solar_spectra_wavelength_meters,
        )
    )

    J_0 = J_0_r / fr

    #  Numerically integrating irradiance over wavelength array
    # elementary charge, not math e!  ### units of A/m**2   W/(V*m**2)
    J_sc = e * simps(
        solar_spectra_photon_flux * absorbed_by_wavelength, solar_spectra_wavelength
    )

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
        plt.plot(V, power(V), linestyle="--")
        plt.savefig("pp.png")
        plt.close()
        # print(power(V_Pmax))

    return efficiency


if __name__ == "__main__":
    path = str(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "vasp",
            "examples",
            "SiOptb88",
            "MAIN-MBJ-bulk@mp_149",
            "vasprun.xml",
        )
    )
    en, abz, dirgap, indirgap = optics(path)
    abz = abz * 100.0
    eff = slme(en, abz, indirgap, indirgap, plot_current_voltage=False)
    print("SLME", 100 * eff)
    eff = calculate_SQ(indirgap)
    print("SQ", 100 * eff)
