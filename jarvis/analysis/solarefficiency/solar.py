"""
Modules for calculating theoretical solar-cell efficiency
See: https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.9b02166
"""
import numpy as np
import os
from jarvis.io.vasp.outputs import Vasprun
from scipy.interpolate import interp1d
from numpy import interp
import scipy.constants as constants
from scipy.integrate import simps
from scipy.constants import physical_constants, speed_of_light
import matplotlib.pyplot as plt




class SolarEfficiency(object):
    def __init__(self, formalism="slme"):
        """
        Calculates theoretical solar cell efficiency using SLME or SQ approach
        See https://pubs.acs.org/doi/10.1021/acs.chemmater.9b02166
        """
        self.formalism = formalism

    def calculate_SQ(
        self, bandgap_ev, temperature=300, fr=1, plot_current_voltage=False
    ):
        """
        
        Args:
        
              bandgap_ev: bandga in electron-volt
              
              temperature: temperature in K
     
                 
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
            / (
                (np.exp(h * c / (solar_spectra_wavelength_meters * k * temperature)))
                - 1.0
            )
        )

        # I've removed a pi in the equation above - Marnik Bercx

        # now to convert the irradiance from Power/m**2(m) into photon#/s*m**2(m)
        blackbody_photon_flux = blackbody_irradiance * (
            solar_spectra_wavelength_meters / (h * c)
        )

        # absorbance interpolation onto each solar spectrum wavelength

        # Get the bandgap in wavelength in meters
        bandgap_wavelength = h_e * c / bandgap_ev

        # Only take the part of the wavelength-dependent solar spectrum and
        # blackbody spectrum below the bandgap wavelength
        bandgap_index = np.searchsorted(
            solar_spectra_wavelength_meters, bandgap_wavelength
        )

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
            (blackbody_photon_flux[:bandgap_index], np.array([bandgap_blackbody])),
            axis=0,
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
        self,
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
            / (
                (np.exp(h * c / (solar_spectra_wavelength_meters * k * temperature)))
                - 1.0
            )
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

        # creates cubic spline interpolating function, set up to use end values
        #  as the guesses if leaving the region where data exists
        material_absorbance_data_function = interp1d(
            material_wavelength_for_absorbance_data,
            material_absorbance_data,
            kind="cubic",
            fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
            bounds_error=False,
        )

        material_interpolated_absorbance = np.zeros(
            len(solar_spectra_wavelength_meters)
        )
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

"""
if __name__ == "__main__":

    v = Vasprun(
        "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO/MAIN-MBJ-bulk@mp-149/vasprun.xml"
    )
    dirgap = v.get_dir_gap
    indirgap = v.get_indir_gap
    en, abz = v.avg_absorption_coefficient
    abz = abz * 100
    eff = SolarEfficiency().slme(
        en, abz, indirgap, indirgap, plot_current_voltage=False
    )
    print("SLME", 100 * eff)
    eff = SolarEfficiency().calculate_SQ(indirgap)
    print("SQ", 100 * eff)
"""
