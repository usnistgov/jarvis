#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 21:34:38 2021

@author: rlg3

Module for post-processing phono3py output
kappa-mxxx.hdf, jdos-mxxx-gy-tz.hdf, phono3py_disp.yaml


Notes for inputs.py:
    tags : generate automatically from the inputs.py. Include info about
    isotope scattering, boundary scattering, etc.

Think I should initialize a Phonopy object to run these.
    
Might need a mesh_dict from phonon.get_mesh_dict() so that I have the phonon
bandstructure properties

Add function for mode gruneisen from third order force constants
"""

import h5py
from jarvis.core.specie import Specie
from jarvis.io.phonopy.outputs import (
    total_dos,
    get_Phonopy_obj,
    get_spectral_heat_capacity,
)
import numpy as np
import spglib
import matplotlib.pyplot as plt

try:
    from phonopy import Phonopy
    from phonopy import load
except Exception as exp:
    print("Phonopy is not installed.", exp)
    pass


"""
Constants 
"""
kB = 1.38064852e-23
hbar = 1.0545718e-34
Na = 6.0221409e23
kappa_unit_conversion = 6.358245562444196


def gruneisen_approximation(vt, vl):
    x = vt / vl
    return (3 / 2) * (3 - 4 * x ** 2) / (1 + 2 * x ** 2)


class Kappa:
    def __init__(
        self,
        filename="",
        mesh=None,
        total_dos_dat="",
        temperatures=None,
        kappa_format="scalar_xx",  # need this?
        composition=None,
        tags=None,
    ):
        """

        Parameters
        ----------
        filename : string,required
            Name of kappa-mxxxx.hdf file. The default is "".
        temperatures : list, optional
            List of temperatures to use for post-processing. Default is None.
            If None, uses all temperatures computed in thermal
            conductivity run.
        kappa_format : string, optional
            Desired kappa output format. The default is 'scalar_xx'.
            Other choices : "tensor" or ...
        composition : Composition object of jarvis.core.composition, optional
            Composition object for material used in property calculation.
            The default is None.
        tags : string, optional
        Descriptors of thermal conductivity calculation. The default is None.
        """
        f = h5py.File(filename, "r")
        self.dict = dict(f)
        if temperatures:
            self.temperatures = temperatures
        else:
            self.temperatures = list(self.dict["temperature"])
        self.composition = composition
        self.tags = tags
        self.n_qpoint = np.shape(f["qpoint"])[0]
        self.total_dos = np.array(total_dos(total_dos_dat))
        self.kappa_format = kappa_format

    def to_dict(self):
        return self.dict

    def kappa(self, T):
        if T not in self.temperatures:
            raise Exception("Kappa not evaluated at this temperature.")
        T_indx = self.temperatures.index(T)
        if self.kappa_format == "tensor":
            return self.dict["kappa"][T_indx]
        if self.kappa_format == "scalar_xx":
            return self.dict["kappa"][T_indx][0]


class JDOS:
    def __init__(self, phonopy_obj, directory, mesh=[1, 1, 1], temperature=None):
        """
        Parameters
        ----------
        directory : string, required
            DESCRIPTION. The default is "".
        qpoints : TYPE, optional
            DESCRIPTION. The default is None.
        temperature : integer, optional
        Temperature for a weighted JDOS calculation. The default is None. When
        None, unweighted JDOS is computed.

        """
        self.mesh = mesh
        self.temperature = temperature
        self.directory = directory
        self.phonopy_obj = phonopy_obj
        phonopy_obj.run_mesh(mesh, with_group_velocities=True)
        self.mesh_dict = phonopy_obj.get_mesh_dict()

    def select_jdos(self):
        """

        Returns
        -------
        None.

        """

        def get_gridpts(self):
            """
            Generates list of gridpoint indices for JDOS calculation

            Parameters
            ----------
            phonopy_obj : TYPE
                DESCRIPTION.

            Returns
            -------
            None.

            """
            latt_vecs = self.phonopy_obj.get_primitive().get_cell()
            positions = self.phonopy_obj.get_primitive().get_scaled_positions()
            atom_type = self.phonopy_obj.get_primitive().get_atomic_numbers()
            cell = (latt_vecs, positions, atom_type)
            mapping, grid = spglib.get_ir_reciprocal_mesh(
                self.mesh, cell, is_shift=[0, 0, 0]
            )
            return mapping

        gridpt_list = get_gridpts(self)
        gridpt_uids = np.unique(gridpt_list)

        """
        Dictionary stores JDOS spectrum for each irreducible q-point
        key : irreducible grid point index
        value : JDOS spectrum [frequency, N1 process, N2 process]
        """
        jdos_ir = np.zeros([len(gridpt_uids), len(self.mesh_dict["frequencies"][0])])

        for g, gp in enumerate(gridpt_uids):
            if self.temperature is not None:
                file = (
                    self.directory
                    + "jdos-m"
                    + "".join(str(m) for m in self.mesh)
                    + "-g"
                    + str(gp)
                    + "-t"
                    + str(self.temperature)
                    + ".dat"
                )
            else:
                file = (
                    self.directory
                    + "jdos-m"
                    + "".join(str(m) for m in self.mesh)
                    + "-g"
                    + str(gp)
                    + ".dat"
                )
            data = np.genfromtxt(file, dtype="float")
            freq = data[:, 0]
            N2 = data[:, 1] + data[:, 2]
            indices = np.digitize(
                self.mesh_dict["frequencies"][g], freq
            )  # select freqeuncy bin of branch frequency
            # JDOS values for each irreducible q-point
            jdos_ir[g] = np.array(
                [(N2[i] + N2[i - 1]) / 2 for i in indices]
            )  # average the ceiling and floor of the bin
            # Map JDOS back to the full q-point mesh
        #     jdos_full = []
        # for i, qpoint in enumerate(gridpt_list):
        #     jdos_full.append(jdos_ir[qpoint])

        return jdos_ir

    # For spectral quantities that already include factor of DOS
    def mode_to_spectral_1(self, mode_jdos):
        self.phonopy_obj.run_total_dos()
        # Get tetrahedron mesh object
        thm = self.phonopy_obj._total_dos._tetrahedron_mesh
        thm.set(
            value="I", frequency_points=self.phonopy_obj._total_dos._frequency_points
        )
        spectral_jdos = np.zeros_like(self.phonopy_obj._total_dos._frequency_points)
        for i, iw in enumerate(thm):
            spectral_jdos += np.sum(iw * mode_jdos[i], axis=1)
        return spectral_jdos

    # For spectral quantities that need to be scaled by DOS
    def mode_to_spectral_2(self, mode_jdos):
        self.phonopy_obj.run_total_dos()
        # Get tetrahedron mesh object
        thm = self.phonopy_obj._total_dos._tetrahedron_mesh
        thm.set(
            value="I", frequency_points=self.phonopy_obj._total_dos._frequency_points
        )
        spectral_jdos = np.zeros_like(self.phonopy_obj._total_dos._frequency_points)
        for i, iw in enumerate(thm):
            spectral_jdos += np.sum(
                iw * mode_jdos[i] * self.phonopy_obj._total_dos._weights[i], axis=1
            )
        return spectral_jdos

    def plot_jdos(self, spectral_jdos):
        freq_pts = self.phonopy_obj._total_dos._frequency_points
        plt.figure()
        plt.plot(freq_pts, spectral_jdos)
        plt.ylabel(r"JDOS (THz$^{-1}$)")
        plt.xlabel(r"Frequency (THz)")

    # Weighted JDOS should work? Maybe write separate method?
    def linewidth_from_jdos(
        self, spectral_jdos, atoms, vs, grun=0.8, T=300, plot=False
    ):
        """
        

        Parameters
        ----------
        spectral_jdos : TYPE
           Currently only takes unweighted jdos values.
        atoms : TYPE
            DESCRIPTION.
        vs : float
            Sound velocity. (Group velocity may be more accurate?)
        gamma : TYPE, optional
            DESCRIPTION. The default is 1.
        T : TYPE, optional
            DESCRIPTION. The default is 300.

        Returns
        -------
        None.
        
        Need to write hbar in terms of THz? Or convert other values to seconds

        """
        prefactor = np.pi * kB * T / 6 / 3  # Added the factor of 3!!
        freq_pts = self.phonopy_obj._total_dos._frequency_points
        species = [Specie(i) for i in atoms.elements]
        N = len(species)
        avgM = sum([species[j].atomic_mass / Na / 1e3 for j in range(N)]) / N
        spectral_2Gamma = (
            prefactor * (grun ** 2 / (avgM * vs ** 2)) * freq_pts ** 2 * spectral_jdos
        )
        if plot:
            plt.figure()
            plt.plot(freq_pts, spectral_2Gamma)
            plt.xlabel("Frequency (THz)")
            plt.ylabel(r"2$\Gamma$ (THz)")
        return spectral_2Gamma

    def kappa_from_linewidth(self, spectral_2Gamma, component=["x", "x"], T=300):
        """
        Currently only works for scalar_xx kappa format
        """
        ij_dict = {"x": 0, "y": 1, "z": 3}
        freq_pts = self.phonopy_obj._total_dos._frequency_points
        find_zeros = np.where(spectral_2Gamma == 0)[0]
        # freq_pts = np.delete(freq_pts, find_zeros)
        mode_vg2_ij = jdos.mesh_dict["group_velocities"][
            :, :, ij_dict[component[0]]
        ]  # *\
        # jdos.mesh_dict['group_velocities'][:, :, ij_dict[component[1]]]
        spectral_vg2 = self.mode_to_spectral_2(
            mode_vg2_ij
        )  # Need to check group velocity units !!
        plt.figure()
        plt.plot(freq_pts, spectral_vg2)
        plt.xlabel("Frequency (THz)")
        plt.ylabel("Squared group velocity")

        spectral_Cp = get_spectral_heat_capacity(self.phonopy_obj, self.mesh, T)
        plt.xlabel("Frequency (THz)")
        plt.ylabel(r"C (eV/K$\cdot$THz)")
        spectral_2Gamma = np.delete(spectral_2Gamma, find_zeros)
        spectral_vg2 = np.delete(spectral_vg2, find_zeros)
        spectral_Cp = np.delete(spectral_Cp, find_zeros)
        spectral_kappa = (
            kappa_unit_conversion * spectral_vg2 * (1 / spectral_2Gamma) * spectral_Cp
        )
        red_freq_pts = np.delete(freq_pts, find_zeros)
        plt.figure()
        plt.plot(red_freq_pts, spectral_kappa)
        plt.xlabel("Frequency (THz)")
        plt.ylabel(r"$\kappa$ (W/m$\cdot$K$\cdot$THz)")


if __name__ == "__main__":
    kappa_Si = Kappa(
        "Si-testing/kappa-m111111.hdf5",
        total_dos_dat="../phonopy/Si-testing/total_dos.dat",
    )
    RT_kappa = kappa_Si.kappa(300.0)
    from jarvis.core.atoms import Atoms

    test_dir = "Si-testing/"
    pos = test_dir + "POSCAR-unitcell"
    atoms = Atoms.from_poscar(pos)
    phonon_obj = get_Phonopy_obj(
        atoms,
        phonopy_yaml="../phonopy/Si-testing/phonopy.yaml",
        FC_file="../phonopy/Si-testing/FORCE_CONSTANTS",
        scell=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    )
    jdos_dir = "Si-testing/unweighted_jdos/"
    jdos = JDOS(phonon_obj, directory=jdos_dir, mesh=[11, 11, 11])
    jdos_ir = jdos.select_jdos()
    spectral_jdos = jdos.mode_to_spectral_2(jdos_ir)
    jdos.plot_jdos(spectral_jdos)

    spectral_2Gamma = jdos.linewidth_from_jdos(spectral_jdos, atoms, vs=6084, plot=True)
    jdos.kappa_from_linewidth(spectral_2Gamma)

    spectral_kappa = jdos.mode_to_spectral_1(kappa_Si.dict["mode_kappa"][30, :, :, 0])
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"$\kappa$ (W/m$\cdot$K$\cdot$THz)")

    freq_pts = jdos.phonopy_obj._total_dos._frequency_points
    plt.figure()
    plt.plot(freq_pts, spectral_kappa)
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"$\kappa$ (W/m$\cdot$K$\cdot$THz)")

    spectral_C = get_spectral_heat_capacity(jdos.phonopy_obj, jdos.mesh, T=30)
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"C (eV/K$\cdot$THz)")

    spectral_vg2_ph3 = jdos.mode_to_spectral_2(kappa_Si.dict["gv_by_gv"][:, :, 0])
    plt.figure()
    plt.plot(freq_pts, spectral_vg2_ph3)
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"Phono3py Group Velocity (THz)")

    spectral_gamma = jdos.mode_to_spectral_2(kappa_Si.dict["gamma"][30, :, :])
    plt.figure()
    plt.plot(freq_pts, 2 * spectral_gamma)
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"2$\Gamma$ (THz)")

    spectral_Cp = jdos.mode_to_spectral_1(kappa_Si.dict["heat_capacity"][3, :, :])
    plt.figure()
    plt.plot(freq_pts, spectral_Cp)
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"C (eV/K$\cdot$THz)")

    # kappa_unit_conversion = kappa_Si.dict['kappa_unit_conversion'][()]
    # kappa_spec = kappa_unit_conversion * spectral_vg2_ph3 *\
    #     (1 / (2 * spectral_gamma)) * spectral_Cp

    # plt.figure()
    # plt.plot(freq_pts, kappa_spec)
    # plt.xlabel('Frequency (THz)')
    # plt.ylabel(r'$\kappa$ (W/m$\cdot$K$\cdot$THz)')

    plt.figure()
    plt.plot(freq_pts, spectral_2Gamma - (2 * spectral_gamma))
    plt.xlabel("Frequency (THz)")
    plt.ylabel(r"2$\Gamma$ Difference")
