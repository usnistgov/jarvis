"""
Modules for analyzing infrared intensities.

Please find more details in https://doi.org/10.1038/s41524-020-0337-2 .
"""
from jarvis.core.spectrum import Spectrum
import numpy as np
import os


def normalize_vecs(phonon_eigenvectors, masses):
    """
    Return the eigenvectors after division of each component by sqrt(mass).

    Adapted from https://github.com/JMSkelton/Phonopy-Spectroscopy/
    TODO: include LO-TO splitting.
    """
    nmodes = len(phonon_eigenvectors)
    nmasses = len(masses)
    natoms = nmasses
    sqrt_masses = np.sqrt(masses)
    eigendisplacements = np.zeros((nmodes, natoms, 3), dtype="float64")
    for i in range(0, nmodes):
        eigenvector = phonon_eigenvectors[i]
        for j in range(0, natoms):
            eigendisplacements[i, j, :] = np.divide(
                eigenvector[j], sqrt_masses[j]
            )
    return eigendisplacements


def ir_intensity(
    phonon_eigenvectors=[],
    phonon_eigenvalues=[],
    masses=[],
    born_charges=[],
    factor=33.35641,
    nac=True,
    epsilon=[],
    enforce_positive_freqs=True,
    smoothen=True,
):
    """Calculate IR intensity using DFPT."""
    # TODO:add non-analytical correction
    eigendisplacements = normalize_vecs(phonon_eigenvectors, masses)
    becDim1, becDim2, becDim3 = np.shape(born_charges)
    freq = []
    ir_ints = []
    for i, v in zip(eigendisplacements, phonon_eigenvalues):
        eigendisplacement = i
        eigDim1, eigDim2 = np.shape(eigendisplacement)
        irIntensity = 0.0
        for a in range(0, 3):
            sumTemp1 = 0.0
            for j in range(0, eigDim1):
                sumTemp2 = 0.0
                for b in range(0, 3):
                    sumTemp2 += born_charges[j][a][b] * eigendisplacement[j][b]
                sumTemp1 += sumTemp2
            irIntensity += sumTemp1 ** 2
        if enforce_positive_freqs:
            if v > 0:
                freq.append(v * factor)  # Thz to cm-1
                ir_ints.append(irIntensity)
        else:
            freq.append(v * factor)  # Thz to cm-1
            ir_ints.append(irIntensity)
    if smoothen:
        freq, ir_ints = Spectrum(x=freq, y=ir_ints).smoothen_spiky_spectrum()
    return freq, ir_ints


def ir_intensity_phonopy(
    run_dir=".",
    vasprun="vasprun.xml",
    BornFileName="BORN",
    PoscarName="POSCAR",
    ForceConstantsName="FORCE_CONSTANTS",
    supercell=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    nac=True,
    symprec=1e-5,
    primitive=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    degeneracy_tolerance=1e-5,
    vector=[0, 0, 0],
    smoothen=False,
):
    """Calculate IR intensity using DFPT and phonopy."""
    from phonopy import Phonopy
    from phonopy.interface.vasp import read_vasp
    from phonopy.file_IO import (
        parse_BORN,
        parse_FORCE_CONSTANTS,
    )
    import shutil
    from phonopy.units import VaspToCm

    # from phonopy.phonon.degeneracy import (
    #    degenerate_sets as get_degenerate_sets,
    # )
    # adapted from https://github.com/JaGeo/IR
    # TODO: Make directory indepndent
    cwd = os.getcwd()
    print("Directory:", cwd)
    os.chdir(run_dir)
    if not os.path.exists(vasprun):
        shutil.copy2(vasprun, "vasprun.xml")
    cmd = str("phonopy --fc vasprun.xml")
    os.system(cmd)
    born_file = os.path.join(os.getcwd(), BornFileName)
    cmd = str("phonopy-vasp-born >  ") + str(born_file)
    os.system(cmd)
    from jarvis.io.vasp.outputs import Vasprun

    v = Vasprun(vasprun)
    strt = v.all_structures[0]
    strt.write_poscar(PoscarName)
    unitcell = read_vasp(PoscarName)
    phonon = Phonopy(
        unitcell,
        supercell_matrix=supercell,
        primitive_matrix=primitive,
        factor=VaspToCm,
        symprec=symprec,
    )
    natoms = phonon.get_primitive().get_number_of_atoms()
    force_constants = parse_FORCE_CONSTANTS(filename=ForceConstantsName)
    phonon.set_force_constants(force_constants)
    masses = phonon.get_primitive().get_masses()
    phonon.set_masses(masses)
    BORN_file = parse_BORN(phonon.get_primitive(), filename=BornFileName)
    BORN_CHARGES = BORN_file["born"]
    # print ('born_charges2',BORN_CHARGES)
    if nac:
        phonon.set_nac_params(BORN_file)
    frequencies, eigvecs = phonon.get_frequencies_with_eigenvectors(vector)
    # frequencies=VaspToTHz*frequencies/VaspToCm
    # x, y = ir_intensity(
    #    phonon_eigenvectors=np.real(eigvecs),
    #    phonon_eigenvalues=frequencies,
    #    masses=masses, #np.ones(len(masses)),
    #    born_charges=born_charges,
    #    smoothen=smoothen,
    # )
    NumberOfBands = len(frequencies)
    EigFormat = {}
    for alpha in range(NumberOfBands):
        laufer = 0
        for beta in range(natoms):
            for xyz in range(0, 3):
                EigFormat[beta, alpha, xyz] = eigvecs[laufer][alpha]
                laufer = laufer + 1
    Intensity = {}
    intensities = []
    for freq in range(len(frequencies)):
        Intensity[freq] = 0
        tmp = 0
        for alpha in range(3):
            asum = 0
            for n in range(natoms):
                for beta in range(3):
                    asum = asum + BORN_CHARGES[n, alpha, beta] * np.real(
                        EigFormat[n, freq, beta]
                    ) / np.sqrt(masses[n])
                    tmp += asum
            Intensity[freq] = Intensity[freq] + np.power(np.absolute(asum), 2)
        intensities.append(Intensity[freq])
    os.chdir(cwd)
    return frequencies, intensities


"""
if __name__ == "__main__":
    from jarvis.io.vasp.outputs import Vasprun, Outcar

    out = Outcar("../../tests/testfiles/io/vasp/OUTCAR.JVASP-39")
    vrun = Vasprun("../../tests/testfiles/io/vasp/vasprun.xml.JVASP-39")
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
        smoothen=False,
    )
    for i, j in zip(x, y):
        if j > 0.1:
            print("first", i, j)
    print()
    x, y = ir_intensity_phonopy()
    for i, j in zip(x, y):
        if j > 0.1:
            print("later", i, j)
    # for i, j in zip(phonon_eigenvalues, vrun_eigs):
    #    print(i, j)
"""
