from jarvis.io.vasp.outputs import Vasprun, Outcar
import numpy as np


def normalize_vecs(phonon_eigenvectors, masses):
    """
    Return the eigenvectors after division of each component by sqrt(mass).
    """

    nmodes = len(phonon_eigenvectors)
    nmasses = len(masses)
    natoms = nmasses

    sqrt_masses = np.sqrt(masses)

    eigendisplacements = np.zeros((nmodes, natoms, 3), dtype=np.float64)

    for i in range(0, nmodes):
        eigenvector = phonon_eigenvectors[i]
        for j in range(0, natoms):
            eigendisplacements[i, j, :] = np.divide(eigenvector[j], sqrt_masses[j])
    return eigendisplacements


def ir_intensity(
    phonon_eigenvectors=[], phonon_eigenvalues=[], masses=[], born_charges=[]
):
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
        freq.append(v * 33.35641)  # Thz to cm-1
        ir_ints.append(irIntensity)
    return freq, ir_ints

"""
if __name__ == "__main__":

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
    )
    for i, j in zip(x, y):
        if j > 0.1:
            print(i, j)
    print ()
    print (round(y[1],2))
    #for i, j in zip(phonon_eigenvalues, vrun_eigs):
    #    print(i, j)
"""
