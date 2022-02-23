#!/usr/bin/env python

from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
import numpy as np
import matplotlib.pyplot as plt


def eigenvector_overlap_at_q(qpt, branch, atm, mesh_dict):
    """
    Return array of eigenvector overlap values at specific q point
    """
    if atm == 0:
        vec0 = mesh_dict["eigenvectors"][qpt][branch][:3].conj()
    elif atm == 1:
        vec0 = mesh_dict["eigenvectors"][qpt][branch][3:6].conj()
    for qpt2 in range(len(mesh_dict["eigenvectors"])):
        for branch in range(len(mesh_dict["eigenvectors"][qpt2])):
            if atm == 0:
                vec = mesh_dict["eigenvectors"][qpt2][branch][:3]
            if atm == 1:
                vec = mesh_dict["eigenvectors"][qpt2][branch][3:6]
            overlap = np.sum(np.abs((vec * vec0).reshape(-1, 3).sum(axis=1)) ** 2)
    return overlap


if __name__ == "__main__":
    unitcell = read_vasp("POSCAR")
    phonon = Phonopy(
        unitcell,
        [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
        primitive_matrix=[[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )

    symmetry = phonon.get_symmetry()
    print("Space group: %s" % symmetry.get_international_table())

    force_sets = parse_FORCE_SETS()
    phonon.set_displacement_dataset(force_sets)
    phonon.produce_force_constants()
    primitive = phonon.get_primitive()

    # Born effective charges and dielectric constants are read from BORN file.
    nac_params = parse_BORN(primitive, filename="BORN")

    m = 6

    phonon.run_mesh([m, m, m], with_eigenvectors=True)

    mesh_dict = phonon.get_mesh_dict()

    overlap = eigenvector_overlap_at_q(qpt=0, branch=0, atm=0, mesh_dict=mesh_dict)

    """
    Overlap for Atom 1
    """
    overlap_atm0 = np.zeros([6, len(mesh_dict["qpoints"])])
    for b in range(len(mesh_dict["eigenvectors"][0])):
        for q in range(len(mesh_dict["qpoints"])):
            overlap_atm0[b, q] = eigenvector_overlap_at_q(
                q, branch=b, atm=0, mesh_dict=mesh_dict
            )

    plt.figure()
    plt.matshow(overlap_atm0, cmap=plt.cm.inferno, vmin=0, vmax=0.25)
    plt.colorbar()
    plt.xlabel("Q-point Index")
    plt.ylabel("Branch Index")  # Change into fraction of q_max
    plt.savefig("NaCl_atm0_overlap_matrix.pdf", bbox_inches="tight")

    """
    Overlap for Atom 2
    """
    overlap_atm1 = np.zeros([6, len(mesh_dict["qpoints"])])
    for b in range(len(mesh_dict["eigenvectors"][0])):
        for q in range(len(mesh_dict["qpoints"])):
            overlap_atm0[b, q] = eigenvector_overlap_at_q(
                q, branch=b, atm=1, mesh_dict=mesh_dict
            )

    plt.figure()
    plt.matshow(overlap_atm0, cmap=plt.cm.inferno, vmin=0, vmax=0.25)
    plt.colorbar()
    plt.xlabel("Q-point Index")
    plt.ylabel("Branch Index")  # Change into fraction of q_max
    plt.savefig("NaCl_atm1_overlap_matrix.pdf", bbox_inches="tight")
    """
    Density of States
    """
    phonon.auto_total_dos(mesh=[m, m, m], plot=True)
