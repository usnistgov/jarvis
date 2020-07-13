"""Module to calculate the dark-matter screening."""
import numpy as np


def peskin_schroeder_constant(bandgap=None, fermi_velocities=[], epsilon=[]):
    """
    Provide Peskin Scroeder constant for a material.

    Just proportionality constant. Doesn't take into account
    g and Lambda.

    Args:
       bandgap: benadgap in eV.

       fermi_velocities: vx, vy, vz

       epslion:  eigenvalues of epsilon tensor


    Return:
       1/(epsion*vF)
    """
    return 1 / (np.array(fermi_velocities) * np.array(epsilon))


"""
if __name__ == "__main__":
    epsilon = [187.5, 9.8, 90.9]
    fermi_velocities = [0.0029, 0.0050, 0.0021]
    print(
        peskin_schroeder_constant(
            bandgap=0.035, fermi_velocities=fermi_velocities, epsilon=epsilon
        )
    )
"""
