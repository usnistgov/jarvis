"""Module for processing force-constant."""
import numpy as np


def qpoint(force_constant=[], qpt=[0.1, 0.1, 0.1]):
    """Get FC as natons x natons x 3 x3."""
    qpt = np.array(qpt)
    exp_iqpt = np.exp(1.0j * qpt)
    dmat = force_constant * exp_iqpt
    vals, vects = np.linalg.eigh(dmat)
    return vals, vects
