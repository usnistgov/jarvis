"""Module for processing force-constant."""
import numpy as np

"""
def dynamical_matrix(
    force_constant=[],
    atoms=None,
    nac_params={},
    num_G_points=300,
    g_rad=100,
    exp_cutoff=1e-10,
):
    G_cutoff = (3 * num_G_points / (4 * np.pi) / atoms.volume) ** (1.0 / 3)
    rec_lat = atoms.lattice.inv_lattice()
    pts = np.arange(-g_rad, g_rad + 1)
    grid = np.meshgrid(pts, pts, pts)
    for i in range(3):
        grid[i] = grid[i].ravel()
    G_vec_list = np.dot(rec_lat, grid).T
    G_norm2 = ((G_vec_list) ** 2).sum(axis=1)
    G_list = np.array(
        G_vec_list[G_norm2 < G_cutoff ** 2], dtype="double", order="C"
    )
    GeG = G_cutoff ** 2 * np.trace(nac_params["epsilon"]) / 3
    Lambda = np.sqrt(-GeG / 4 / np.log(exp_cutoff))
"""


def read_fc(filename="FORCE_CONSTANTS"):
    """Read force-constants."""
    from jarvis.io.phonopy.outputs import read_fc

    fc = read_fc(filename)
    return fc


def qpoint(force_constant=[], qpt=[0.0, 0.0, 0.0]):
    """Get FC as natons x natons x 3 x3."""
    qpt = np.array(qpt)
    exp_iqpt = np.exp(1.0j * qpt)
    dmat = force_constant * exp_iqpt
    vals, vects = np.linalg.eigh(dmat)
    return vals, vects
