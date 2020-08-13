from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.core.atoms import Atoms
import numpy as np
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')

def test_nbors():
    box = [[5.493642, 0, 0], [0, 5.493642, 0], [0, 0, 5.493642]]
    elements = ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"]
    coords = [
        [0, 0, 0],
        [0.25, 0.25, 0.25],
        [0.000000, 0.500000, 0.500000],
        [0.250000, 0.750000, 0.750000],
        [0.500000, 0.000000, 0.500000],
        [0.750000, 0.250000, 0.750000],
        [0.500000, 0.500000, 0.000000],
        [0.750000, 0.750000, 0.250000],
    ]
    coords = np.array(coords)
    # box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    # coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    # elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    nbr = NeighborsAnalysis(Si)
    nb = nbr.get_all_distributions
    tmp = round((nb["rdf"][-3]), 2)
    assert (tmp) == (4.08)
    nbr.get_rdf(plot=True)
    # nbr.ang_dist(nbor_info=info,plot=True)
    nbr.ang_dist_first(plot=True)
    nbr.ang_dist_second(plot=True)
    nbr.get_ddf(plot=True)
    angs = nbr.atomwise_angle_dist()
    ardf = nbr.atomwise_radial_dist()
    cmd = "rm *.png"
    os.system(cmd)
