from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.core.atoms import Atoms
import numpy as np


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
    nb = NeighborsAnalysis(Si).get_all_distributions
    tmp = round((nb["rdf"][-3]), 2)
    assert (tmp) == (4.08)


# test_nbors()
