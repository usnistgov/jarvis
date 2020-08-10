"""Module to generate networkx graphs."""
from jarvis.core.atoms import get_supercell_dims
import networkx as nx
from jarvis.core.specie import Specie
from jarvis.core.utils import random_colors


def get_networkx_graph(atoms=None, rcut=5.0, enforce_c_size=5.0):
    """Get Networkx graph. Requires Networkx installation."""
    color_dict = random_colors()
    graph = nx.Graph()
    dim = get_supercell_dims(atoms=atoms, enforce_c_size=enforce_c_size)
    atoms = atoms.make_supercell(dim)
    # nbr = NeighborsAnalysis(atom)
    # rcut1=nbr.rcut1
    # rcut2=nbr.rcut2
    # print ('rcut1,rcut2',rcut1,rcut2)
    adj = atoms.raw_distance_matrix
    adj[adj > rcut] = 0
    color_map = []
    mapping = {}
    for ii, i in enumerate(atoms.elements):
        color_map.append(color_dict[Specie(i).Z])
        mapping[str(ii)] = i
    for ii, i in enumerate(atoms.elements):
        graph.add_node(ii, weight=Specie(i).Z)

    for ii, i in enumerate(atoms.elements):
        for jj, j in enumerate(atoms.elements):
            graph.add_edge(ii, jj, weight=adj[ii, jj])
    return graph, color_map


"""
from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
atom = Atoms.from_dict(get_jid_data("JVASP-664")["atoms"])
graph, color_map = get_networkx_graph(atom)
nx.draw(graph, node_color=color_map, with_labels=True)
"""
