"""Module to generate networkx graphs."""
from jarvis.core.atoms import get_supercell_dims
from jarvis.core.specie import Specie
from jarvis.core.utils import random_colors
import numpy as np
from collections import OrderedDict
from jarvis.analysis.structure.neighbors import NeighborsAnalysis


class Graph(object):
    """Generate a graph object."""

    def __init__(
        self,
        nodes=[],
        node_attributes=[],
        edges=[],
        edge_attributes=[],
        color_map=None,
        labels=None,
    ):
        """
        Initialize the graph object.

        Args:
            nodes: IDs of the graph nodes as integer array.

            node_attributes: node features as multi-dimensional array.

            edges: connectivity as a (u,v) pair where u is
                   the source index and v the destination ID.

            edge_attributes: attributes for each connectivity.
                             as simple as euclidean distances.
        """
        self.nodes = nodes
        self.node_attributes = node_attributes
        self.edges = edges
        self.edge_attributes = edge_attributes
        self.color_map = color_map
        self.labels = labels

    @staticmethod
    def from_atoms(
        atoms=None,
        lengthscale=0.5,
        variance=1.0,
        get_prim=False,
        zero_diag=False,
        node_atomwise_angle_dist=False,
        node_atomwise_rdf=False,
        features="basic",
        enforce_c_size=10.0,
        max_n=100,
        max_cut=5.0,
        verbose=False,
        make_colormap=True,
    ):
        """
        Get Networkx graph. Requires Networkx installation.

        Args:
             atoms: jarvis.core.Atoms object.

             rcut: cut-off after which distance will be set to zero
                   in the adjacency matrix.

             features: Node features.
                       'atomic_number': graph with atomic numbers only.
                       'cfid': 438 chemical descriptors from CFID.
                       'basic':10 features
                       array: array with CFID chemical descriptor names.
                       See: jarvis/core/specie.py

             enforce_c_size: minimum size of the simulation cell in Angst.
        """
        if get_prim:
            atoms = atoms.get_primitive_atoms
        dim = get_supercell_dims(atoms=atoms, enforce_c_size=enforce_c_size)
        atoms = atoms.make_supercell(dim)
        raw_data = np.array(atoms.raw_distance_matrix)
        adj = variance * np.exp(-raw_data / lengthscale)
        if zero_diag:
            np.fill_diagonal(adj, 0.0)
        nodes = np.arange(atoms.num_atoms)
        if features == "atomic_number":
            node_attributes = np.array(
                [[np.array(Specie(i).Z)] for i in atoms.elements],
                dtype="float",
            )
        elif features == "basic":
            feats = [
                "Z",
                "coulmn",
                "row",
                "X",
                "atom_rad",
                "nsvalence",
                "npvalence",
                "ndvalence",
                "nfvalence",
                "first_ion_en",
                "elec_aff",
            ]
            node_attributes = []
            for i in atoms.elements:
                tmp = []
                for j in feats:
                    tmp.append(Specie(i).element_property(j))
                node_attributes.append(tmp)
            node_attributes = np.array(node_attributes, dtype="float")
        elif features == "cfid":
            node_attributes = np.array(
                [np.array(Specie(i).get_descrp_arr) for i in atoms.elements],
                dtype="float",
            )
        elif isinstance(features, list):
            node_attributes = []
            for i in atoms.elements:
                tmp = []
                for j in features:
                    tmp.append(Specie(i).element_property(j))
                node_attributes.append(tmp)
            node_attributes = np.array(node_attributes, dtype="float")
        else:
            raise ("Please check the input options.")
        if node_atomwise_rdf or node_atomwise_angle_dist:
            nbr = NeighborsAnalysis(
                atoms, max_n=max_n, verbose=verbose, max_cut=max_cut
            )
        if node_atomwise_rdf:
            node_attributes = np.concatenate(
                (node_attributes, nbr.atomwise_radial_dist()), axis=1
            )
            node_attributes = np.array(node_attributes, dtype="float")
        if node_atomwise_angle_dist:
            node_attributes = np.concatenate(
                (node_attributes, nbr.atomwise_angle_dist()), axis=1
            )
            node_attributes = np.array(node_attributes, dtype="float")
        uv = []
        edge_features = []
        for ii, i in enumerate(atoms.elements):
            for jj, j in enumerate(atoms.elements):
                uv.append((ii, jj))
                edge_features.append(adj[ii, jj])
        edge_attributes = edge_features
        if make_colormap:
            sps = atoms.uniq_species
            color_dict = random_colors(number_of_colors=len(sps))
            new_colors = {}
            for i, j in color_dict.items():
                new_colors[sps[i]] = j
            color_map = []
            for ii, i in enumerate(atoms.elements):
                color_map.append(new_colors[i])
        return Graph(
            nodes=nodes,
            edges=uv,
            node_attributes=np.array(node_attributes),
            edge_attributes=np.array(edge_attributes),
            color_map=color_map,
        )

    def to_networkx(self):
        """Get networkx representation."""
        import networkx as nx

        graph = nx.DiGraph()
        graph.add_nodes_from(self.nodes)
        graph.add_edges_from(self.edges)
        for i, j in zip(self.edges, self.edge_attributes):
            graph.add_edge(i[0], i[1], weight=j)
        return graph

    @property
    def num_nodes(self):
        """Return number of nodes in the graph."""
        return len(self.nodes)

    @property
    def num_edges(self):
        """Return number of edges in the graph."""
        return len(self.edges)

    @classmethod
    def from_dict(self, d={}):
        """Constuct class from a dictionary."""
        return Graph(
            nodes=d["nodes"],
            edges=d["edges"],
            node_attributes=d["node_attributes"],
            edge_attributes=d["edge_attributes"],
            color_map=d["color_map"],
            labels=d["labels"],
        )

    def to_dict(self):
        """Provide dictionary representation of the Graph object."""
        info = OrderedDict()
        info["nodes"] = np.array(self.nodes).tolist()
        info["edges"] = np.array(self.edges).tolist()
        info["node_attributes"] = np.array(self.node_attributes).tolist()
        info["edge_attributes"] = np.array(self.edge_attributes).tolist()
        info["color_map"] = np.array(self.color_map).tolist()
        info["labels"] = np.array(self.labels).tolist()
        return info

    def __repr__(self):
        """Provide representation during print statements."""
        return "Graph({})".format(self.to_dict())

    @property
    def adjacency_matrix(self):
        """Provide adjacency_matrix of graph."""
        return np.array(self.edge_attributes).reshape(
            self.num_nodes, self.num_nodes
        )


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    from jarvis.db.figshare import get_jid_data

    atoms = Atoms.from_dict(get_jid_data("JVASP-664")["atoms"])
    g = Graph.from_atoms(
        atoms=atoms,
        features="basic",
        get_prim=True,
        zero_diag=True,
        node_atomwise_angle_dist=True,
        node_atomwise_rdf=True,
    )
    g = Graph.from_atoms(
        atoms=atoms,
        features="cfid",
        get_prim=True,
        zero_diag=True,
        node_atomwise_angle_dist=True,
        node_atomwise_rdf=True,
    )
    g = Graph.from_atoms(
        atoms=atoms,
        features="atomic_number",
        get_prim=True,
        zero_diag=True,
        node_atomwise_angle_dist=True,
        node_atomwise_rdf=True,
    )
    g = Graph.from_atoms(atoms=atoms, features="basic")
    g = Graph.from_atoms(
        atoms=atoms, features=["Z", "atom_mass", "max_oxid_s"]
    )
    g = Graph.from_atoms(atoms=atoms, features="cfid")
    # print(g)
    d = g.to_dict()
    g = Graph.from_dict(d)
    num_nodes = g.num_nodes
    num_edges = g.num_edges
    print(num_nodes, num_edges)
    assert num_nodes == 48
    assert num_edges == 2304
    assert len(g.adjacency_matrix) == 2304
    # graph, color_map = get_networkx_graph(atoms)
    # nx.draw(graph, node_color=color_map, with_labels=True)
    # from jarvis.analysis.structure.neighbors import NeighborsAnalysis
"""
