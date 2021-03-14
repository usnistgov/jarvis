"""Module to generate networkx graphs."""
from jarvis.core.atoms import get_supercell_dims
from jarvis.core.specie import Specie
from jarvis.core.utils import random_colors
import numpy as np
from collections import OrderedDict
from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.core.specie import get_node_attributes
import warnings
from typing import List, Tuple
from jarvis.core.atoms import Atoms
from collections import defaultdict

try:
    import torch
    from tqdm import tqdm
    import dgl
except Exception as exp:
    warnings.warn("dgl/torch/tqdm is not installed.", exp)
    pass


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
    def atom_dgl_multigraph(
        atoms=None,
        cutoff=8.0,
        max_neighbors=12,
        atom_features="cgcnn",
        enforce_undirected=False,
        max_attempts=3,
        id=None,
    ):
        """Obtain a DGLGraph for Atoms object."""
        all_neighbors = atoms.get_all_neighbors(r=cutoff)
        # if a site has too few neighbors, increase the cutoff radius
        min_nbrs = min(len(neighborlist) for neighborlist in all_neighbors)
        # print('min_nbrs,max_neighbors=',min_nbrs,max_neighbors)
        attempt = 0
        while min_nbrs < max_neighbors:
            print("extending cutoff radius!", attempt, cutoff, id)
            lat = atoms.lattice
            r_cut = max(cutoff, lat.a, lat.b, lat.c)
            attempt += 1
            if attempt >= max_attempts:
                atoms = atoms.make_supercell([2, 2, 2])
                print(
                    "Making supercell, exceeded,attempts",
                    max_attempts,
                    "cutoff",
                    r_cut,
                    id,
                )
            cutoff = r_cut
            all_neighbors = atoms.get_all_neighbors(r=cutoff)
            min_nbrs = min(len(neighborlist) for neighborlist in all_neighbors)
            # return Graph.atom_dgl_multigraph(
            #    atoms, r_cut, max_neighbors, atom_features
            # )

        # build up edge list
        # Currently there's no guarantee that this creates undirected graphs
        # An undirected solution would build the full edge list where nodes are
        # keyed by (index,image), and ensure each edge has a complementary edge

        # indeed,JVASP-59628 is an example of a calculation where this produces
        # a graph where one site has no incident edges!

        # build an edge dictionary u -> v
        # so later we can run through the dictionary
        # and remove all pairs of edges
        # so what's left is the odd ones out
        edges = defaultdict(list)

        u, v, r = [], [], []
        for site_idx, neighborlist in enumerate(all_neighbors):

            # sort on distance
            neighborlist = sorted(neighborlist, key=lambda x: x[2])

            ids = np.array([nbr[1] for nbr in neighborlist])
            distances = np.array([nbr[2] for nbr in neighborlist])
            c = np.array([nbr[3] for nbr in neighborlist])

            # find the distance to the k-th nearest neighbor
            max_dist = distances[max_neighbors - 1]

            # keep all edges out to the neighbor shell of the k-th neighbor
            ids = ids[distances <= max_dist]
            c = c[distances <= max_dist]
            distances = distances[distances <= max_dist]

            u.append([site_idx] * len(ids))
            v.append(ids)
            r.append(distances)

            # keep track of cell-resolved edges
            # to enforce undirected graph construction
            for dst, cell_id in zip(ids, c):
                u_key = f"{site_idx}-(0.0, 0.0, 0.0)"
                v_key = f"{dst}-{tuple(cell_id)}"
                edge_key = tuple(sorted((u_key, v_key)))
                edges[edge_key].append((site_idx, dst))

        if enforce_undirected:
            # add complementary edges to unpaired edges
            for edge_pair in edges.values():
                if len(edge_pair) == 1:
                    src, dst = edge_pair[0]
                    u.append(dst)  # swap the order!
                    v.append(src)
                    r.append(atoms.raw_distance_matrix[src, dst])

        u = torch.tensor(np.hstack(u))
        v = torch.tensor(np.hstack(v))
        r = torch.tensor(np.hstack(r)).type(torch.get_default_dtype())

        # build up atom attribute tensor
        species = atoms.elements
        node_features = torch.tensor(
            [
                get_node_attributes(s, atom_features=atom_features)
                for s in species
            ]
        ).type(torch.get_default_dtype())

        g = dgl.graph((u, v))
        g.ndata["atom_features"] = node_features
        g.edata["bondlength"] = r

        return g

    @staticmethod
    def from_atoms(
        atoms=None,
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
                       'atomic_fraction': graph with atomic fractions
                                         in 103 elements.
                       array: array with CFID chemical descriptor names.
                       See: jarvis/core/specie.py

             enforce_c_size: minimum size of the simulation cell in Angst.
        """
        if get_prim:
            atoms = atoms.get_primitive_atoms
        dim = get_supercell_dims(atoms=atoms, enforce_c_size=enforce_c_size)
        atoms = atoms.make_supercell(dim)

        adj = np.array(atoms.raw_distance_matrix.copy())

        # zero out edges with bond length greater than threshold
        adj[adj >= max_cut] = 0

        if zero_diag:
            np.fill_diagonal(adj, 0.0)
        nodes = np.arange(atoms.num_atoms)
        if features == "atomic_number":
            node_attributes = np.array(
                [[np.array(Specie(i).Z)] for i in atoms.elements],
                dtype="float",
            )
        if features == "atomic_fraction":
            node_attributes = []
            fracs = atoms.composition.atomic_fraction_array
            for i in fracs:
                node_attributes.append(np.array([float(i)]))
            node_attributes = np.array(node_attributes)

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
            print("Please check the input options.")
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

        # construct edge list
        uv = []
        edge_features = []
        for ii, i in enumerate(atoms.elements):
            for jj, j in enumerate(atoms.elements):
                bondlength = adj[ii, jj]
                if bondlength > 0:
                    uv.append((ii, jj))
                    edge_features.append(bondlength)

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
        A = np.zeros((self.num_nodes, self.num_nodes))
        for edge, a in zip(self.edges, self.edge_attributes):
            A[edge] = a
        return A


class Standardize(torch.nn.Module):
    """Standardize atom_features: subtract mean and divide by std."""

    def __init__(self, mean: torch.Tensor, std: torch.Tensor):
        """Register featurewise mean and standard deviation."""
        super().__init__()
        self.mean = mean
        self.std = std

    def forward(self, g: dgl.DGLGraph):
        """Apply standardization to atom_features."""
        g = g.local_var()
        h = g.ndata.pop("atom_features")
        g.ndata["atom_features"] = (h - self.mean) / self.std
        return g


class StructureDataset(torch.utils.data.Dataset):
    """Dataset of crystal DGLGraphs."""

    def __init__(
        self,
        structures,
        targets,
        ids=None,
        cutoff=8.0,
        maxrows=np.inf,
        atom_features="atomic_number",
        transform=None,
        enforce_undirected=False,
        max_neighbors=12,
    ):
        """Initialize the class."""
        self.graphs = []
        self.labels = []
        self.ids = []

        for idx, (structure, target, id) in enumerate(
            tqdm(zip(structures, targets, ids))
        ):

            if idx >= maxrows:
                break
            a = Atoms.from_dict(structure)
            g = Graph.atom_dgl_multigraph(
                a,
                atom_features=atom_features,
                enforce_undirected=enforce_undirected,
                cutoff=cutoff,
                max_neighbors=max_neighbors,
                id=id,
            )

            self.graphs.append(g)
            self.labels.append(target)
            self.ids.append(id)

        self.labels = torch.tensor(self.labels).type(torch.get_default_dtype())
        self.transform = transform

    def __len__(self):
        """Get length."""
        return self.labels.shape[0]

    def __getitem__(self, idx):
        """Get StructureDataset sample."""
        g = self.graphs[idx]
        label = self.labels[idx]

        if self.transform:
            g = self.transform(g)

        return g, label

    def setup_standardizer(self):
        """Atom-wise feature standardization transform."""
        x = torch.cat([g.ndata["atom_features"] for g in self.graphs])
        self.atom_feature_mean = x.mean(0)
        self.atom_feature_std = x.std(0)

        self.transform = Standardize(
            self.atom_feature_mean, self.atom_feature_std
        )

    @staticmethod
    def collate(samples: List[Tuple[dgl.DGLGraph, torch.Tensor]]):
        """Dataloader helper to batch graphs cross `samples`."""
        graphs, labels = map(list, zip(*samples))
        batched_graph = dgl.batch(graphs)
        return batched_graph, torch.tensor(labels)


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
