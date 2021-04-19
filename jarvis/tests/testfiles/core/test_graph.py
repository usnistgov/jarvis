from jarvis.core.graphs import StructureDataset, Graph
from jarvis.db.figshare import data
import os

test_pos = os.path.join(
    os.path.dirname(__file__),
    "POSCAR-JVASP-7577",
)


def test_graph():
    from jarvis.core.atoms import Atoms
    from jarvis.db.figshare import get_jid_data

    atoms = Atoms.from_poscar(test_pos)
    g = g = Graph.atom_dgl_multigraph(atoms=atoms, atom_features="cgcnn")
    atoms = Atoms.from_dict(get_jid_data("JVASP-664")["atoms"])
    feature_sets = ("atomic_number", "basic", "cfid", "cgcnn")
    for i in feature_sets:
        g = Graph.atom_dgl_multigraph(atoms=atoms, atom_features=i)
        g = Graph.atom_dgl_multigraph(atoms=atoms, atom_features=i)
        print(i, g)
    g = Graph.from_atoms(atoms=atoms, features="atomic_number")
    g = Graph.from_atoms(atoms=atoms, features="atomic_number")
    g = Graph.from_atoms(atoms=atoms, features="atomic_fraction")
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
    g = Graph.from_atoms(atoms=atoms, features="cfid", max_cut=10000)
    print(g)
    d = g.to_dict()
    g = Graph.from_dict(d)
    num_nodes = g.num_nodes
    num_edges = g.num_edges
    print(num_nodes, num_edges)
    assert num_nodes == 48
    assert num_edges == 2256
    assert (g.adjacency_matrix.shape) == (48, 48)


def test_dataset():
    d = data("dft_2d")
    x = []
    y = []
    z = []
    for i in d[0:100]:
        if i["formation_energy_peratom"] != "na":
            x.append(i["atoms"])
            y.append(i["formation_energy_peratom"])
            z.append(i["jid"])
    s = StructureDataset(x, y, ids=z)
    col = s.collate


# test_graph()
