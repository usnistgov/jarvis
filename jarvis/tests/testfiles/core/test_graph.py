from jarvis.core.graphs import Graph


def test_graph():
    from jarvis.core.atoms import Atoms
    from jarvis.db.figshare import get_jid_data

    atoms = Atoms.from_dict(get_jid_data("JVASP-664")["atoms"])
    g = Graph.from_atoms(atoms=atoms, features="atomic_number")
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
    print(g)
    d = g.to_dict()
    g = Graph.from_dict(d)
    num_nodes = g.num_nodes
    num_edges = g.num_edges
    print(num_nodes, num_edges)
    assert num_nodes == 48
    assert num_edges == 2304
    assert (g.adjacency_matrix.shape) == (48,48)
