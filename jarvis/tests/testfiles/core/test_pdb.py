from jarvis.core.pdb_atoms import read_pdb
import os

pdb_path = os.path.join(os.path.dirname(__file__), ".", "pdb101d.ent")


def test_pdb():
    pdb = read_pdb(pdb_path)
    # print (pdb.num_atoms,round(pdb.density,2))
    assert (pdb.num_atoms, round(pdb.density, 2)) == (448, 0.09)


# test_pdb()
