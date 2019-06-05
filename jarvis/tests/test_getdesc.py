import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.sklearn.get_desc import get_comp_descp


def test_desc():
    fname = os.path.join(
        os.path.dirname(__file__), "..", "sklearn", "examples", "POSCAR"
    )
    poscar = Structure.from_file(fname)
    desc = get_comp_descp(poscar)
    assert len(desc) == 1557
