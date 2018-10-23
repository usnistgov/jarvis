import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.sklearn.get_desc import get_comp_descp

def test_desc():
    poscar=Structure.from_file(str('../sklearn/examples/POSCAR'))
    desc=get_comp_descp(poscar)
    assert len(desc)== 1557

