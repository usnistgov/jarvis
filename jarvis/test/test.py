import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.vasp.joptb88vdw import Auto_Kpoints
from jarvis.lammps.jlammps import read_data
from jarvis.sklearn.get_desc import get_comp_descp

class jvasp(unittest.TestCase):
    def test_read_poscar(self):
        poscar=Structure.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
        self.assertEqual(len(poscar), 2)

    def test_auto_kpoints(self):
        poscar=Poscar.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
        kp=list(Auto_Kpoints(mat=poscar,length=20).kpts[0])
        self.assertEqual(kp, [6,6,6])

class jlammps(unittest.TestCase):
    def test_read_poscar(self):
        poscar=Structure.from_file(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/POSCAR'))
        self.assertEqual(len(poscar), 1)

    def test_read_data(self):
        dat=(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/data'))
        ff=str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/potential.mod')
        data= (read_data(data=dat,ff=ff))
        self.assertEqual(len(data), 1)

class jsklearn(unittest.TestCase):
    def test_read_poscar(self):
        poscar=Structure.from_file(str('../sklearn/examples/POSCAR'))
        self.assertEqual(len(poscar), 6)
    def test_desc(self):
        poscar=Structure.from_file(str('../sklearn/examples/POSCAR'))
        desc=get_comp_descp(poscar)
        self.assertEqual(len(desc), 1557)

if __name__ == '__main__':
    unittest.main()
