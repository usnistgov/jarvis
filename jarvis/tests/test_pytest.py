import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.vasp.joptb88vdw import Auto_Kpoints
from jarvis.lammps.jlammps import read_data
from jarvis.sklearn.get_desc import get_comp_descp
from jarvis.lammps.Surf_Def import vac_antisite_def_struct_gen

def test_read_poscar():
    poscar=Structure.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
    assert len(poscar)== 2

def test_Auto_kpoints():
    poscar=Poscar.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
    kp=list(Auto_Kpoints(mat=poscar,length=20).kpts[0])
    assert kp == [6,6,6]

def test_read_poscar():
    poscar=Structure.from_file(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/POSCAR'))
    assert len(poscar)== 1

def test_read_data():
    dat=(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/data'))
    ff=str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/potential.mod')
    data= (read_data(data=dat,ff=ff))
    assert len(data)== 1

def test_vac_antisite_def_struct_gen():
    dat=(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/data'))
    ff=str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/potential.mod')
    data= (read_data(data=dat,ff=ff))
    vacs=vac_antisite_def_struct_gen(struct=data,c_size=0)
    assert len(vacs)== 2

def test_read_poscar():
    poscar=Structure.from_file(str('../sklearn/examples/POSCAR'))
    assert len(poscar) == 6
def test_desc():
    poscar=Structure.from_file(str('../sklearn/examples/POSCAR'))
    desc=get_comp_descp(poscar)
    assert len(desc)== 1557

