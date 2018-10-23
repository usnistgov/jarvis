import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.vasp.joptb88vdw import  * #Auto_Kpoints
from jarvis.lammps.jlammps import read_data
from jarvis.sklearn.get_desc import get_comp_descp
from jarvis.lammps.Surf_Def import vac_antisite_def_struct_gen,surfer

def sample_structure():
    poscar=Poscar.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
    return poscar

def test_read_poscar():
    poscar=Structure.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
    assert len(poscar)== 2

def test_Auto_kpoints():
    poscar=Poscar.from_file(str('../vasp/examples/SiOptb88/POSCAR'))
    kp=list(Auto_Kpoints(mat=poscar,length=20).kpts[0])
    assert kp == [6,6,6]

def test_check_polar():
    p=sample_structure().structure
    surfaces=surfer(mat=p)
    s=(surfaces[0]).structure
    print (check_polar(s))
    assert check_polar(s)==False

def test_get_lowest_en_from_mp():
    en_atom=get_lowest_en_from_mp('Al')
    assert round(float(en_atom),4) == round(float(-3.74810318),4)

def test_make_big():
    p=sample_structure()
    strt=p.structure.copy()
    big_p=make_big(poscar=p)
    assert len(strt)*3**3==len(big_p.structure)

def test_converg_encut():
    p=sample_structure()
    cwd=str(os.getcwd())
    os.chdir('../vasp/examples/SiOptb88')
    final_enc=converg_encut(mat=p)    
    os.chdir(cwd)
    assert final_enc == 500

def test_converg_kpoints():
    p=sample_structure()
    cwd=str(os.getcwd())
    os.chdir('../vasp/examples/SiOptb88')
    final_kp=converg_kpoints(mat=p)    
    os.chdir(cwd)
    assert final_kp == 35

def test_smart_converge():
    p=sample_structure()
    cwd=str(os.getcwd())
    os.chdir('../vasp/examples/SiOptb88')
    en,mat_f=smart_converge(mat=p)
    assert float(en)==float(-8.3381172)

