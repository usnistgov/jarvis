import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.vasp.joptb88vdw import  * #Auto_Kpoints
from jarvis.lammps.jlammps import read_data
from jarvis.sklearn.get_desc import get_comp_descp
from jarvis.lammps.Surf_Def import vac_antisite_def_struct_gen,surfer
from jarvis.tools.vasp import *

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


def test_bandgap1():
    vrun='../vasp/examples/SiOptb88/MAIN-BAND-bulk@mp_149/vasprun.xml'
    eg,dir=bandgap(vrun=vrun)
    assert eg==0.006000000000000227

def test_bandgap2():
    vrun='../vasp/examples/SiOptb88/MAIN-RELAX-bulk@mp_149/vasprun.xml'
    eg,dir=bandgap(vrun=vrun)
    #print (eg,dir)
    assert eg==0.14800000000000058

def test_bandstr():
    vrun='../vasp/examples/SiOptb88/MAIN-BAND-bulk@mp_149/vasprun.xml'
    kpfile='../vasp/examples/SiOptb88/MAIN-BAND-bulk@mp_149/KPOINTS'
    bplot=bandstr(vrun=vrun,kpfile=kpfile)
    #bplot.savefig('pp.png')
    assert (bplot.xlim()[0])==0
    bplot.close()

def test_plot_dos():
   vrun='../vasp/examples/SiOptb88/MAIN-BAND-bulk@mp_149/vasprun.xml'
   plt1,plt2,plt3= plot_dos(vrun=vrun)
   assert (plt1.xlim())==-5
   plt1.close()
   plt2.close()
   plt3.close()

def test_plot_kp_convergence():
    dir='../vasp/examples/SiOptb88'
    _,kp=plot_kp_convergence(dir)
    assert kp=='19x19x19'

def test_plot_enc_convergence():
    dir='../vasp/examples/SiOptb88'
    _,enc=plot_enc_convergence(dir)
    assert enc==500

def test_IP_optics():
   vrun='../vasp/examples/SiOptb88/MAIN-OPTICS-bulk@mp_149/vasprun.xml'
   op=optics(vrun)
   val=op['real_part'][0][0]
   assert val==13.4865

def test_elastic_props():
   out='../vasp/examples/SiOptb88/MAIN-ELASTIC-bulk@mp_149/OUTCAR'
   info=elastic_props(out)
   assert  (info['KV'])==87.267

def test_magnetic_moment():
   osz='../vasp/examples/SiOptb88/MAIN-ELASTIC-bulk@mp_149/OSZICAR'
   m=magnetic_moment(osz)
   assert m==0

def test_get_spacegroup():
   strt=Structure.from_file('../vasp/examples/SiOptb88/MAIN-ELASTIC-bulk@mp_149/POSCAR')
   num,symb=get_spacegroup(strt)
   assert symb=='Fd-3m'

def test_get_aea():
   strt=Structure.from_file('../vasp/examples/SiOptb88/MAIN-ELASTIC-bulk@mp_149/POSCAR')
   ar=get_area(strt)
   assert ar==30.180025513225004

