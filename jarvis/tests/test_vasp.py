import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.vasp.joptb88vdw import *  # Auto_Kpoints
from jarvis.lammps.jlammps import read_data
from jarvis.sklearn.get_desc import get_comp_descp
from jarvis.tools.vasp import *
from jarvis.slme.slme import *
from jarvis.boltztrap.boltztrap import *
from jarvis.elastic_tens.vasp import *
from jarvis.ip_optics.freq_dielectric import *
from jarvis.ip_optics.freq_dielectric import ip_optics

pos = os.path.join(
    os.path.dirname(__file__), "..", "vasp", "examples", "SiOptb88", "POSCAR"
)
dir = os.path.join(os.path.dirname(__file__), "..", "vasp", "examples", "SiOptb88")
run = os.path.join(
    os.path.dirname(__file__),
    "..",
    "vasp",
    "examples",
    "SiOptb88",
    "MAIN-BAND-bulk@mp_149",
    "vasprun.xml",
)
kpfile = os.path.join(
    os.path.dirname(__file__),
    "..",
    "vasp",
    "examples",
    "SiOptb88",
    "MAIN-BAND-bulk@mp_149",
    "KPOINTS",
)
out = os.path.join(
    os.path.dirname(__file__),
    "..",
    "vasp",
    "examples",
    "SiOptb88",
    "MAIN-ELASTIC-bulk@mp_149",
    "OUTCAR",
)
mbjrun = os.path.join(
    os.path.dirname(__file__),
    "..",
    "vasp",
    "examples",
    "SiOptb88",
    "MAIN-MBJ-bulk@mp_149",
    "vasprun.xml",
)
mainrun = os.path.join(
    os.path.dirname(__file__),
    "..",
    "vasp",
    "examples",
    "SiOptb88",
    "MAIN-RELAX-bulk@mp_149",
    "vasprun.xml",
)


def sample_structure():
    poscar = Poscar.from_file(pos)
    return poscar


def test_read_poscar():
    poscar = Structure.from_file(pos)
    assert len(poscar) == 2


def test_Auto_kpoints():
    poscar = Poscar.from_file(pos)
    kp = list(Auto_Kpoints(mat=poscar, length=20).kpts[0])
    assert kp == [6, 6, 6]


def test_check_polar():
    p = sample_structure().structure
    surfaces = surfer(mat=p, write_file=False)
    s = (surfaces[0]).structure
    print(check_polar(s))
    assert check_polar(s) == False


# def test_get_lowest_en_from_mp():
#    en_atom=get_lowest_en_from_mp('Al')
#    assert round(float(en_atom),4) == round(float(-3.74810318),4)


def test_make_big():
    p = sample_structure()
    strt = p.structure.copy()
    big_p = make_big(poscar=p)
    assert len(strt) * 3 ** 3 == len(big_p.structure)


def test_converg_encut():
    p = sample_structure()
    cwd = str(os.getcwd())
    os.chdir(dir)
    final_enc = converg_encut(mat=p)
    os.chdir(cwd)
    assert final_enc == 500


def test_converg_kpoints():
    p = sample_structure()
    cwd = str(os.getcwd())
    os.chdir(dir)
    final_kp = converg_kpoints(mat=p)
    os.chdir(cwd)
    assert final_kp == 35


def test_smart_converge():
    p = sample_structure()
    cwd = str(os.getcwd())
    os.chdir(dir)
    en, mat_f = smart_converge(mat=p)
    assert float(en) == float(-8.3381172)


def test_bandgap1():
    eg, dir = bandgap(vrun=run)
    assert eg == 0.006000000000000227


def test_bandgap2():
    vrun = str(run).replace("BAND", "RELAX")
    eg, dir = bandgap(vrun=vrun)
    # print (eg,dir)
    assert eg == 0.14800000000000058


def test_bandstr():
    bplot = bandstr(vrun=run, kpfile=kpfile)
    # bplot.savefig('pp.png')
    assert (bplot.xlim()[0]) == 0
    bplot.close()


def test_plot_dos():
    plt1, plt2, plt3 = plot_dos(vrun=run)
    assert (plt1.xlim()[0]) == -5
    plt1.close()
    plt2.close()
    plt3.close()


# def test_plot_kp_convergence():
#    _,kp=plot_kp_convergence(dir)
#    assert kp=='11x11x11'


def test_plot_enc_convergence():
    _, enc = plot_enc_convergence(dir)
    assert enc == 500


def test_ip_optics():
    vrun = str(run).replace("BAND", "OPTICS")
    op = ip_optics(vrun)
    val = op["real_part"][0][0]
    assert val == 13.4865


def test_elastic_props():
    info = elastic_props(outcar=out)
    assert (info["KV"]) == 87.267


def test_magnetic_moment():
    osz = str(run).replace("BAND", "RELAX").replace("vasprun.xml", "OSZICAR")
    m = magnetic_moment(osz)
    assert m == 0


def test_get_spacegroup():
    strt = Structure.from_file(pos)
    num, symb = get_spacegroup(strt)
    assert symb == "Fd-3m"


def test_get_area():
    strt = Structure.from_file(pos)
    ar = get_area(strt)
    assert ar == 12.950104933895442


def test_slme():
    en, abz, dirgap, indirgap = optics(mbjrun)
    abz = abz * 100.0
    SLME = slme(en, abz, indirgap, indirgap, plot_current_voltage=False)
    SQ = calculate_SQ(indirgap)
    assert SLME == 0.3323125002699776
def test_ip_opt():
  ip = ip_optics(mbjrun)  
  assert ip["absorption"][0][0]==0
# def test_boltztrap():
#   b=boltz_run(mainrun)
#   val=get_prop(b,prop='zt')
#   assert val[0]==1.5913528701329165
