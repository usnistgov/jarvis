import os, builtins, io, pytest
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar

# from jarvis.lammps.jlammps import *
from jarvis.plane_defect.surface import surfer
from jarvis.point_defect.vacancy import vac_antisite_def_struct_gen

from jarvis.heterostruct.interface import *

pos = os.path.join(
    os.path.dirname(__file__), "..", "vasp", "examples", "SiOptb88", "POSCAR"
)


def test_get_ase_surf():
    s = Structure.from_file(pos)
    surf = get_ase_surf(s)
    assert (surf.lattice.matrix[2][2]) == 24.78934299449321


def test_mismatch_strts():
    s = Structure.from_file(pos)
    surf = get_ase_surf(s)
    info = mismatch_strts(film=surf, subs=surf)
    assert info["area1"] == 12.95010493389544


def test_get_hetero_type():
    A = {}
    B = {}
    A["scf_vbm"] = 1.0
    A["scf_cbm"] = 3.0
    A["avg_max"] = 2.5

    B["scf_vbm"] = 1.5
    B["scf_cbm"] = 2.0
    B["avg_max"] = 4.5
    int_type, stack = get_hetero_type(A=A, B=B)
    assert (int_type, stack) == ("III", "AB")


def test_get_direct_hetero():
    s = Structure.from_file(pos)
    x = get_direct_hetero(bulk_film=s, bulk_subs=s)
    print(len(x))
    assert len(x) == 24
