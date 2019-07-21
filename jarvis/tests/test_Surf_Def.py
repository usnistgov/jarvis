import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.lammps.jlammps import read_data
from jarvis.plane_defect.surface import surfer, pmg_surfer
from jarvis.point_defect.vacancy import vac_antisite_def_struct_gen

dat = os.path.join(
    os.path.dirname(__file__),
    "..",
    "lammps",
    "examples",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold",
    "mp-134",
    "data",
)
ff = os.path.join(
    os.path.dirname(__file__),
    "..",
    "lammps",
    "examples",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold",
    "mp-134",
    "potential.mod",
)


def sample_strt():
    data = read_data(data=dat, ff=ff)
    return data


def test_vac_antisite_def_struct_gen():
    data = read_data(data=dat, ff=ff)
    vacs = vac_antisite_def_struct_gen(struct=data, c_size=0, write_file=False)
    assert len(vacs) == 2


def test_pmg_surfer():
    s = sample_strt()
    surf = len(pmg_surfer(mat=s, write_file=False))
    assert surf == 4
