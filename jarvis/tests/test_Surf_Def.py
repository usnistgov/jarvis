import unittest
import os
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.lammps.jlammps import read_data
from jarvis.lammps.Surf_Def import * #vac_antisite_def_struct_gen,surfer


def sample_strt():
    dat=(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/data'))
    ff=str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/potential.mod')
    data= (read_data(data=dat,ff=ff))
    return (data)

def test_vac_antisite_def_struct_gen():
    dat=(str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/data'))
    ff=str('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/mp-134/potential.mod')
    data= (read_data(data=dat,ff=ff))
    vacs=vac_antisite_def_struct_gen(struct=data,c_size=0)
    assert len(vacs)== 2

def test_pmg_surfer():
    s=sample_strt()
    surf=len(pmg_surfer(mat=s))
    assert surf==4

