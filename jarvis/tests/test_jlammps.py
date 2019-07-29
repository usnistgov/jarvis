import os, builtins, io, pytest
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Poscar
from jarvis.lammps.jlammps import *
from jarvis.plane_defect.surface import surfer
from jarvis.point_defect.vacancy import vac_antisite_def_struct_gen
import pytest

parameters = {
    "exec": "mpirun /cluster/bin/lmp_ctcms-14439-knc6 <in.elastic >out",
    "pair_coeff": "/users/knc6/Software/jarvis/jarvis/lammps/examples/Mishin-Ni-Al-2009.eam.alloy",
    "control_file": "/users/knc6/inelast.mod",
    "pair_style": "eam/alloy",
    "atom_style": "charge",
    "cluster": "pbs",
}


init = os.path.join(
    os.path.dirname(__file__),
    "..",
    "lammps",
    "examples",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold",
    "mp-134",
    "init.mod",
)


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
log = os.path.join(
    os.path.dirname(__file__),
    "..",
    "lammps",
    "examples",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold",
    "mp-134",
    "log.lammps",
)

def test_read_data():
    data = read_data(data=dat, ff=ff)
    assert len(data) == 1

def test_write_lammps_in():
    data = read_data(data=dat, ff=ff)
    fold = os.path.join(os.path.dirname(__file__), "tmp")
    if not os.path.exists(fold):
        os.makedirs(fold)
    file = os.path.join(os.path.dirname(__file__), fold, "tmp")
    write_lammps_in(
        structure=data,
        lammps_in=init,
        lammps_in1=ff,
        lammps_in2=file,
        parameters=parameters,
    )

def test_vac_antisite_def_struct_gen(tmpdir):
    data = read_data(data=dat, ff=ff)
    vacs = vac_antisite_def_struct_gen(struct=data, c_size=0, write_file=False)
    assert len(vacs) == 2

def sample_strt():
    data = read_data(data=dat, ff=ff)
    return data


# def test_get_get_phonopy_atoms():
#    s=sample_strt()
#    phn=get_phonopy_atoms(mat=s)
#    assert phn.symbols==['Al']

def test_write_lammps_data():
    success = False
    s = sample_strt()
    fold = dat = os.path.join(os.path.dirname(__file__), "tmp")
    if not os.path.exists(fold):
        os.makedirs(fold)
    file = dat = os.path.join(os.path.dirname(__file__), fold, "tmp")
    write_lammps_data(structure=s, file=file, write_tmp_file=False)
    success = True
    assert success == True

def test_analyz_loge():
    x = len(analyz_loge(log))
    assert x == 23


# def test_smart_surf():
#    s=sample_strt()
#    parameters = {'phonon_control_file':'/users/knc6/in.phonon','surf_control_file':'/users/knc6/inelast_nobox.mod','def_control_file':'/users/knc6/inelast_nobox.mod','json_dat':'/rk2/knc6/JARVIS-DFT/MDCS/all_mp.json','c_size':3,'vac_size':25,'surf_size':25,'phon_size':20,'cluster':'pbs','exec':'/users/knc6/Software/LAMMPS/lammps-master/src/lmp_serial <in.main >out','pair_style':'eam/alloy','pair_coeff':'../lammps/examples/Al03.eam.alloy','atom_style': 'charge' ,'control_file':'/users/knc6/inelast.mod'}
#    cwd=str(os.getcwd())
#    os.chdir('../lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold')
#    sl,sh=smart_surf(strt=s,parameters=parameters)
#    print (sl,sh)
# test_smart_surf()
