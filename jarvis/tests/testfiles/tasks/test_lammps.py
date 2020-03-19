import os

lammps_path = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "lammps",
    "Al03.eam.alloy_nist",
    "bulk@mp-134_fold"
)

cwd = str(os.getcwd())

def test_lammps():
    os.chdir(lammps_path)
    cmd = str('python master.py')
    os.system(cmd)

os.chdir(cwd)    
