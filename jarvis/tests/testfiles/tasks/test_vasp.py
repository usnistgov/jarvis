import os

vasp_path = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "..",
    "examples",
    "vasp",
    "SiOptb88"
)

cwd = str(os.getcwd())

def test_vasp():
    os.chdir(vasp_path)
    cmd = str('python master.py')
    os.system(cmd)

os.chdir(cwd)    
