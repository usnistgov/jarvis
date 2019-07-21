from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.io.vasp.inputs import Potcar
import yaml
import os

# os.environ["VASP_PSP_DIR"] = "/home/knc6/Software/VASP-POTENTIAL"
pot_yaml = str(
    os.path.join(os.path.dirname(__file__), "..", "vasp", "Special_POTCAR.yaml")
)
with open(pot_yaml, "r") as f:
    doc = yaml.load(f)
    pots = doc["POTCAR"]
# mat = Poscar.from_file("POSCAR")
# new_symb = []
# for el in mat.site_symbols:
#    new_symb.append(pots[el])
# potcar = Potcar(symbols=new_symb, functional="PBE")
# potcar.write_file("POTCAR")
