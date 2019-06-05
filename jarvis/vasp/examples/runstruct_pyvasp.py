from jarvis.vasp.joptb88vdw import smart_converge
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure

import os, time

s = Structure.from_file("POSCAR")
p = Poscar(s)
p.comment = "bulk@ATAT"
en, final = smart_converge(mat=p, elast_prop=False)
print("en,finel", en, final)
