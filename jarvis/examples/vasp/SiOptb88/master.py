from jarvis.io.vasp.inputs import Poscar
from jarvis.tasks.vasp.vasp import JobFactory

p = Poscar.from_file("POSCAR")
print (p)
JobFactory().all_optb88vdw_props(mat=p)
