from pymatgen.io.vasp.inputs import Poscar
from jarvis.vasp.joptb88vdw import smart_converge
p=Poscar.from_file('POSCAR')
smart_converge(mat=p)


