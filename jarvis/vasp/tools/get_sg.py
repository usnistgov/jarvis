##http://pymatgen.org/_static/Basic%20functionality.html
import pymatgen as mg
from pymatgen.io.aseio import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
#si = mg.Element("Si")
#print("Atomic mass of Si is {}".format(si.atomic_mass))
#print("Si has a melting point of {}".format(si.melting_point))
#print("Ionic radii for Si: {}".format(si.ionic_radii))
######################################################
structure = mg.Structure.from_file("POSCAR")
finder = SpacegroupAnalyzer(structure)
num=finder.get_spacegroup_number()
print(num)
ase_atoms = AseAtomsAdaptor().get_atoms(structure)
import ase.io.vasp
ase.io.vasp.write_vasp('ase_atoms',ase_atoms,direct=False)
