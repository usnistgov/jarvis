from jarvis.core.atoms import Atoms
import spglib
from jarvis.core.specie import Specie

class Spacegroup3D(object):
    def __init__(self,atoms=[],dataset={},symprec=1e-2,angle_tolerance=5):
        self._dataset = dataset
        self._atoms=atoms
        self._symprec=symprec
        self._angle_tolerance=angle_tolerance

    def spacegroup_data(self):
            phonopy_atoms = (self._atoms.lattice_mat, self._atoms.frac_coords, self._atoms.atomic_numbers)
            dataset = spglib.get_symmetry_dataset(
            phonopy_atoms, symprec=self._symprec, angle_tolerance=self._angle_tolerance)
            """    
            keys = ('number',
            'hall_number',
            'international',
            'hall',
            'choice',
            'transformation_matrix',
            'origin_shift',
            'rotations',
            'translations',
            'wyckoffs',
            'site_symmetry_symbols',
            'equivalent_atoms',
            'mapping_to_primitive',
            'std_lattice',
            'std_types',
            'std_positions',
            'std_rotation_matrix',
            'std_mapping_to_primitive',
            'pointgroup')
            """
            return Spacegroup3D( atoms=self._atoms,symprec=self._symprec, angle_tolerance=self._angle_tolerance,dataset=dataset)
    @property
    def get_space_group_symbol(self):
         spg = self.spacegroup_data()
         return spg.spacegroup_data()._dataset["international"]


    @property
    def get_space_group_number(self):
         spg = self.spacegroup_data()
         return spg.spacegroup_data()._dataset["number"]

    @property
    def primitive_atoms(self):
         phonopy_atoms = (self._atoms.lattice_mat, self._atoms.frac_coords, self._atoms.atomic_numbers)
         lattice, scaled_positions, numbers = spglib.find_primitive(phonopy_atoms,symprec=self._symprec)
         elements = self._atoms.elements
         el_dict = {}
         for i in elements:
           el_dict.setdefault(Specie(i).Z,i)     
         prim_elements =  [el_dict[i] for i in numbers]
         prim_atoms=Atoms(lattice_mat=lattice,elements=prim_elements,coords=scaled_positions,cartesian=False)
         return prim_atoms

    @property
    def get_crystal_system(self):

        spg = self.spacegroup_data()
        n = spg.spacegroup_data()._dataset["number"]

        def f(i, j):
            return i <= n <= j

        cs = {"triclinic": (1, 2), "monoclinic": (3, 15),
              "orthorhombic": (16, 74), "tetragonal": (75, 142),
              "trigonal": (143, 167), "hexagonal": (168, 194),
              "cubic": (195, 230)}

        crystal_sytem = None

        for k, v in cs.items():
            if f(*v):
                crystal_sytem = k
                break
        return crystal_sytem


if __name__=='__main__':
   box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
   coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
   elements = ["Si", "Si"]
   Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
   spg=Spacegroup3D(atoms=Si)#.spacegroup_data()
   print (spg.get_space_group_symbol)
   print (spg.get_space_group_number)
   primt = spg.primitive_atoms
   print ('primt',primt)
   print ('cryst_sys',spg.get_crystal_system)
   #pmg=Si.pymatgen_converter()
   #from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
   #spg=SpacegroupAnalyzer(pmg)
   #print (spg.get_space_group_symbol(),spg.get_space_group_number())
   #print (pmg.get_primitive_structure())
