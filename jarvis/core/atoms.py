#https://github.com/acadien/matcalc/blob/037d7080c07856877d7a3f5f9dcbb2dec5f38dd1/analysis/rdf.py
#/rk2/knc6/UFHPC_3_16_2016/scratch-uf/LOREN/spicykamal/TaoLammps/COMB3_v18/src/analysis_tools.f
"""
This module provides classes to specify atomic structure
"""
from collections import Counter
import numpy as np
from jarvis.core.composition import Composition
from jarvis.core.specie import Specie
from jarvis.core.lattice import Lattice
import importlib.util
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import math
amu_gm = 1.66054e-24
ang_cm = 1e-8


class Atoms(object):
    def __init__(
        self, lattice_mat=None, coords=None, elements=None, cartesian: bool = False
    ):
        """
        Create atomic structure with lattice, coordinates, atom type and other information
        >>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
        >>> coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
        >>> elements = ["Si", "Si"]
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
        >>> print(round(Si.volume,2))
        40.03
        >>> Si.composition
        {'Si': 2}
        >>> round(Si.density,2)
        2.33
        >>> round(Si.packing_fraction,2)
        0.28
        >>> Si.atomic_numbers
        [14, 14]
        >>> Si.num_atoms
        2
        >>> Si.frac_coords[0][0]
        0
        >>> Si.cart_coords[0][0]
        0.0
        >>> coords = [[0, 0, 0], [1.3575 , 1.22175, 1.22175]]
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements,cartesian=True)
        >>> round(Si.density,2)
        2.33
        >>> Si.spacegroup()
        'C2/m (12)'
        >>> Si.pymatgen_converter()!={}
        True
        """

        self.lattice_mat = np.array(lattice_mat)
        self.lattice = Lattice(lattice_mat)
        self.coords = coords
        self.elements = elements
        if cartesian == True:
            self.cart_coords = self.coords
            self.frac_coords = self.lattice.frac_coords(self.coords)
            #print ('TRUE')
        else:
            self.frac_coords = self.coords
            self.cart_coords = self.lattice.cart_coords(self.coords)
            #print ('FALSE')

    @property
    def volume(self):
        m = self.lattice_mat
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol

    @property
    def composition(self):
        comp = {}
        for i in self.elements:
            comp[i] = comp.setdefault(i, 0) + 1
        return Composition(comp)

    @property
    def density(self):
        den = float(self.composition.weight * amu_gm) / (
            float(self.volume) * (ang_cm) ** 3
        )
        return den

    @property
    def atomic_numbers(self):
        numbers = []
        for i in self.elements:
            numbers.append(Specie(i).Z)
        return numbers

    @property
    def num_atoms(self):
        return len(self.coords)

    def pymatgen_converter(self):
        try:
            from pymatgen.core.structure import Structure

            return Structure(
                self.lattice_mat,
                self.elements,
                self.frac_coords,
                coords_are_cartesian=False,
            ) 
        except:pass

    def spacegroup(self, symprec=1e-3):
        #try:
            import spglib
            sg = spglib.get_spacegroup(
                (self.lattice_mat, self.frac_coords, self.atomic_numbers),
                symprec=symprec,
            )
            return sg
        #except:
        #    pass
    @property
    def packing_fraction(self):
        total_rad = 0
        for i in self.elements:
            total_rad = total_rad + Specie(i).atomic_rad**3
        pf = np.array([4 * np.pi * total_rad / (3 * self.volume)])
        return round(pf[0],5)
    #def write(self,filename='POSCAR;,format='poscar'):

    def make_supercell(self,dim=[2,2,2]):
     dim = np.array(dim)
     if dim.shape==(3,3):
          dim= np.array([int(np.linalg.norm(v)) for v in dim])
     coords = self.frac_coords
     all_symbs = self.elements #[i.symbol for i in s.species]
     nat = len(coords)

     new_nat = nat * dim[0] * dim[1] * dim[2]
     print ('new_nat,dim',new_nat,dim)
     new_coords = np.zeros((new_nat, 3))
     new_symbs = []  # np.chararray((new_nat))

     count = 0
     for i in range(nat):
        for j in range(dim[0]):
            for k in range(dim[1]):
                for l in range(dim[2]):
                    new_coords[count][0] = (coords[i][0] + j) / float(dim[0])
                    new_coords[count][1] = (coords[i][1] + k) / float(dim[1])
                    new_coords[count][2] = (coords[i][2] + l) / float(dim[2])
                    new_symbs.append(all_symbs[i])
                    count = count + 1

     nat = new_nat

     nat = len(coords)  # int(s.composition.num_atoms)
     lat = np.zeros((3, 3))
     box = self.lattice_mat
     lat[0][0] = dim[0] * box[0][0]
     lat[0][1] = dim[0] * box[0][1]
     lat[0][2] = dim[0] * box[0][2]
     lat[1][0] = dim[1] * box[1][0]
     lat[1][1] = dim[1] * box[1][1]
     lat[1][2] = dim[1] * box[1][2]
     lat[2][0] = dim[2] * box[2][0]
     lat[2][1] = dim[2] * box[2][1]
     lat[2][2] = dim[2] * box[2][2]

     super_cell=Atoms(lattice_mat=lat, coords=new_coords, elements=new_symbs, cartesian=False)
     return super_cell


    def __repr__(self):
        header= str('\nSystem\n1.0\n')+str(self.lattice_mat[0][0])+' '+str(self.lattice_mat[0][1])+' '+str(self.lattice_mat[0][2])+'\n'+ str(self.lattice_mat[1][0])+' '+str(self.lattice_mat[1][1])+' '+str(self.lattice_mat[1][2])+'\n'+str(self.lattice_mat[2][0])+' '+str(self.lattice_mat[2][1])+' '+str(self.lattice_mat[2][2])+'\n'
        middle = ' '.join(map(str,Counter(self.elements).keys()))+'\n'+' '.join(map(str,Counter(self.elements).values()))+'\ndirect\n'
        rest = ''
        for i in self.frac_coords:
            rest=rest+' '.join(map(str,i))+'\n'
        result = header+middle+rest
        return result
       
     
if __name__=='__main__':
   box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]] 
   coords = [[0, 0, 0], [0.25, 0.2, 0.25]] 
   elements = ["Si", "Si"]
   Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
   print ('pf',Si.packing_fraction,Si.make_supercell())
   pmg = Si.pymatgen_converter()
   pmg.make_supercell([2,2,2])
   print (pmg)
# if __name__=='__main__':
#    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]] 
#    coords = [[0, 0, 0], [0.25, 0.2, 0.25]] 
#    elements = ["Si", "Si"]
#    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
#    comp=Si.composition
#    #print (comp,Si.density)
#    print (Si.atomic_numbers)
#   print (Si.pymatgen_converter().composition.weight,Si.composition.weight,Si.density)
