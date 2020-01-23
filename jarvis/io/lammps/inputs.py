from jarvis.io.vasp.inputs import Poscar
from jarvis.core.atoms import Atoms
import numpy as np
from jarvis.analysis.structure.spacegroup import Spacegroup3D
class LammpsData(object):
     def __init__(self,atoms=None, lammps_box=[],species=[],charges=[],cart_coords=[],element_order = []):

        self._lammps_box = lammps_box
        self._species = species
        self._charges = charges
        self._cart_coords = cart_coords
        self._atoms = atoms

     def atoms_to_lammps(self,origin=(0,0,0), write_to_file='lmp'):
         atoms = self._atoms
         elements = atoms.elements
         num_atoms = atoms.num_atoms
         a,b,c = atoms.lattice.abc
         xlo, ylo, zlo = origin
         xhi = a + xlo
         m = atoms.lattice._lat
         xy = np.dot(m[1], m[0] / a)
         yhi = np.sqrt(b ** 2 - xy ** 2) + ylo
         xz = np.dot(m[2], m[0] / a)
         yz = (np.dot(m[1], m[2]) - xy * xz) / (yhi - ylo)
         zhi = np.sqrt(c ** 2 - xz ** 2 - yz ** 2) + zlo
         tilt =  [xy, xz, yz]
         rot_matrix = np.linalg.solve([[xhi - xlo, 0, 0],
                                  [xy, yhi - ylo, 0],
                                  [xz, yz, zhi - zlo]], m)
         box = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
         new_coords =  (np.dot(rot_matrix, atoms.cart_coords.T)).T
         #new_coords =  (np.dot(rot_matrix, atoms.frac_coords.T)).T
         unique_elements = list(set(elements))
         n_atom_types = len(unique_elements)
         n_atoms = len(elements)
         f = open(write_to_file, "w")
         f.write("datafile (written by JARVIS-FF) \n\n")
         f.write("%d \t atoms \n" % n_atoms)
         #print ('unique_elements', unique_elements)
         species = []
         if self._charges == []:
            self._charges=np.zeros(len(elements))
         f.write("%d  atom types\n" % n_atom_types)
         f.write("0.0 %s  xlo xhi\n" % xhi)
         f.write("0.0 %s  ylo yhi\n" % yhi)
         f.write("0.0 %s  zlo zhi\n" % zhi)
         f.write("%s %s %s  xy xz yz\n" % (xy, xz, yz))
         f.write("\n\n")
         f.write("Atoms \n\n")
         #print ('lens',len(elements),len(new_coords),new_coords.shape)
         for i,ii in enumerate(elements):
            s = unique_elements.index(ii)+1
            r = new_coords[i]
            charge = self._charges[i]
            #print (i+1,s,charge,r)
            f.write("%6d %3d %6f %s %s %s\n" % (i + 1, s, charge, r[0], r[1], r[2]))
            #f.write("%6d %3d %6f %s %s %s\n" % ((i + 1, s, charge) + tuple(r)))
         #return LammpsData(box=box,species=


         f.close()
         alpha, beta, gamma = atoms.lattice.angles
         print (alpha, beta, gamma,a,b,c)
         h11 = a
         h12 = 0.0
         h13 = 0.0

         h21 = np.cos(gamma) * b
         h22 = np.sin(gamma)*b #np.sqrt(b**2-h21**2)#np.sin(gamma) * b #np.sqrt(b**2-h21**2)
         h23 = 0.0

         h31 = np.cos(beta) * c
         h32 = (b * c * np.cos(alpha) - h31 * h21) / float(h22)
         h33 = (c ** 2 - h31 ** 2 - h32 ** 2)
         #print ('h11,h12,h13,h21,h22,h23,h31,h32,h33',h11,h12,h13,h21,h22,h23,h31,h32,h33)

if __name__ == "__main__":
   box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
   coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
   elements = ["Si", "Si"]
   Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
   spg=Spacegroup3D(atoms=Si) 
   Si = spg.conventional_standard_structure
   p =Poscar.from_file('/rk2/knc6/JARVIS-FF/FS/Al1.eam.fs_nist/bulk@mp-134_fold/mp-134/new_pymatgen_slab.vasp'  )
   p =Poscar.from_file('/rk2/knc6/JARVIS-FF/COMB/ffield.comb3.NiAlO_nist/bulk@mp-1143_fold/bulk@mp-1143/new_pymatgen_slab.vasp'  )
   lmp = LammpsData(atoms=p.atoms)
   lat = lmp.atoms_to_lammps()
   #print (p)
