from ase.lattice.surface import surface
from pymatgen.io.ase import AseAtomsAdaptor
from jarvis.io.vasp.inputs import Poscar
from jarvis.core.lattice import Lattice
from jarvis.core.atoms import Atoms
import numpy as np
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from numpy.linalg import norm, solve
from numpy import gcd 
def ext_gcd(a, b):
    if b == 0:
        return 1, 0
    elif a % b == 0:
        return 0, 1
    else:
        x, y = self.ext_gcd(b, a % b)
        return y, x - y * (a // b)


class Surface(object):
        def __init__(self,atoms=None, indices=[0,0,1],layers=3, vacuum=18.0, tol=1e-10,from_conventional_structure=True):
            self.indices =  np.array(indices)
            if from_conventional_structure:
                self.atoms = Spacegroup3D(atoms).conventional_standard_structure
            else:
                  self.atoms = atoms
            self.tol=tol
            self.vacuum=vacuum
            self.layers=layers

        def make_surface(self):
            atoms=self.atoms
            h, k, l = self.indices
            h0, k0, l0 = (self.indices == 0)

            if h0 and k0 or h0 and l0 or k0 and l0:  # if two indices are zero
                if not h0:
                    c1, c2, c3 = [(0, 1, 0), (0, 0, 1), (1, 0, 0)]
                if not k0:
                    c1, c2, c3 = [(0, 0, 1), (1, 0, 0), (0, 1, 0)]
                if not l0:
                    c1, c2, c3 = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
            else:
                p, q = ext_gcd(k, l)
                a1, a2, a3 = self.atoms.lattice_mat#.lat_lengths()

                # constants describing the dot product of basis c1 and c2:
                # dot(c1,c2) = k1+i*k2, i in Z
                k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3),
                            l * a2 - k * a3)
                k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3),
                            l * a2 - k * a3)

                if abs(k2) > self.tol:
                    i = -int(round(k1 / k2))  # i corresponding to the optimal basis
                    p, q = p + i * l, q - i * k

                a, b = ext_gcd(p * k + q * l, h)

                c1 = (p * k + q * l, -p * h, -q * h)
                c2 = np.array((0, l, -k)) // abs(gcd(l, k))
                c3 = (b, a * p, a * q)
            #print ('c1c2c3',c1,c2,c3)
            lattice = atoms.lattice_mat#.lat_lengths()
            basis = np.array([c1,c2,c3])
            print ('basis',basis)
            scaled = np.dot(basis, np.array(atoms.frac_coords).T).T
            #scaled = np.linalg.solve(basis, np.array(atoms.frac_coords).T).T
            scaled -= np.floor(scaled + self.tol)
            print ('scaled_tol',scaled)
            #atoms.frac_coords=scaled
            new_coords = scaled
            new_lattice = np.dot(basis, lattice)
            #print ('new_lattice',new_lattice)
            surf_atoms = Atoms(lattice_mat=new_lattice, elements=self.atoms.elements, coords=new_coords, cartesian=True)
            surf_atoms=surf_atoms.make_supercell([1,1,self.layers])
            print ('scaled_111',surf_atoms.frac_coords)
            new_lat = surf_atoms.lattice_mat#lat_lengths()
            a1=new_lat[0]
            a2=new_lat[1]
            a3=new_lat[2]
            new_lat = np.array([a1,a2,np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /norm(np.cross(a1, a2))**2] )
            print ('new_lat1',new_lat)
            #surf_atoms.lattice_mat = new_lat
            new_coords = np.dot(new_lat, np.array(surf_atoms.frac_coords.T)).T
            #new_coords = np.linalg.solve(new_lat.T, np.array(surf_atoms.frac_coords).T).T
            surf_atoms = Atoms(lattice_mat=new_lat, elements=surf_atoms.elements, coords=new_coords, cartesian=True)
            print ('scaled_222',surf_atoms.frac_coords)
            #print ('new_lat1',new_lat)
            #tmp  = np.array(new_lat, dtype=np.float64)#.reshape((3, 3))
            a1=new_lat[0]
            a2=new_lat[1]
            a3=new_lat[2]
            #print ('a1a2a3',a1,a2,a3)
            #print ('a3',np.linalg.norm(a3))
            new_lat = np.array([(np.linalg.norm(a1), 0, 0),(np.dot(a1, a2) / np.linalg.norm(a1),np.sqrt(np.linalg.norm(a2)**2 - (np.dot(a1, a2) / np.linalg.norm(a1))**2), 0),(0, 0, np.linalg.norm(a3))])
            #print ('new_lat2',new_lat)
            new_coords = np.dot(new_lat, np.array(surf_atoms.frac_coords.T)).T
            surf_atoms = Atoms(lattice_mat=new_lat, elements=surf_atoms.elements, coords=new_coords, cartesian=False)
            print ('scaled',surf_atoms.frac_coords)
            new_coords = surf_atoms.frac_coords
            #print ('scaled',new_coords)
            new_coords[:,:2]%=1 
            #print ('new_lat',new_lat)
            surf_atoms=Atoms(lattice_mat=new_lat, elements=surf_atoms.elements,coords=new_coords,cartesian=False).center()
            print (surf_atoms)


if __name__=='__main__':
   box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
   coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
   elements = ["Si", "Si"]
   Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
   #Si = Poscar.from_file('/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-134_bulk_PBEBO/MAIN-RELAX-bulk@mp-134/POSCAR').atoms
   tmp = Spacegroup3D(Si).conventional_standard_structure
   pmg = tmp.pymatgen_converter()
   ase_atoms = AseAtomsAdaptor().get_atoms(pmg)
   print ('cell=',ase_atoms.cell)
   ase_slab = surface(ase_atoms, [1,1,1], 3)
   ase_slab.center(vacuum=18.0, axis=2)
   slab_pymatgen = AseAtomsAdaptor().get_structure(ase_slab)
   slab_pymatgen.to(filename='POSCAR-111',fmt='poscar')

   s= Surface(atoms=Si,indices=[1,1,1]).make_surface()
