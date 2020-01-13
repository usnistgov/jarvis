import numpy as np
from jarvis.core.atoms import Atoms

class Poscar(object):
    def __init__(self,atoms,comment='System'):
          self.atoms = atoms
    @staticmethod
    def from_file(filename:str='POSCAR'):
    
        with open(filename, "r") as f:
          return Poscar.from_string(f.read())

    @staticmethod
    def from_string(lines):
       text = lines.splitlines()
       comment=text[0]
       scale = float(text[1])       
       lattice_mat = []
       lattice_mat.append([float(i) for i in text[2].split()])
       lattice_mat.append([float(i) for i in text[3].split()])
       lattice_mat.append([float(i) for i in text[4].split()])
       lattice_mat = np.array(lattice_mat)

       uniq_elements = text[5].split()
       element_count = np.array([int(i) for i in text[6].split()])
       elements = []
       for i,ii in enumerate(element_count):
          for j in range(ii):
            elements.append(uniq_elements[i])
       cartesian = True
       if 'd' in text[7] or 'D' in text[7]:
           cartesian = False
       #print ('cartesian poscar=',cartesian,text[7])
       num_atoms = int(np.sum(element_count))
       coords = []
       for i in range(num_atoms):
          coords.append( [float(i) for i in text[8+i].split()[0:3]])
       coords = np.array(coords)
       atoms = Atoms(lattice_mat=lattice_mat, coords=coords, elements=elements,cartesian=cartesian)
       #print (atoms)
       formula = atoms.composition.formula
       return Poscar(atoms,comment=formula)
    def __repr__(self):
        return str(self.atoms)
if __name__=='__main__':
    p = Poscar.from_file(filename='/rk2/knc6/JARVIS-DFT/2D-1L/POSCAR-mp-2815-1L.vasp_PBEBO/MAIN-RELAX-Surf-mp-2815/POSCAR')
    print (p.atoms*[3,3,3])
