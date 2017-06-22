from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar
from pymatgen.io.vaspio_set import  MPNonSCFVaspInputSet
def get_kpoints(structure=None,nbands=30,encut=600):
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = encut,
           NBANDS=int(nbands)+10,
           ISMEAR = 0,
           EDIFF = '1E-7',


           LCHARG = '.FALSE.',
           NEDOS = 5000,
           ISPIN = 2,
           ISIF = 2,
           IBRION = 1,
           NELM = 400,
           LORBIT = 11,
           NPAR = 4,
           LWAVE = '.FALSE.' )
       incar = Incar.from_dict(incar_dict)
       user_incar_settings={"EDIFF":1E-6,"ISIF":2,"NSW":0,"LORBIT":11,"ENCUT":encut,"LWAVE":'.FALSE.',"PREC":'Accurate'}
       mpvis = MPNonSCFVaspInputSet(user_incar_settings=incar)


       kpoints=mpvis.get_kpoints(structure)
       kpoints.write_file('KP')

strt=Structure.from_file('POSCAR')
get_kpoints(structure=strt)
