try:
 from pymatgen.io.vaspio_set import MPBSHSEVaspInputSet,MITHSEVaspInputSet,MPNonSCFVaspInputSet
except:
   pass
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput
def hse():
       allline=''
       s=Structure.from_file("POSCAR")
       incar_dict = dict(
           PREC = 'Accurate',
           ENCUT = 45678,
           NBANDS=456789,
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
       user_incar_settings={"EDIFF":1E-6,"ISIF":2,"NSW":0,"LORBIT":11,"ENCUT":34,"LWAVE":'.FALSE.',"PREC":'Accurate'}
       mpvis = MPBSHSEVaspInputSet(user_incar_settings=incar)
       kpoints=mpvis.get_kpoints(s)
       kpoints.write_file("KPOINTS")
#hse()
