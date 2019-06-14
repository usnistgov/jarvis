from pymatgen.io.vaspio_set import MPBSHSEVaspInputSet,MITHSEVaspInputSet,MPNonSCFVaspInputSet
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, VaspInput
def fff():
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
       #mpvis = MPNonSCFVaspInputSet(user_incar_settings=incar)
       kpoints=mpvis.get_kpoints(s)
       kpoints.write_file("KPOINTS")
       print (kpoints)
       import sys
       sys.exit()
       all_kp=kpoints.kpts
       labels=kpoints.labels
       all_labels=''
       for l in labels:
           all_labels=all_labels+str(l)+str(' ')
       for k in all_kp:
           #allline=allline+str(k)
           allline=allline+ str(k[0])+str(' ')+str(k[1])+str(' ')+str(k[2])+str(' ')
       print (allline)
       file=open('bandd.conf','w')
       line=str('FREQUENCY_CONVERSION_FACTOR = 521.471')+'\n'
       file.write(line)
       line=str('ATOM_NAME = Si O')+'\n'
       file.write(line)
       line=str('DIM = 1 1 1')+'\n'
       file.write(line)
       line=str('FORCE_CONSTANTS = READ')+'\n'
       file.write(line)
       line=str("BAND= ")+str(allline)+'\n'
       file.write(line)
       line=str("BAND_LABELS= ")+str(all_labels)+'\n'
       file.write(line)
       file.close()
       # bandplot -o PBAND.png
#fff()
