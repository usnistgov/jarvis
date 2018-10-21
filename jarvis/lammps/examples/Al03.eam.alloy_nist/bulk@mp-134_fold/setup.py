from jarvis.lammps.NEW_LAMMPS import main_func
from pymatgen.io.vasp.inputs import  Poscar
p=Poscar.from_file("/users/knc6/Software/jarvis/jarvis/lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/POSCAR")
main_func(mat=p,parameters={'pair_coeff': '/users/knc6/Software/jarvis/jarvis/lammps/examples/Al03.eam.alloy', 'control_file': '/users/knc6/inelast.mod', 'atom_style': 'charge', 'pair_style': 'eam/alloy', 'exec': 'mpirun /cluster/bin/lmp_ctcms-14439-knc6-2 <in.elastic >out'})
