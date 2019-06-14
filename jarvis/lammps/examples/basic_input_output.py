from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from jarvis.lammps.jlammps import main_func, write_lammps_data, write_lammps_in
from pymatgen.io.vasp.inputs import Poscar
import os
"""
Step-1: writing LAMMPS data file
Import VASP's POSCAR or cif format file,make it 15x15x15 with conventional cell
"""
pos = os.path.join(
    os.path.dirname(__file__), "..", "..", "vasp", "examples", "SiOptb88", "POSCAR"
)

c_size = 15
p = Structure.from_file(pos)
sg_mat = SpacegroupAnalyzer(p)
mat_cvn = sg_mat.get_conventional_standard_structure()
# mat_cvn.to(fmt='poscar',filename='POSCAR_cvb.vasp') #can be visualized using VESTA
dim1 = int((float(c_size) / float(max(abs(mat_cvn.lattice.matrix[0]))))) + 1
dim2 = int(float(c_size) / float(max(abs(mat_cvn.lattice.matrix[1])))) + 1
dim3 = int(float(c_size) / float(max(abs(mat_cvn.lattice.matrix[2])))) + 1
p.make_supercell([dim1, dim2, dim3])
write_lammps_data(p, file="dat.dat")  # can be visulaized using Ovito

"""
Step-2: writing LAMMPS input file
Adjust parameters according to LAMMPS executable, force-field file, main input file that you'll be using, pair style and atom style
"""
parameters = {
    "exec": "mpirun /cluster/bin/lmp_ctcms-14439-knc6 <in.elastic >out",
    "pair_coeff": "/users/knc6/Software/jarvis/jarvis/lammps/examples/Mishin-Ni-Al-2009.eam.alloy",
    "control_file": "/users/knc6/inelast.mod",
    "pair_style": "eam/alloy",
    "atom_style": "charge",
    "cluster": "pbs",
}
write_lammps_in(
    structure=p,
    lammps_in="init.mod",
    lammps_in1="potential.mod",
    lammps_in2="in.elastic",
    parameters=parameters,
)

"""
Step-3
If you survive the above two steps, lets submit job to computer-cluster


pos = Poscar(p)
pos.comment = "bulk@Al"
main_func(mat=pos, parameters=parameters)
"""
