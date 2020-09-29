pair_style eam/alloy 
pair_coeff * * /users/knc6/Software/LAMMPS/lammps-master/potentials/Al_zhou.eam.alloy  Al 
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes
min_style  cg
min_modify           dmax ${dmax} line quadratic
thermo          1
thermo_style custom step temp press cpu pxx pyy pzz             pxy pxz pyz ke pe etotal vol lx ly lz atoms
thermo_modify norm no
