pair_style eam/alloy 
pair_coeff * * /users/knc6/Software/jarvis/jarvis/lammps/examples/Al03.eam.alloy  Al 
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes
thermo          1
thermo_style custom step temp press cpu pxx pyy pzz pxy pxz pyz ke pe etotal vol lx ly lz atoms
thermo_modify norm no
