# Compute elastic constant tensor for a crystal
#
# Written by Aidan Thompson (Sandia, athomps@sandia.gov)
#
#  This script uses the following three include files.
#
#   init.mod      (must be modified for different crystal structures)
# 	       	  Define units, deformation parameters and initial
#		  configuration of the atoms and simulation cell.  
#
#
#   potential.mod    (must be modified for different pair styles)
# 		     Define pair style and other attributes 
#		     not stored in restart file
#
#
#   displace.mod    (displace.mod should not need to be modified)
# 		    Perform positive and negative box displacements 
# 		    in direction ${dir} and size ${up}. 
# 		    It uses the resultant changes 
#		    in stress to compute one
# 		    row of the elastic stiffness tensor
#		    
#		    Inputs variables:
#		    	   dir = the Voigt deformation component 
#		    		    (1,2,3,4,5,6)  
#		    Global constants:
#       	    	   up = the deformation magnitude (strain units)
#       		   cfac = conversion from LAMMPS pressure units to 
#               	   output units for elastic constants 
#
#
#  To run this on a different system, it should only be necessary to 
#  modify the files init.mod and potential.mod. In order to calculate
#  the elastic constants correctly, care must be taken to specify
#  the correct units in init.mod (units, cfac and cunits). It is also
#  important to verify that the minimization of energy w.r.t atom
#  positions in the deformed cell is fully converged.
#  One indication of this is that the elastic constants are insensitive
#  to the choice of the variable ${up} in init.mod. Another is to check
#  the final max and two-norm forces reported in the log file. If you know
#  that minimization is not required, you can set maxiter = 0.0 in 
#  init.mod. 
#
#  There are two alternate versions of displace.mod provided.
#  They are displace_restart.mod and displace_reverse.mod. 
#  The former resets the box using a restart file while 
#  the latter reverses the deformation. Copy whichever
#  one you like best to displace.mod.
#

include init.mod
include potential.mod
#fix 3 all nph aniso 0.0 0.0 1.0
#fix langevin all langevin 0.0 0.0 1.0 12233 

#fix 3 all box/relax tri 0.0 vmax 0.001
#fix            3 all box/relax aniso 0.0 vmax 0.001
#fix 3 all box/relax  aniso 0.0
#dump            2a all cfg 100 *.cfg mass type xs ys zs  vx vy vz fx fy fz
#dump_modify     2a element  Mo S

#fix 1a all qeq/comb 1 0.0001 file fq.out
#fix             2 all qeq/reax 1 0.0 10.0 1e-6 /scratch/lfs/kamal/JARVIS/All2/ReaxFF/param.qeq
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
#minimize ${etol} ${ftol} ${maxiter} ${maxeval}
#minimize ${etol} ${ftol} ${maxiter} ${maxeval}

variable tmp equal pxx
variable pxx0 equal ${tmp}
variable tmp equal pyy
variable pyy0 equal ${tmp}
variable tmp equal pzz
variable pzz0 equal ${tmp}
variable tmp equal pyz
variable pyz0 equal ${tmp}
variable tmp equal pxz
variable pxz0 equal ${tmp}
variable tmp equal pxy
variable pxy0 equal ${tmp}

variable tmp equal lx
variable lx0 equal ${tmp}
variable tmp equal ly
variable ly0 equal ${tmp}
variable tmp equal lz
variable lz0 equal ${tmp}

# These formulas define the derivatives w.r.t. strain components
# Constants uses $, variables use v_ 
variable d1 equal -(v_pxx1-${pxx0})/(v_delta/v_len0)*${cfac}
variable d2 equal -(v_pyy1-${pyy0})/(v_delta/v_len0)*${cfac}
variable d3 equal -(v_pzz1-${pzz0})/(v_delta/v_len0)*${cfac}
variable d4 equal -(v_pyz1-${pyz0})/(v_delta/v_len0)*${cfac}
variable d5 equal -(v_pxz1-${pxz0})/(v_delta/v_len0)*${cfac}
variable d6 equal -(v_pxy1-${pxy0})/(v_delta/v_len0)*${cfac}

# Write restart
#unfix 3
minimize ${etol} ${ftol} ${maxiter} ${maxeval}
#fix 1a all qeq/comb 1 0.0001 file fq.out
write_restart restart.equil
write_data data0

# uxx Perturbation

variable dir equal 1
include /users/knc6/displace.mod

# uyy Perturbation

variable dir equal 2
include /users/knc6/displace.mod

# uzz Perturbation

variable dir equal 3
include /users/knc6/displace.mod

# uyz Perturbation

variable dir equal 4
include /users/knc6/displace.mod

# uxz Perturbation

variable dir equal 5
include /users/knc6/displace.mod

# uxy Perturbation

variable dir equal 6
include /users/knc6/displace.mod

# Output final values

variable C11all equal ${C11}
variable C22all equal ${C22}
variable C33all equal ${C33}

variable C12all equal 0.5*(${C12}+${C21})
variable C13all equal 0.5*(${C13}+${C31})
variable C23all equal 0.5*(${C23}+${C32})

variable C44all equal ${C44}
variable C55all equal ${C55}
variable C66all equal ${C66}

variable C14all equal 0.5*(${C14}+${C41})
variable C15all equal 0.5*(${C15}+${C51})
variable C16all equal 0.5*(${C16}+${C61})

variable C24all equal 0.5*(${C24}+${C42})
variable C25all equal 0.5*(${C25}+${C52})
variable C26all equal 0.5*(${C26}+${C62})

variable C34all equal 0.5*(${C34}+${C43})
variable C35all equal 0.5*(${C35}+${C53})
variable C36all equal 0.5*(${C36}+${C63})

variable C45all equal 0.5*(${C45}+${C54})
variable C46all equal 0.5*(${C46}+${C64})
variable C56all equal 0.5*(${C56}+${C65})

# For Stillinger-Weber silicon, the analytical results
# are known to be (E. R. Cowley, 1988):
#               C11 = 151.4 GPa
#               C12 = 76.4 GPa
#               C44 = 56.4 GPa

print "Elastic Constant C11all = ${C11all} ${cunits}"
print "Elastic Constant C22all = ${C22all} ${cunits}"
print "Elastic Constant C33all = ${C33all} ${cunits}"

print "Elastic Constant C12all = ${C12all} ${cunits}"
print "Elastic Constant C13all = ${C13all} ${cunits}"
print "Elastic Constant C23all = ${C23all} ${cunits}"

print "Elastic Constant C44all = ${C44all} ${cunits}"
print "Elastic Constant C55all = ${C55all} ${cunits}"
print "Elastic Constant C66all = ${C66all} ${cunits}"

print "Elastic Constant C14all = ${C14all} ${cunits}"
print "Elastic Constant C15all = ${C15all} ${cunits}"
print "Elastic Constant C16all = ${C16all} ${cunits}"

print "Elastic Constant C24all = ${C24all} ${cunits}"
print "Elastic Constant C25all = ${C25all} ${cunits}"
print "Elastic Constant C26all = ${C26all} ${cunits}"

print "Elastic Constant C34all = ${C34all} ${cunits}"
print "Elastic Constant C35all = ${C35all} ${cunits}"
print "Elastic Constant C36all = ${C36all} ${cunits}"

print "Elastic Constant C45all = ${C45all} ${cunits}"
print "Elastic Constant C46all = ${C46all} ${cunits}"
print "Elastic Constant C56all = ${C56all} ${cunits}"

