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
run 0
#fix 1a all qeq/comb 1 0.0001 file fq.out
write_data data0

