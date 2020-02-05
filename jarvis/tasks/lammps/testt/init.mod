variable dump_file string "None"
variable up  equal 1.0e-6
variable cfac  equal 1.0e-4
variable cunits  string GPa
variable etol  equal 0
variable ftol  equal 1.0e-10
variable maxiter equal 1000
variable maxeval equal 10000
variable dmax equal 1.0e-2
variable data_file string "data"
units metal 
atom_style charge 
boundary p p p 
atom_modify sort 0 0.0 

read_data data

### interactions 
mass 1 26.9815386
