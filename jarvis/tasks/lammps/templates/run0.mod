
include init.mod
include potential.mod
#fix 1a all qeq/comb 1 0.0001 file fq.out
dump 1 all custom 1 *.dump id  fx fy fz
write_data data0
run 0
