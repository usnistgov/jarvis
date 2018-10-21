#!/bin/bash
#PBS -N bulk@cellmax4
#PBS -l walltime=0:30:00
#PBS -o test.log
#PBS -m abe
#PBS -j oe
#PBS -r n
#PBS -l nodes=1:ppn=1
cd /users/knc6/Software/jarvis/jarvis/lammps/examples/Al03.eam.alloy_nist/bulk@mp-134_fold/bulk@cellmax4
mpirun /cluster/bin/lmp_ctcms-14439-knc6-2 <in.elastic >out
