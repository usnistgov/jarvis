#! /usr/bin/env python

from nexus import settings,job,run_project,obj
from nexus import generate_physical_system
from nexus import generate_pwscf, generate_pw2qmcpack
from nexus import generate_qmcpack,vmc,loop,linear,dmc
from nexus import read_structure
from numpy import mod, sqrt, array
from qmcpack_input import spindensity
import numpy as np
from structure import optimal_tilematrix
from numpy.linalg import det
settings(
    results = './results',
    pseudo_dir = './pseudopotentials', #location of pseudopotential directory
    sleep   = 1,
    runs    = './runs',
    machine = 'kisir', #machine, defined in the machine.py file in Nexus lib
    )

structure = read_structure('POSCAR', format='poscar') #poscar to be read from JARVIS-ID

boundaries = 'ppp'
supercells = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]] #Matrix used to tile to different sized supercells in kspace (for finite-size convergence)
    
shared_qe = obj(        #QE parameters, if desired can loop over values of ecut and kgrid to find converged values at this stage
    occupations = 'smearing',
    smearing    = 'gaussian',
    degauss     = 0.005,
    input_dft   = 'PBE',
    ecut        = 100,
    conv_thr    = 1.0e-7,
    mixing_beta = 0.2,
    nosym       = True,
    use_folded  = True,
    spin_polarized = True
)

qe_presub = 'module load intel impi openmpi-3.0.1/intel qe-6.4.1' #loaded module for QE, can modify this for NIST cluster 

qe_job = job(nodes=2,threads=40,app='pw.x', presub=qe_presub) #submission details specific to kisir

scales = [1.00] #can modify to scale the lattice 

for scale in scales: 
#Self-consistent DFT calculation to obtain charge density of structure in POSCAR
    temp = structure.copy()
    temp.stretch(scale, scale , scale)
    system = generate_physical_system(
        structure = temp,
        Si = 4,  # This is number of valence electrons for each element 
        kshift   = (0,0,0),
        net_spin = 0
    )

    scf = generate_pwscf(
        identifier = 'scf',                      # log output goes to scf.out
        path       = 'scf-{}'.format(scale),      # directory to run in
        job        = qe_job,
        system     = system,
        input_type = 'scf',
        pseudos    = ['Si.ccECP.upf'], #DFT pseudopotential from the directory
        kgrid      = (4,4,4),
        wf_collect = False,
        **shared_qe
    )
    #Nonself-consistent DFT calculation to read in charge density and tile to appropriate supercell for QMC
    for supercell in supercells:
        scell_vol = det(supercell)
        nscf_kgrid_k = int(np.ceil(2/sqrt(scell_vol))) #The kgrid for the nscf calculation used to tile to different sized supercells, modify this based on system
        nscf_grid = (nscf_kgrid_k, nscf_kgrid_k, nscf_kgrid_k)
        nscf_system = generate_physical_system(
            structure = temp,
            Si = 4,
            tiling = supercell,
            kgrid  = nscf_grid,
            kshift = (0,0,0),
            net_spin = 0,
        )
        nscf = generate_pwscf(
            identifier  = 'nscf',
            path        = 'nscf-{}-{}'.format(scale,scell_vol),
            job         = qe_job ,
            system      = nscf_system,
            input_type  = 'nscf',
            pseudos     = ['Si.ccECP.upf'], #DFT pseudopotential from the directory
            wf_collect   = True,
            dependencies = (scf,'charge_density'),
            **shared_qe
            )
#Wavefunction conversion (p2qmcpack.x run on top of pw.x in same folder as NSCF calculation)
        p2q = generate_pw2qmcpack(
            identifier   = 'p2q',
            path         = 'nscf-{}-{}'.format(scale,scell_vol),
            job          = job(cores= 1,threads= 1,app='pw2qmcpack.x',presub=qe_presub),
            write_psir   = False,
            dependencies = (nscf,'orbitals'),
            )
run_project()