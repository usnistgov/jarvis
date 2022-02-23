#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 23:06:59 2022

@author: rlg3


Plot Si Gruneisen parameter
"""

import jarvis.io

import h5py
import matplotlib.pyplot as plt

from jarvis.io.phonopy.outputs import get_Phonopy_obj
from jarvis.io.phono3py.inputs import prepare_jdos
from jarvis.io.phono3py.outputs import JDOS
from jarvis.core.atoms import Atoms
import numpy as np


f = h5py.File("gruneisen.hdf5", 'r')

grun_dict = dict(f)


plt.scatter(grun_dict['frequency'], grun_dict['gruneisen'], s = 2)
plt.xlabel('Frequency (THz)')
plt.ylabel('Gruneisen')


'''
Trying to get spectral JDOS...
'''

test_dir = ''
pos = test_dir + "POSCAR-unitcell"
atoms = Atoms.from_poscar(pos)
phonon_obj = get_Phonopy_obj(
    atoms,
    phonopy_yaml=test_dir + "phono3py.yaml",
    FC_file=test_dir + "fc2.hdf5",
    scell=np.array([[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
)

jdos_dir = "../unweighted_jdos/"
jdos = JDOS(phonon_obj, directory=jdos_dir, mesh=[11, 11, 11])
spectral_jdos = jdos.mode_to_spectral(grun_dict['gruneisen'])