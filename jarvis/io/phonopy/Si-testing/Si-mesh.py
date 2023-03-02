#!/usr/bin/env python

from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
import yaml
import numpy as np

unitcell = read_vasp("POSCAR")
phonon = Phonopy(
    unitcell,
    [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    primitive_matrix=[[0.875, 0.875, 0.875], [0.125, 0.125, 0.125]],
)

symmetry = phonon.get_symmetry()
print("Space group: %s" % symmetry.get_international_table())

force_sets = parse_FORCE_SETS()
phonon.set_displacement_dataset(force_sets)
phonon.produce_force_constants()
primitive = phonon.get_primitive()

# Born effective charges and dielectric constants are read from BORN file.
nac_params = parse_BORN(primitive, filename="BORN")

phonon.run_mesh([20, 20, 20], with_eigenvectors=True)

mesh_dict = phonon.get_mesh_dict()
