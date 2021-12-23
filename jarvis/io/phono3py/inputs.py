import numpy as np
import h5py


class Phono3pyInputs(object):
	'''
	Generate inputs for post-processing.

	Obtain fc2.hdf, fc3.hdf, phono3py_disp.yaml prior to execution
	'''
	def __init__(self, atoms):
