#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 21:34:38 2021

@author: rlg3

Module for post-processing phono3py output
kappa-mxxx.hdf, jdos-mxxx-gy-tz.hdf, phono3py_disp.yaml


Notes for inputs.py:
    tags : generate automatically from the inputs.py. Include info about
    isotope scattering, boundary scattering, etc.

Think I should initialize a Phonopy object to run these.
    
Might need a mesh_dict from phonon.get_mesh_dict() so that I have the phonon
bandstructure properties
"""

import h5py
from jarvis.core.kpoints import Kpoints3D
from jarvis.core.composition import Composition
from jarvis.io.phonopy.outputs import total_dos, get_Phonopy_obj
import numpy as np
import spglib

try:
    from phonopy import Phonopy
except Exception as exp:
    print("Phonopy is not installed.", exp)
    pass

class Kappa():
    
    def __init__(self,
            filename = "",
            mesh = None,
            total_dos_dat = "",
            temperatures = None,
            kappa_format = 'scalar_xx', #need this?
            composition = None,
            tags = None):
        '''

        Parameters
        ----------
        filename : string,required
            Name of kappa-mxxxx.hdf file. The default is "".
        temperatures : list, optional
            List of temperatures to use for post-processing. The default is None.
            If None, uses all temperatures computed in thermal conductivity run.
        kappa_format : string, optional
            Desired kappa output format. The default is 'scalar_xx'.
            Other choices : "tensor" or ...
        composition : Composition object of jarvis.core.composition, optional
            Composition object for material used in property calculation. The default is None.
        tags : string, optional
        Descriptors of thermal conductivity calculation. The default is None.
        '''
        f = h5py.File(filename, 'r')
        self.dict = dict(f)
        if temperatures:
            self.temperatures = temperatures
        else:
            self.temperatures = list(self.dict['temperature'])
        self.composition = composition
        self.tags = tags
        self.n_qpoint = np.shape(f['qpoint'])[0]
        self.total_dos = np.array(total_dos(total_dos_dat))
        self.kappa_format = kappa_format
        
    def to_dict(self):
        return self.dict
    
    
    def kappa(self, T):
        if T not in self.temperatures:
            raise Exception('Thermal conductivity not evaluated at this temperature.')
        T_indx = self.temperatures.index(T)
        if self.kappa_format == 'tensor':
            return self.dict['kappa'][T_indx]
        if self.kappa_format == 'scalar_xx':
            return self.dict['kappa'][T_indx][0]
        
    
    
    # def generate_spectral_property(self,
    #                                property_key = 'kappa',
    #                                T = 300):
    #     '''

    #     Parameters
    #     ----------
    #     property_key : TYPE, optional
    #         DESCRIPTION. The default is 'kappa'.
    #     T : TYPE, optional
    #         DESCRIPTION. The default is 300.

    #     Returns
    #     -------
    #     None.

    #     '''
        

class JDOS():
    
    def __init__(self,
                 phonopy_obj,
                  directory,
                  mesh = [1, 1, 1],
                  temperature = None):
        '''
        

        Parameters
        ----------
        directory : string, required
            DESCRIPTION. The default is "".
        qpoints : TYPE, optional
            DESCRIPTION. The default is None.
        temperature : integer, optional
        Temperature for a weighted JDOS calculation. The default is None. When
        None, unweighted JDOS is computed.

        '''
        self.mesh = mesh
        self.temperature = temperature
        self.directory = directory
        self.phonopy_obj = phonopy_obj
        
    def select_jdos(self):
        '''

        Returns
        -------
        None.

        '''
        def get_gridpt_indices(self):
            '''
            Generates list of gridpoint indices for JDOS calculation

            Parameters
            ----------
            phonopy_obj : TYPE
                DESCRIPTION.

            Returns
            -------
            None.

            '''
            latt_vecs = self.phonopy_obj.get_primitive.get_cell()
            positions = self.phonopy_obj.get_primitive().get_scaled_positions()
            atom_type = self.phonopy_obj.get_primitive().get_atomic_numbers()
            cell = (latt_vecs, positions, atom_type)
            mapping, grid = spglib.get_ir_reciprocal_mesh(self.mesh, cell, is_shift=[0, 0, 0])
            gridpt_list = np.unique(mapping)
            return gridpt_list
        gridpt_list = get_gridpt_indices(self)
        return gridpt_list
            
        
        
        
        
        
        
        

if __name__ == '__main__':  
    kappa_Si = Kappa('Si-testing/kappa-m111111.hdf5', total_dos_dat =\
                     '../phonopy/Si-testing/total_dos.dat')
    RT_kappa = kappa_Si.kappa(300.)
    from jarvis.core.atoms import Atoms
    test_dir = 'Si-testing/'
    pos = test_dir + "POSCAR-unitcell"
    atoms = Atoms.from_poscar(pos)
    phonon_obj = get_Phonopy_obj(atoms, scell = np.array([[2, 0, 0],[0, 2, 0],[0, 0, 2]]))
    jdos_dir = 'Si-testing/jdos_output'
    jdos = JDOS(phonon_obj, directory = jdos_dir, mesh = [11, 11, 11], temperature = 300)
    
    