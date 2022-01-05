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
from jarvis.io.phonopy.outputs import total_dos
from jarvis.io.phonopy.inputs import PhonopyInputs
import numpy as np
import spglib

try:
    from phonopy import Phonopy
    from phonopy import load
except Exception as exp:
    print("Phonopy is not installed.", exp)
    pass


def get_Phonopy_obj(atoms, phonopy_yaml = None, FC_file = None, factor = None, symprec = 1e-05,
                    scell = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]]),
                    ):
    """
    Create separate method for generating a Phonopy object
    """
    if phonopy_yaml is not None:
        phonon = load(phonopy_yaml, force_constants_filename = FC_file)
    else:
        unitcell = atoms.phonopy_converter()
        prim_mat = np.array(
            PhonopyInputs(atoms).prim_axis().split("=")[1].split(), dtype="float"
        ).reshape(3, 3)
        if factor is None:
            from phonopy.units import VaspToCm
    
            factor = VaspToCm    
        phonon = Phonopy(
            unitcell,
            scell,
            primitive_matrix=prim_mat,
            factor=factor,
            dynamical_matrix_decimals=None,
            force_constants_decimals=None,
            symprec=symprec,
            is_symmetry=True,
            use_lapack_solver=False,
            log_level=1,
        )  
    return phonon

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
        phonopy_obj.run_mesh(mesh)
        self.mesh_dict = phonopy_obj.get_mesh_dict()
        
        
    def select_jdos(self):
        '''

        Returns
        -------
        None.

        '''
        def get_gridpts(self):
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
            latt_vecs = self.phonopy_obj.get_primitive().get_cell()
            positions = self.phonopy_obj.get_primitive().get_scaled_positions()
            atom_type = self.phonopy_obj.get_primitive().get_atomic_numbers()
            cell = (latt_vecs, positions, atom_type)
            mapping, grid = spglib.get_ir_reciprocal_mesh(self.mesh, cell, is_shift=[0, 0, 0])
            return mapping
        
        gridpt_list = get_gridpts(self)
        gridpt_uids = np.unique(gridpt_list)

        """
        Dictionary stores JDOS spectrum for each irreducible q-point
        key : irreducible grid point index
        value : JDOS spectrum [frequency, N1 process, N2 process]
        """
        jdos_ir = np.zeros([len(gridpt_uids), len(self.mesh_dict['frequencies'][0])])
        
        for g, gp in enumerate(gridpt_uids):
            print(gp)
            if self.temperature is not None:
                file = self.directory + 'jdos-m' + "".join(str(m) for m in self.mesh) + '-g' + str(gp) +\
                '-t' + str(self.temperature) + '.dat'
            else:
                file = self.directory + 'jdos-m' + "".join(str(m) for m in self.mesh) + '-g' + str(gp) + '.dat'                
            data = np.genfromtxt(file, dtype = 'float')
            freq = data[:,0]
            N2 = data[:,1] + data[:,2]
            indices = np.digitize(self.mesh_dict['frequencies'][g], freq) #select freqeuncy bin of branch frequency
            print(indices)
            #JDOS values for each irreducible q-point
            jdos_ir[g] = np.array([ (N2[i] + N2[i-1])/2 for i in indices ]) #average the ceiling and floor of the bin
            #Map JDOS back to the full q-point mesh
        #     jdos_full = []
        # for i, qpoint in enumerate(gridpt_list):
        #     jdos_full.append(jdos_ir[qpoint])
        
        return jdos_ir
            
        
    def mode_to_spectral(self, mode_jdos):
        self.phonopy_obj.run_total_dos()
        # Get tetrahedron mesh object
        thm = self.phonopy_obj._total_dos._tetrahedron_mesh
        thm.set(value="I",\
                frequency_points=self.phonopy_obj._total_dos._frequency_points)
        spectral_jdos = np.zeros_like(self.phonopy_obj._total_dos._frequency_points)
        for i, iw in enumerate(thm):
            spectral_jdos += np.sum(iw * mode_jdos[i] * self.phonopy_obj._total_dos._weights[i], axis=1)
        return spectral_jdos
        
            
        
        
        

if __name__ == '__main__':  
    kappa_Si = Kappa('Si-testing/kappa-m111111.hdf5', total_dos_dat =\
                     '../phonopy/Si-testing/total_dos.dat')
    RT_kappa = kappa_Si.kappa(300.)
    from jarvis.core.atoms import Atoms
    test_dir = 'Si-testing/'
    pos = test_dir + "POSCAR-unitcell"
    atoms = Atoms.from_poscar(pos)
    phonon_obj = get_Phonopy_obj(atoms, phonopy_yaml = '../phonopy/Si-testing/phonopy.yaml',\
                                 FC_file = '../phonopy/Si-testing/FORCE_CONSTANTS', scell = np.array([[2, 0, 0],[0, 2, 0],[0, 0, 2]]))
#    phonon_obj = get_Phonopy_obj(atoms, fc, scell = np.array([[2, 0, 0],[0, 2, 0],[0, 0, 2]]))
    jdos_dir = 'Si-testing/jdos_output/'
    jdos = JDOS(phonon_obj, directory = jdos_dir, mesh = [11, 11, 11], temperature = 300)
    jdos_ir = jdos.select_jdos()
    spectral_jdos = jdos.mode_to_spectral(jdos_ir)
    
    