#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 21:34:38 2021

@author: rlg3

Module for post-processing phono3py output
kappa-mxxx.hdf, jdos-mxxx-gy-tz.hdf


Notes for inputs.py:
    tags : generate automatically from the inputs.py. Include info about
    isotope scattering, boundary scattering, etc.
"""

import h5py
from jarvis.core.kpoints import Kpoints3D
from jarvis.core.composition import Composition


class Kappa():
    
    def __init__(self,
            filename = "",
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
            Composition object for material used in proeprty calculation. The default is None.
        tags : string, optional
        Descriptors of thermal conductivity calculation. The default is None.

        Returns
        -------
        Kappa object generated from output hdf file

        '''
        f = h5py.File(filename, 'r')
        self.dict = dict(f)
        self.temperatures = temperatures
        self.composition = composition
        self.tags = tags
        
    def to_dict(self):
        return self.dict
            
        
        

# class JDOS():
    
#     def __init__():
        

if __name__ == '__main__':  
    kappa_Si = Kappa('Si-testing/kappa-m111111.hdf5')
    