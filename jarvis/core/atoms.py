"""
This module provides classes to specify atomic structure
"""
import numpy as np


class Atoms(object):
    def __init__(self, box=None, coords=None, elements=None, cartesian=False):
        """
     Create atomic structure with lattice, coordinates, atom type and other information

     """
        self.box = box
        self.coords = coords
        self.elements = elements

    def volume(self):
        """
        >>> box = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]
        >>> coords = [[0, 0, 0], [0, 0, 0.5]]
        >>> elements = ["Si", "Si"]
        >>> Si = Atoms(box=box, coords=coords, elements=elements)
        >>> print(Si.volume())
        1000.0
        """
        m = self.box
        vol = float(abs(np.dot(np.cross(m[0], m[1]), m[2])))
        return vol
