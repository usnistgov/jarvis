import numpy as np


class Lattice(object):
    def __init__(self, box=None):
        """
        >>> box=[[10,0,0],[0,10,0],[0,0,10]]
        >>> lat=Lattice(box)
        >>> lat.lat_lengths()
        [10.0, 10.0, 10.0]
        >>> Lattice(box).lattice()[0][0]
        10.0
        >>> lat.lat_angles()
        [90.0, 90.0, 90.0]
        >>> round(lat.inv_lattice()[0][0],2)
        0.1
        >>> [round(i,2) for i in lat.lat_angles(radians=True)]
        [1.57, 1.57, 1.57]
        >>> frac_coords=[[0,0,0],[0.5,0.5,0.5]]
        >>> lat.cart_coords(frac_coords)[1][1]
        5.0
        >>> cart_coords=[[0,0,0],[5,5,5]]
        >>> lat.frac_coords(cart_coords)[1][1]
        0.5
        >>> box=[[1e-8,0,0],[0,1e-8,0],[0,0,1e-8]]
        >>> lat=Lattice(box)
        >>> [round(i,2) for i in lat.lat_angles()]
        [90.0, 90.0, 90.0]
        """
        self._lat = np.array(box, dtype=np.float64).reshape((3, 3))
        self._inv_lat = None

    def lat_lengths(self):
        return [np.linalg.norm(v) for v in self._lat]

    def lat_angles(self, tol=1e-2, radians=False):
        lengths = self.lat_lengths()
        angles = []
        for i in range(3):
            j = i - 1
            k = i - 2
            leng2 = lengths[j] * lengths[k]
            if leng2 > tol:
                tmp = np.dot(self._lat[j], self._lat[k]) / leng2
                angle = 180.0 * np.arccos(tmp) / np.pi
            else:
                angle = 90.0
            angles.append(angle)
        if radians:
            angles = [angle * np.pi / 180.0 for angle in angles]
        return angles

    def lattice(self):
        return self._lat

    def inv_lattice(self):
        # if self._inv_lat == None:
        self._inv_lat = np.linalg.inv(self._lat)
        return self._inv_lat

    def cart_coords(self, frac_coords):
        return np.dot(np.array(frac_coords), self._lat)

    def frac_coords(self, cart_coords):
        return np.dot(np.array(cart_coords), self.inv_lattice())


"""
if __name__=='__main__':

        box=[[10,0,0],[0,10,0],[0,0,10]]
        lat=Lattice(box)
        frac_coords=[[0,0,0],[0.5,0.5,0.5]]
        print (lat.cart_coords(frac_coords)[1][1])
        
        cart_coords=[[0,0,0],[5,5,5]]
        print (lat.frac_coords(cart_coords)[1][1])
        
"""
