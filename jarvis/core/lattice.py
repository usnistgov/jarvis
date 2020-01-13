import numpy as np


class Lattice(object):
    def __init__(self, lattice_mat=None):
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
        self._lat = np.array(lattice_mat, dtype=np.float64).reshape((3, 3))
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

    def reciprocal_lattice(self):
        return Lattice(2*np.pi*np.linalg.inv(self._lat).T)

    def get_points_in_sphere(self,frac_points,center,r):
        recp_len = np.array(self.reciprocal_lattice().lat_lengths()) / (2 * np.pi)
        nmax = float(r) * recp_len + 0.01

        # Get the fractional coordinates of the center of the sphere
        pcoords = self.frac_coords(center)
        center = np.array(center)

        # Prepare the list of output atoms
        n = len(frac_points)
        fcoords = np.array(frac_points) % 1
        indices = np.arange(n)

        # Generate all possible images that could be within `r` of `center`
        mins = np.floor(pcoords - nmax)
        maxes = np.ceil(pcoords + nmax)
        arange = np.arange(start=mins[0], stop=maxes[0], dtype=np.int)
        brange = np.arange(start=mins[1], stop=maxes[1], dtype=np.int)
        crange = np.arange(start=mins[2], stop=maxes[2], dtype=np.int)
        arange = arange[:, None] * np.array([1, 0, 0], dtype=np.int)[None, :]
        brange = brange[:, None] * np.array([0, 1, 0], dtype=np.int)[None, :]
        crange = crange[:, None] * np.array([0, 0, 1], dtype=np.int)[None, :]
        images = arange[:, None, None] + brange[None, :, None] + crange[None, None, :]

        # Generate the coordinates of all atoms within these images
        shifted_coords = fcoords[:, None, None, None, :] + images[None, :, :, :, :]

        # Determine distance from `center`
        cart_coords = self.cart_coords(fcoords)
        cart_images = self.cart_coords(images)
        coords = cart_coords[:, None, None, None, :] + cart_images[None, :, :, :, :]
        coords -= center[None, None, None, None, :]
        coords **= 2
        d_2 = np.sum(coords, axis=4)

        # Determine which points are within `r` of `center`
        within_r = np.where(d_2 <= r ** 2)        
        return  (
                shifted_coords[within_r],
                np.sqrt(d_2[within_r]),
                indices[within_r[0]],
                images[within_r[1:]],
            )

 
    def find_all_matches(self,other_lattice,ltol=1e-5,atol=1):
        """
        Adapted from pymatgen, which is available under MIT license:
        """
        lengths = other_lattice.lat_lengths()
        angles = other_lattice.lat_angles()
        print ('angles',angles)
        alpha = angles[0]
        beta = angles[1]
        gamma = angles[2]
        frac, dist, _, _ = self.get_points_in_sphere(
            [[0, 0, 0]], [0, 0, 0], max(lengths) * (1 + ltol)
        )
        cart = self.cart_coords(frac)
        # this can't be broadcast because they're different lengths
        inds = [
            np.logical_and(dist / l < 1 + ltol, dist / l > 1 / (1 + ltol))  # type: ignore
            for l in lengths
        ]
        c_a, c_b, c_c = (cart[i] for i in inds)
        f_a, f_b, f_c = (frac[i] for i in inds)
        l_a, l_b, l_c = (np.sum(c ** 2, axis=-1) ** 0.5 for c in (c_a, c_b, c_c))

        def get_angles(v1, v2, l1, l2):
            x = np.inner(v1, v2) / l1[:, None] / l2
            x[x > 1] = 1
            x[x < -1] = -1
            angles = np.arccos(x) * 180.0 / np.pi
            return angles

        alphab = np.abs(get_angles(c_b, c_c, l_b, l_c) - alpha) < atol
        betab = np.abs(get_angles(c_a, c_c, l_a, l_c) - beta) < atol
        gammab = np.abs(get_angles(c_a, c_b, l_a, l_b) - gamma) < atol

        for i, all_j in enumerate(gammab):
            inds = np.logical_and(
                all_j[:, None], np.logical_and(alphab, betab[i][None, :])
            )
            for j, k in np.argwhere(inds):
                scale_m = np.array((f_a[i], f_b[j], f_c[k]), dtype=np.int)  # type: ignore
                if abs(np.linalg.det(scale_m)) < 1e-8:
                    continue

                aligned_m = np.array((c_a[i], c_b[j], c_c[k]))

                rotation_m = np.linalg.solve(aligned_m, other_lattice._lat)
                yield Lattice(aligned_m), rotation_m, scale_m



    def find_matches(self,other_lattice,ltol=1e-5,atol=1):
           for x in self.find_all_matches(
                other_lattice, ltol, atol):
                return x
           return None

# if __name__=='__main__':
#
#        box=[[10,0,0],[0,10,0],[0,0,10]]
#        lat=Lattice(box)
#        frac_coords=[[0,0,0],[0.5,0.5,0.5]]
#        print (lat.cart_coords(frac_coords)[1][1])
#
#        cart_coords=[[0,0,0],[5,5,5]]
#        print (lat.frac_coords(cart_coords)[1][1])
