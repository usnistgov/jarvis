"""This module provides classes to specify atomic structure."""
# import line_profiler
import numpy as np
import matplotlib.pyplot as plt
import math
from toolz.curried import pipe


def check_array(arr, type_=None, shape=None):
    """
    Check type and shape of np.ndarray.

    Follwing arguments are needed.
    Args:
      x: array to check

      type_: the type of the array (int, float), None to not check

      shape: the shape of the array, None to not check, -1 to ignore a
        dimension

    """
    if not isinstance(arr, np.ndarray):
        raise TypeError("not an array", arr)


def special_arange(dims):
    """
    Multiple dimensional arange.

    This could be extended to any dimension.
    Args:

      dims: sequence of ints of length 3

    Returns:
      an array of shape (np.prod(dims), 3) with each index having a
      different arrangement of arange.
    This function implements the equivalent of a multidimensional for
    loop, which feels like multidimensional arange (see the test)

    >>> dim0, dim1, dim2 = 2, 4, 3
    >>> arr_test = np.zeros((dim0 * dim1 * dim2, 3), dtype=int)
    >>> counter = 0
    >>> for i in np.arange(dim0):
    ...     for j in np.arange(dim1):
    ...         for k in np.arange(dim2):
    ...             arr_test[counter, 0] = i
    ...             arr_test[counter, 1] = j
    ...             arr_test[counter, 2] = k
    ...             counter += 1

    >>> arr_actual = special_arange((dim0, dim1, dim2))
    >>> print(arr_actual)
    [[0 0 0]
     [0 0 1]
     [0 0 2]
     [0 1 0]
     [0 1 1]
     [0 1 2]
     [0 2 0]
     [0 2 1]
     [0 2 2]
     [0 3 0]
     [0 3 1]
     [0 3 2]
     [1 0 0]
     [1 0 1]
     [1 0 2]
     [1 1 0]
     [1 1 1]
     [1 1 2]
     [1 2 0]
     [1 2 1]
     [1 2 2]
     [1 3 0]
     [1 3 1]
     [1 3 2]]
    >>> assert np.all(arr_actual == arr_test)

    """
    # This can be made much more functional and for arbitrary
    # dimensions if necessary
    dim0, dim1, dim2 = dims
    i = np.repeat(np.repeat(np.arange(dim0), dim2), dim1)
    j = np.tile(np.repeat(np.arange(dim1), dim2), dim0)
    k = np.tile(np.tile(np.arange(dim2), dim1), dim0)
    return np.concatenate((i[:, None], j[:, None], k[:, None]), axis=-1)


def calc_structure_data(coords, box, all_symbs, c_size):
    """
    Calcuate a dictionary of structure data.

    Args:

      coords: the coordinates for each element

      box: the lattic matrix

      all_symbs: the elements

      c_size: the c size

    Returns:
      a set of structure data

    >>> coords = np.array([[0, 0, 0], [0.25, 0.2, 0.25]])
    >>> lat = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    >>> box = np.array(lat)
    >>> elements = ["Si", "Si"]
    >>> c_size = 10.0
    >>> data = calc_structure_data(coords, box, elements, c_size)
    >>> assert len(data['coords']) == 128
    >>> assert np.allclose(data['coords'][9], [0.    , 0.5   , 0.25  ])
    >>> assert np.all(data['dim'] == [4, 4, 4])
    >>> assert len(data['new_symbs']) == 128

    """
    check_array(coords, float, (len(all_symbs), 3))
    check_array(box, float, (3, 3))

    def make_coords(dim):
        return pipe(
            dim,
            special_arange,
            lambda x: (coords[:, None, :] + x[None, :, :])
            / dim[None, None, :],
            lambda x: np.reshape(x, (-1, 3)),
            lambda x: dict(
                coords=x,
                dim=dim,
                nat=len(x),
                lat=dim[:, None] * box,
                new_symbs=np.repeat(all_symbs, np.prod(dim)),
            ),
        )

    return make_coords((c_size / np.max(np.abs(box), axis=1)).astype(int) + 1)


class NeighborsAnalysis(object):
    """Get neighbor informations (RDF,ADF,DDF) for Atoms object."""

    def __init__(
        self,
        atoms=None,
        max_n=500,
        rcut1=None,
        max_cut=10.0,
        rcut2=None,
        verbose=False,
    ):
        """
        Initialize the function.

        >>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
        >>> coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
        >>> elements = ["Si", "Si"]
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
        >>> distributions = NeighborsAnalysis(Si).get_all_distributions
        >>> distributions['rdf']
        1
        """
        self._atoms = atoms
        self.max_n = max_n
        self.max_cut = max_cut
        self.nb_warn = ""
        if rcut1 is None or rcut2 is None:
            rcut1, rcut2 = self.get_dist_cutoffs()
        self.rcut1 = rcut1
        self.rcut2 = rcut2
        self.verbose = verbose
        if self.nb_warn != "" and self.verbose:
            print(self.nb_warn)
            print(
                "Try setting higher max_n in the NeighborsAnalysis module"
                + atoms.get_string()
            )

    def get_structure_data(self, c_size=10.0):
        """Provide non-repetitive structure information."""
        return calc_structure_data(
            self._atoms.frac_coords,
            self._atoms.lattice_mat,
            self._atoms.elements,
            c_size,
        )

    # @profile
    # kernprof -v -l neighbors.py
    def nbor_list(self, rcut=10.0, c_size=12.0):
        """Generate neighbor info."""
        max_n = self.max_n
        nbor_info = {}
        struct_info = self.get_structure_data(c_size)
        coords = np.array(struct_info["coords"])
        nat = struct_info["nat"]
        new_symbs = struct_info["new_symbs"]
        lat = np.array(struct_info["lat"])
        different_bond = {}
        nn = np.zeros((nat), dtype="int")
        # print ('max_n, nat',max_n, nat)
        dist = np.zeros((max_n, nat))
        nn_id = np.zeros((max_n, nat), dtype="int")
        bondx = np.zeros((max_n, nat))
        bondy = np.zeros((max_n, nat))
        bondz = np.zeros((max_n, nat))
        for i in range(nat):
            for j in range(i + 1, nat):
                diff = coords[i] - coords[j]
                ind = np.where(np.fabs(diff) > np.array([0.5, 0.5, 0.5]))
                diff[ind] -= np.sign(diff[ind])
                new_diff = np.dot(diff, lat)
                dd = np.linalg.norm(new_diff)
                if dd < rcut and dd >= 0.1:
                    sumb_i = new_symbs[i]
                    sumb_j = new_symbs[j]
                    comb = "_".join(
                        sorted(str(sumb_i + "_" + sumb_j).split("_"))
                    )
                    different_bond.setdefault(comb, []).append(dd)

                    # print ('dd',dd)
                    nn_index = nn[i]  # index of the neighbor
                    nn_index1 = nn[j]  # index of the neighbor
                    if nn_index < max_n and nn_index1 < max_n:
                        nn[i] = nn[i] + 1
                        dist[nn_index][i] = dd  # nn_index counter id
                        nn_id[nn_index][i] = j  # exact id
                        bondx[nn_index][i] = new_diff[0]
                        bondy[nn_index][i] = new_diff[1]
                        bondz[nn_index][i] = new_diff[2]
                        nn[j] = nn[j] + 1
                        dist[nn_index1][j] = dd  # nn_index counter id
                        nn_id[nn_index1][j] = i  # exact id
                        bondx[nn_index1][j] = -new_diff[0]
                        bondy[nn_index1][j] = -new_diff[1]
                        bondz[nn_index1][j] = -new_diff[2]
                    else:
                        self.nb_warn = (
                            "Very large nearest neighbors observed "
                            + str(nn_index)
                        )
        nbor_info["dist"] = dist
        nbor_info["nat"] = nat
        nbor_info["nn_id"] = nn_id
        nbor_info["nn"] = nn
        nbor_info["bondx"] = bondx
        nbor_info["bondy"] = bondy
        nbor_info["bondz"] = bondz
        # print ('nat',nat)

        return nbor_info

    def get_rdf(self, plot=False):
        """Calculate radial distribution function."""
        nbor_info = self.nbor_list(c_size=2 * self.max_cut + 1)
        # nbor_info = self.nbor_list(c_size=21.0)
        n_zero_d = nbor_info["dist"][np.nonzero(nbor_info["dist"])]
        hist, bins = np.histogram(
            n_zero_d.ravel(), bins=np.arange(0.1, 10.2, 0.1)
        )
        const = float(nbor_info["nat"]) / float(self._atoms.num_atoms)
        hist = hist / float(const)
        shell_vol = (
            4.0
            / 3.0
            * np.pi
            * (np.power(bins[1:], 3) - np.power(bins[:-1], 3))
        )
        number_density = self._atoms.num_atoms / self._atoms.volume
        rdf = (
            hist / shell_vol / number_density / self._atoms.num_atoms
        )  # /len(n_zero_d)
        nn = rdf / self._atoms.num_atoms
        if plot:
            plt.plot(bins[:-1], rdf)
            plt.savefig("rdf.png")
            plt.close()
        return bins[:-1], rdf, nn

    def get_dist_cutoffs(self):
        """
        Get different distance cut-offs.

        Args:

            s: Structure object

        Returns:
               rcut: max-cutoff to ensure all the element-combinations
               are included, used in calculating angluar distribution
               upto first neighbor

               rcut1: decide first cut-off based on total RDF and a buffer
               (previously used in dihedrals, but not used now in the code)

               rcut2: second neighbor cut-off

               rcut_dihed: specialized cut-off for dihedrals to avaoid large
               bonds such as N-N, uses average bond-distance and standard
               deviations
        """
        max_cut = self.max_cut
        x, y, z = self.get_rdf()
        arr = []
        for i, j in zip(x, z):
            if j > 0.0:
                arr.append(i)
        rcut_buffer = 0.11
        io1 = 0
        io2 = 1
        io3 = 2
        try:
            delta = arr[io2] - arr[io1]
            while delta < rcut_buffer and arr[io2] < max_cut:
                io1 = io1 + 1
                io2 = io2 + 1
                io3 = io3 + 1
                delta = arr[io2] - arr[io1]
            rcut1 = (arr[io2] + arr[io1]) / float(2.0)
        except Exception:
            print("Warning:Setting first nbr cut-off as minimum bond-angle")
            rcut1 = arr[0]
            pass
        try:
            delta = arr[io3] - arr[io2]
            while (
                delta < rcut_buffer
                and arr[io3] < max_cut
                and arr[io2] < max_cut
            ):
                io2 = io2 + 1
                io3 = io3 + 1
                delta = arr[io3] - arr[io2]
            rcut2 = float(arr[io3] + arr[io2]) / float(2.0)
        except Exception:
            print("Warning:Setting first and second nbr cut-off equal")
            print("You might consider increasing max_n parameter")
            rcut2 = rcut1
            pass
        # rcut, rcut_dihed = get_prdf(s=s)
        # rcut_dihed=min(rcut_dihed,max_dihed)
        return rcut1, rcut2

    def ang_dist(self, nbor_info={}, plot=False):
        """
        Get  angular distribution function upto first neighbor.

        Args:
            struct_info: struct information

            max_n: maximum number of neigbors

            c_size: max. cell size

            plot: whether to plot distributions


        Retruns:

            ang_hist1: Angular distribution upto first cut-off

            ang_bins1: angle bins
        """
        # rcut1,rcut2=self.get_dist_cutoffs()
        # rcut=rcut1
        # nbor_info=self.nbor_list()
        znm = 0
        nat = nbor_info["nat"]
        dist = nbor_info["dist"]
        # nn_id = nbor_info["nn_id"]
        nn = nbor_info["nn"]
        bondx = nbor_info["bondx"]
        bondy = nbor_info["bondy"]
        bondz = nbor_info["bondz"]

        ang_at = {}

        for i in range(nat):
            for in1 in range(nn[i]):
                # j1 = nn_id[in1][i]
                for in2 in range(in1 + 1, nn[i]):
                    # j2 = nn_id[in2][i]
                    nm = dist[in1][i] * dist[in2][i]
                    if nm != 0:
                        rrx = bondx[in1][i] * bondx[in2][i]
                        rry = bondy[in1][i] * bondy[in2][i]
                        rrz = bondz[in1][i] * bondz[in2][i]
                        cos = float(rrx + rry + rrz) / float(nm)
                        if cos <= -1.0:
                            cos = cos + 0.000001
                        if cos >= 1.0:
                            cos = cos - 0.000001
                        deg = math.degrees(math.acos(cos))
                        ang_at.setdefault(round(deg, 3), []).append(i)
                    else:
                        znm = znm + 1
        angs = np.array([float(i) for i in ang_at.keys()])
        norm = np.array(
            [float(len(i)) / float(len(set(i))) for i in ang_at.values()]
        )
        ang_hist, ang_bins = np.histogram(
            angs, weights=norm, bins=np.arange(1, 181.0, 1), density=False,
        )
        # if plot == True:
        #    plt.bar(ang_bins[:-1], ang_hist)
        #    plt.savefig("ang1.png")
        #    plt.close()

        return ang_hist, ang_bins

    def ang_dist_first(self, plot=False):
        """Get angular distribution upto first neighbor."""
        rcut1 = self.rcut1
        # rcut2 = self.rcut2
        nbor_info = self.nbor_list(rcut=rcut1)
        ang_hist, ang_bins = self.ang_dist(nbor_info=nbor_info)
        # print ('anghist',ang_hist)
        if plot:

            plt.plot(ang_bins[:-1], ang_hist)
            plt.savefig("adf1.png")
            plt.close()
        return ang_hist, ang_bins

    def ang_dist_second(self, plot=False):
        """Get angular distribution upto second neighbor."""
        # rcut1, rcut2 = self.get_dist_cutoffs()
        # rcut1 = self.rcut1
        rcut2 = self.rcut2
        # print ('rcut1,rcut2',rcut1,rcut2)
        nbor_info = self.nbor_list(rcut=rcut2)
        ang_hist, ang_bins = self.ang_dist(nbor_info=nbor_info)
        if plot:
            plt.plot(ang_bins[:-1], ang_hist)
            plt.savefig("adf2.png")
            plt.close()
        return ang_hist, ang_bins

    def get_ddf(self, plot=False):
        """Get dihedral distribution upto first neighbor."""
        # rcut1, rcut2 = self.get_dist_cutoffs()
        rcut1 = self.rcut1
        # rcut2 = self.rcut2
        nbor_info = self.nbor_list(rcut=rcut1)
        nat = nbor_info["nat"]
        # dist = nbor_info["dist"]
        nn_id = nbor_info["nn_id"]
        nn = nbor_info["nn"]
        bondx = nbor_info["bondx"]
        bondy = nbor_info["bondy"]
        bondz = nbor_info["bondz"]
        dih_at = {}
        for i in range(nat):
            for in1 in range(nn[i]):
                j1 = nn_id[in1][i]
                if j1 > i:
                    for in2 in range(nn[i]):
                        # all other nn of i that are not j
                        j2 = nn_id[in2][i]
                        if j2 != j1:
                            for in3 in range(
                                nn[j1]
                            ):  # all other nn of j that are not i
                                j3 = nn_id[in3][j1]
                                if j3 != i:
                                    v1 = []
                                    v2 = []
                                    v3 = []
                                    v1.append(bondx[in2][i])
                                    v1.append(bondy[in2][i])
                                    v1.append(bondz[in2][i])
                                    v2.append(-bondx[in1][i])
                                    v2.append(-bondy[in1][i])
                                    v2.append(-bondz[in1][i])
                                    v3.append(-bondx[in3][j1])
                                    v3.append(-bondy[in3][j1])
                                    v3.append(-bondz[in3][j1])
                                    v23 = np.cross(v2, v3)
                                    v12 = np.cross(v1, v2)
                                    theta = math.degrees(
                                        math.atan2(
                                            np.linalg.norm(v2)
                                            * np.dot(v1, v23),
                                            np.dot(v12, v23),
                                        )
                                    )
                                    if theta < 0.00001:
                                        theta = -theta
                                    # print "theta=",theta
                                    dih_at.setdefault(
                                        round(theta, 3), []
                                    ).append(i)
        dih = np.array([float(i) for i in dih_at.keys()])
        # dih1 = set(dih)
        # print "dih",dih1
        norm = np.array(
            [float(len(i)) / float(len(set(i))) for i in dih_at.values()]
        )

        dih_hist1, dih_bins1 = np.histogram(
            dih, weights=norm, bins=np.arange(1, 181.0, 1), density=False,
        )
        if plot:
            plt.plot(dih_bins1[:-1], dih_hist1)
            plt.savefig("dihedrals.png")
            plt.close()
        return dih_hist1, dih_bins1

    @property
    def get_all_distributions(self):
        """Get all distributions."""
        distributions = {}
        _, rdf, nn = self.get_rdf()
        ddf, _ = self.get_ddf()
        adfa, _ = self.ang_dist_first()
        adfb, _ = self.ang_dist_second()
        distributions["rdf"] = rdf
        distributions["ddf"] = ddf
        distributions["adfa"] = adfa
        distributions["adfb"] = adfb
        distributions["nn"] = nn
        return distributions

    def atomwise_radial_dist(self, rcut=10.0, c_size=0):
        """Get pair/radial distribution for each atom."""
        if rcut < c_size:
            rcut = c_size + 1
        nbor_info = self.nbor_list(rcut=rcut, c_size=c_size)
        nat = nbor_info["nat"]
        dist = nbor_info["dist"]
        atom_rdfs = []

        for i in range(nat):
            hist, bins = np.histogram(
                dist[:, i], bins=np.arange(0.1, rcut + 0.2, 0.1)
            )
            atom_rdfs.append(hist.tolist())
            if self.verbose:
                exact_dists = np.arange(0.1, rcut + 0.2, 0.1)[hist.nonzero()]
                print("exact_dists", exact_dists)
        return np.array(atom_rdfs)

    def atomwise_angle_dist(self, rcut=None, nbins=180, c_size=0):
        """Get angle distribution for each atom."""

        def angle(
            dist1, dist2, bondx1, bondx2, bondy1, bondy2, bondz1, bondz2,
        ):
            """Get an angle."""
            nm = dist1 * dist2
            rrx = bondx1 * bondx2
            rry = bondy1 * bondy2
            rrz = bondz1 * bondz2
            cos = (rrx + rry + rrz) / (nm)
            if cos <= -1.0:
                cos = cos + 0.000001
            if cos >= 1.0:
                cos = cos - 0.000001
            deg = math.degrees(math.acos(cos))
            return deg

        if rcut is None:
            rcut = self.rcut1
        nbor_info = self.nbor_list(rcut=rcut, c_size=c_size)
        atom_angles = []
        for i in range(nbor_info["nat"]):
            angles = [
                angle(
                    nbor_info["dist"][in1][i],
                    nbor_info["dist"][in2][i],
                    nbor_info["bondx"][in1][i],
                    nbor_info["bondx"][in2][i],
                    nbor_info["bondy"][in1][i],
                    nbor_info["bondy"][in2][i],
                    nbor_info["bondz"][in1][i],
                    nbor_info["bondz"][in2][i],
                )
                for in1 in range(nbor_info["nn"][i])
                for in2 in range(nbor_info["nn"][i])
                if in2 > in1
                and nbor_info["dist"][in1][i] * nbor_info["dist"][in2][i] != 0
            ]
            ang_hist, ang_bins = np.histogram(
                angles, bins=np.arange(1, nbins + 2, 1), density=False,
            )
            atom_angles.append(ang_hist)
            # print("fff",i,ang_hist)
            if self.verbose:
                exact_angles = np.arange(1, nbins + 2, 1)[ang_hist.nonzero()]
                print("exact_angles", exact_angles)
        # return (atom_angles)#/nbor_info['nat']
        return np.array(atom_angles)  # /nbor_info['nat']


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms

    box = [[5.493642, 0, 0], [0, 5.493642, 0], [0, 0, 5.493642]]
    elements = ["Si", "Si", "Si", "Si", "Si", "Si", "Si", "Si"]
    coords = [
        [0, 0, 0],
        [0.25, 0.25, 0.25],
        [0.000000, 0.500000, 0.500000],
        [0.250000, 0.750000, 0.750000],
        [0.500000, 0.000000, 0.500000],
        [0.750000, 0.250000, 0.750000],
        [0.500000, 0.500000, 0.000000],
        [0.750000, 0.750000, 0.250000],
    ]
    # box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    # coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    # elements = ["Si", "Si"]
    Si = Atoms(
        lattice_mat=box, coords=coords, elements=elements
    )  # .get_primitive_atoms
    ver = False  # True
    a = NeighborsAnalysis(Si, verbose=ver, max_cut=5).atomwise_radial_dist(
        c_size=5
    )

    print(a[1], len(a), a[1].nonzero())
    a = NeighborsAnalysis(Si, verbose=ver, max_cut=5).atomwise_radial_dist(
        c_size=10
    )
    print(a[1], len(a), a[1].nonzero())
    a = NeighborsAnalysis(Si, verbose=ver, max_cut=5).atomwise_radial_dist(
        c_size=20
    )
    print(a[1], len(a), a[1].nonzero())
"""
