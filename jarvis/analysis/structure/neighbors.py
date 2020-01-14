# https://github.com/acadien/matcalc/blob/037d7080c07856877d7a3f5f9dcbb2dec5f38dd1/analysis/rdf.py
# /rk2/knc6/UFHPC_3_16_2016/scratch-uf/LOREN/spicykamal/TaoLammps/COMB3_v18/src/analysis_tools.f
"""
This module provides classes to specify atomic structure
"""
import numpy as np
from jarvis.core.atoms import Atoms
import matplotlib.pyplot as plt
from operator import itemgetter
from jarvis.io.vasp.inputs import Poscar
plt.switch_backend("agg")
import math


class NeighborsAnalysis(object):
    def __init__(self, atoms=None):

        self._atoms = atoms
        """
        >>> box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
        >>> coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
        >>> elements = ["Si", "Si"]
        >>> Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
        >>> distributions = NeighborsAnalysis(Si).get_all_distributions
        >>> distributions['rdf']
        1
        """

    def get_structure_data(self, c_size=10.0):

        info = {}
        coords = self._atoms.frac_coords
        box = self._atoms.lattice_mat
        # print ('coords',coords)
        # print ('box',box)
        all_symbs = self._atoms.elements
        dim1 = int(float(c_size) / float(max(abs(box[0])))) + 1
        dim2 = int(float(c_size) / float(max(abs(box[1])))) + 1
        dim3 = int(float(c_size) / float(max(abs(box[2])))) + 1
        dim = [dim1, dim2, dim3]
        dim = np.array(dim)
        nat = len(coords)

        new_nat = nat * dim[0] * dim[1] * dim[2]
        new_coords = np.zeros((new_nat, 3))
        new_symbs = []  # np.chararray((new_nat))

        count = 0
        for i in range(nat):
            for j in range(dim[0]):
                for k in range(dim[1]):
                    for l in range(dim[2]):
                        new_coords[count][0] = (coords[i][0] + j) / float(dim[0])
                        new_coords[count][1] = (coords[i][1] + k) / float(dim[1])
                        new_coords[count][2] = (coords[i][2] + l) / float(dim[2])
                        new_symbs.append(all_symbs[i])
                        count = count + 1

        # print ('here4')
        nat = new_nat
        coords = new_coords

        nat = len(coords)  # int(s.composition.num_atoms)
        lat = np.zeros((3, 3))
        lat[0][0] = dim[0] * box[0][0]
        lat[0][1] = dim[0] * box[0][1]
        lat[0][2] = dim[0] * box[0][2]
        lat[1][0] = dim[1] * box[1][0]
        lat[1][1] = dim[1] * box[1][1]
        lat[1][2] = dim[1] * box[1][2]
        lat[2][0] = dim[2] * box[2][0]
        lat[2][1] = dim[2] * box[2][1]
        lat[2][2] = dim[2] * box[2][2]

        info["coords"] = coords
        info["dim"] = dim
        info["nat"] = nat
        info["lat"] = lat
        info["new_symbs"] = new_symbs
        return info

    def nbor_list(self, max_n=500, rcut=10.0, c_size=12.0):
        nbor_info = {}
        struct_info = self.get_structure_data(c_size)
        coords = np.array(struct_info["coords"])
        dim = struct_info["dim"]
        nat = struct_info["nat"]
        new_symbs = struct_info['new_symbs']
        lat = np.array(struct_info["lat"])
        znm = 0
        bond_arr = []
        deg_arr = []
        different_bond={}
        nn = np.zeros((nat), dtype="int")
        # print ('max_n, nat',max_n, nat)
        dist = np.zeros((max_n, nat))
        nn_id = np.zeros((max_n, nat), dtype="int")
        bondx = np.zeros((max_n, nat))
        bondy = np.zeros((max_n, nat))
        bondz = np.zeros((max_n, nat))
        dim05 = [float(1 / 2.0) for i in dim]
        for i in range(nat):
            for j in range(i + 1, nat):
                diff = coords[i] - coords[j]
                for v in range(3):
                    # if np.fabs(diff[v]) > dim05[v]:
                    if np.fabs(diff[v]) >= dim05[v]:
                        diff[v] = diff[v] - np.sign(diff[v])
                new_diff = np.dot(diff, lat)
                dd = np.linalg.norm(new_diff)
                if dd < rcut and dd >= 0.1:
                    sumb_i=new_symbs[i]
                    sumb_j=new_symbs[j]
                    comb = '_'.join(sorted(str(sumb_i+'_'+sumb_j).split('_')))
                    different_bond.setdefault(comb,[]).append(dd)

                    # print ('dd',dd)
                    nn_index = nn[i]  # index of the neighbor
                    nn[i] = nn[i] + 1
                    dist[nn_index][i] = dd  # nn_index counter id
                    nn_id[nn_index][i] = j  # exact id
                    bondx[nn_index][i] = new_diff[0]
                    bondy[nn_index][i] = new_diff[1]
                    bondz[nn_index][i] = new_diff[2]
                    nn_index1 = nn[j]  # index of the neighbor
                    nn[j] = nn[j] + 1
                    dist[nn_index1][j] = dd  # nn_index counter id
                    nn_id[nn_index1][j] = i  # exact id
                    bondx[nn_index1][j] = -new_diff[0]
                    bondy[nn_index1][j] = -new_diff[1]
                    bondz[nn_index1][j] = -new_diff[2]

        nbor_info["dist"] = dist
        nbor_info["nat"] = nat
        nbor_info["nn_id"] = nn_id
        nbor_info["nn"] = nn
        nbor_info["bondx"] = bondx
        nbor_info["bondy"] = bondy
        nbor_info["bondz"] = bondz
        print ('nat',nat)

        return nbor_info

    def get_rdf(self, plot=True):
        nbor_info = self.nbor_list(c_size=21.0)
        # print (nbor_info['dist'].tolist())
        n_zero_d = nbor_info["dist"][np.nonzero(nbor_info["dist"])]
        #print ('n_zero_d',n_zero_d)
        hist, bins = np.histogram(
            n_zero_d.ravel(), bins=np.arange(0.1, 10.1, 0.1)
        )
        const = float(nbor_info['nat'])/float(self._atoms.num_atoms)
        hist = hist /float(const)
        print ('nbor_info,num_atoms',nbor_info['nat'],self._atoms.num_atoms,const)
        print ('our_hist',hist)
        #print("bins[:-1]", bins[:-1])
        #print("bins[1:]", bins[1:])
        shell_vol = 4.0 / 3.0 * np.pi * (np.power(bins[1:], 3) - np.power(bins[:-1], 3))
        number_density = self._atoms.num_atoms / self._atoms.volume
        rdf =  hist / shell_vol / number_density/self._atoms.num_atoms#/len(n_zero_d)
        #rdf = 2*len(bins) * hist / np.sum(hist) / shell_vol / number_density
        #rdf = 2*len(bins) * hist / np.sum(hist) / shell_vol / number_density
        nn = rdf /len(n_zero_d)# self._atoms.num_atoms
        if plot:
            plt.bar(bins[:-1], rdf, width=0.1)
            plt.savefig("rdf.png")
            plt.close()
        return bins[:-1], rdf, nn

    def get_dist_cutoffs(self, max_cut=5.0):
        """
    Helper function to get dufferent type of distance cut-offs
    Args:
        s: Structure object
    Returns:
           rcut: max-cutoff to ensure all the element-combinations are included, used in calculating angluar distribution upto first neighbor
           rcut1: decide first cut-off based on total RDF and a buffer (previously used in dihedrals, but not used now in the code)
           rcut2: second neighbor cut-off
           rcut_dihed: specialized cut-off for dihedrals to avaoid large bonds such as N-N, uses average bond-distance and standard deviations
    """

        x, y, z = self.get_rdf()
        arr = []
        for i, j in zip(x, z):
            if j > 0.0:
                arr.append(i)
        rcut_buffer = 0.11
        io1 = 0
        io2 = 1
        io3 = 2
        delta = arr[io2] - arr[io1]
        while delta < rcut_buffer and arr[io2] < max_cut:
            io1 = io1 + 1
            io2 = io2 + 1
            io3 = io3 + 1
            delta = arr[io2] - arr[io1]
        rcut1 = (arr[io2] + arr[io1]) / float(2.0)

        delta = arr[io3] - arr[io2]
        while delta < rcut_buffer and arr[io3] < max_cut and arr[io2] < max_cut:
            io2 = io2 + 1
            io3 = io3 + 1
            delta = arr[io3] - arr[io2]
        rcut2 = float(arr[io3] + arr[io2]) / float(2.0)
        # rcut, rcut_dihed = get_prdf(s=s)
        # rcut_dihed=min(rcut_dihed,max_dihed)
        return rcut1, rcut2

    def ang_dist(self, nbor_info={}, plot=False):
        """
    Get  angular distribution function upto first neighbor
    Args:
        struct_info: struct information
        max_n: maximum number of neigbors
        c_size: max. cell size
        plot: whether to plot distributions
        max_cut: max. bond cut-off for angular distribution
    Retruns:
         ang_hist1: Angular distribution upto first cut-off
         ang_bins1: angle bins
    """

        # rcut1,rcut2=self.get_dist_cutoffs()
        # rcut=rcut1
        # nbor_info=self.nbor_list()
        nat = nbor_info["nat"]
        dist = nbor_info["dist"]
        nn_id = nbor_info["nn_id"]
        nn = nbor_info["nn"]
        bondx = nbor_info["bondx"]
        bondy = nbor_info["bondy"]
        bondz = nbor_info["bondz"]

        ang_at = {}

        for i in range(nat):
            for in1 in range(nn[i]):
                j1 = nn_id[in1][i]
                for in2 in range(in1 + 1, nn[i]):
                    j2 = nn_id[in2][i]
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
        norm = np.array([float(len(i)) / float(len(set(i))) for i in ang_at.values()])
        ang_hist, ang_bins = np.histogram(
            angs, weights=norm, bins=np.arange(1, 181.0, 1), density=False
        )
        # if plot == True:
        #    plt.bar(ang_bins[:-1], ang_hist)
        #    plt.savefig("ang1.png")
        #    plt.close()

        return ang_hist, ang_bins

    def ang_dist_first(self, plot=True):
        rcut1, rcut2 = self.get_dist_cutoffs()
        nbor_info = self.nbor_list(rcut=rcut1)
        ang_hist, ang_bins = self.ang_dist(nbor_info=nbor_info)
        # print ('anghist',ang_hist)
        if plot:

            plt.bar(ang_bins[:-1], ang_hist)
            plt.savefig("adf1.png")
            plt.close()
        return ang_hist, ang_bins

    def ang_dist_second(self, plot=True):
        rcut1, rcut2 = self.get_dist_cutoffs()
        # print ('rcut1,rcut2',rcut1,rcut2)
        nbor_info = self.nbor_list(rcut=rcut2)
        ang_hist, ang_bins = self.ang_dist(nbor_info=nbor_info)
        if plot:
            plt.bar(ang_bins[:-1], ang_hist)
            plt.savefig("adf2.png")
            plt.close()
        return ang_hist, ang_bins

    def get_ddf(self, plot=True):
        rcut1, rcut2 = self.get_dist_cutoffs()
        nbor_info = self.nbor_list(rcut=rcut1)
        nat = nbor_info["nat"]
        dist = nbor_info["dist"]
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
                    for in2 in range(nn[i]):  # all other nn of i that are not j
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
                                            np.linalg.norm(v2) * np.dot(v1, v23),
                                            np.dot(v12, v23),
                                        )
                                    )
                                    if theta < 0.00001:
                                        theta = -theta
                                    # print "theta=",theta
                                    dih_at.setdefault(round(theta, 3), []).append(i)
        dih = np.array([float(i) for i in dih_at.keys()])
        dih1 = set(dih)
        # print "dih",dih1
        norm = np.array([float(len(i)) / float(len(set(i))) for i in dih_at.values()])

        dih_hist1, dih_bins1 = np.histogram(
            dih, weights=norm, bins=np.arange(1, 181.0, 1), density=False
        )
        if plot == True:
            plt.bar(dih_bins1[:-1], dih_hist1)
            plt.savefig("dihedrals.png")
            plt.close()
        return dih_hist1, dih_bins1

    @property
    def get_all_distributions(self):
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


if __name__ == "__main__":
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
    #box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    #coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    #elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    #Si = Poscar.from_file('/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO/MAIN-RELAX-bulk@mp_541837/POSCAR').atoms
    Si = Poscar.from_file('/rk2/knc6/JARVIS-DFT/Mr-French/mp-672232.vasp_PBEBO/MAIN-ELASTIC-bulk@mp-672232/POSCAR').atoms
    x = NeighborsAnalysis(Si).get_rdf()
    #print(x)
    #import sys
    # distributions = NeighborsAnalysis(Si).get_all_distributions
    s = Si.pymatgen_converter()
    neighbors_list = s.get_all_neighbors(12.0)
    all_distances = np.concatenate(
        tuple(map(lambda x: [itemgetter(1)(e) for e in x], neighbors_list))
    )
    rdf_dict = {}
    cutoff = 10.0
    intvl = 0.1
    dist_hist, dist_bins = np.histogram(
        all_distances, bins=np.arange(0, cutoff + intvl, intvl), density=False
    )  # equivalent to bon-order
    shell_vol = (
        4.0 / 3.0 * np.pi * (np.power(dist_bins[1:], 3) - np.power(dist_bins[:-1], 3))
    )
    print ('pmg',dist_hist)
    number_density = s.num_sites / s.volume
    rdf = dist_hist / shell_vol / number_density / len(neighbors_list)
    plt.plot(dist_bins[:-1],rdf)
    plt.savefig('pmgrdf.png')
    plt.close()
    sys.exit()
    # print ('shell_vol',shell_vol)
    # print ('all_distances',all_distances)
    pmg=tuple(map(lambda x:[itemgetter(1)(e) for e in x],neighbors_list))
    our=NeighborsAnalysis(Si).nbor_list()['dist']
    print (pmg,len(pmg))
    print ()
    print (our,len(our))
    # print(distributions['rdf'])
    _, Nb = NeighborsAnalysis(Si).ang_dist_first()
    _, Nb = NeighborsAnalysis(Si).ang_dist_second()
    _, Nb = NeighborsAnalysis(Si).get_ddf()
"""
# https://github.com/caspervdw/rdf/blob/master/rdfmain.py

# if __name__=='__main__':
#    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
#    coords = [[0, 0, 0], [0.25, 0.2, 0.25]]
#    elements = ["Si", "Si"]
#    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
#    comp=Si.composition
#    #print (comp,Si.density)
#    print (Si.atomic_numbers)
#   print (Si.pymatgen_converter().composition.weight,Si.composition.weight,Si.density)
"""
