"""
Classical Force-field Inspired Descriptors (CFID)
Find details in:
https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.083801
"""
from __future__ import unicode_literals, print_function
import matplotlib.pyplot as plt

plt.switch_backend("agg")
from collections import defaultdict
import itertools
from scipy.stats import gaussian_kde
from math import pi
from operator import itemgetter
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.outputs import Vasprun

# from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
import collections, math, os

# from pymatgen.analysis.defects.point_defects import ValenceIonicRadiusEvaluator
import numpy as np
from math import log
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from monty.json import MontyEncoder, MontyDecoder
from pymatgen.core.periodic_table import Element
import json, sys
from pymatgen.core.structure import Structure
import numpy as np
import xml.etree.ElementTree as ET
from monty.serialization import loadfn
from pandas import DataFrame as pd
from pymatgen.core.lattice import Lattice
from pymatgen.core import Composition
from interruptingcow import timeout

# Note time limit for angular part is hardcoded

el_chrg_json = str(os.path.join(os.path.dirname(__file__), "element_charge.json"))
el_chem_json = str(os.path.join(os.path.dirname(__file__), "Elements.json"))
timelimit = 300  # 5 minutes


def get_effective_structure(s=None, tol=8.0):
    """
  Check if there is vacuum, if so get actual size of the structure
  and the add vaccum of size tol to make sure structures
  are independent of user defined vacuum

  Args:
       s: Structure object
       tol: vacuum tolerance
  Returns:
         s: re-structure structure with tol vacuum
  """
    coords = s.cart_coords
    range_x = abs(max(coords[:, 0]) - min(coords[:, 0]))
    range_y = abs(max(coords[:, 1]) - min(coords[:, 1]))
    range_z = abs(max(coords[:, 2]) - min(coords[:, 2]))
    a = abs(s.lattice.a)  # matrix[0][0]
    b = abs(s.lattice.b)  # matrix[1][1]
    c = abs(s.lattice.c)  # matrix[2][2]
    if abs(a - range_x) > tol:
        a = range_x + tol
        print("Find vaccum in x-direction")
    if abs(b - range_y) > tol:
        b = range_y + tol
        print("Find vaccum in y-direction")
    if abs(c - range_z) > tol:
        c = range_z + tol
        print("Find vaccum in z-direction")
    arr = Lattice(
        [
            [s.lattice.matrix[0][0], s.lattice.matrix[0][1], s.lattice.matrix[0][2]],
            [s.lattice.matrix[1][0], s.lattice.matrix[1][1], s.lattice.matrix[1][2]],
            [s.lattice.matrix[2][0], s.lattice.matrix[2][1], s.lattice.matrix[2][2]],
        ]
    )
    s = Structure(arr, s.species, coords, coords_are_cartesian=True)
    return s


def el_combs(s=[]):
    """
  Get element combinations for a Structure object

  Args:
      s: Structure object
  Returns: 
         comb: combinations
    """
    symb = s.symbol_set
    tmp = map("-".join, itertools.product(symb, repeat=2))
    comb = list(set([str("-".join(sorted(i.split("-")))) for i in tmp]))
    return comb


def flatten_out(arr=[], tol=0.1):
    """
    Determine first cut-off
 
    Args:
         arr: array
         tol: toelrance
    Return:
          rcut: cut-off for a given tolerance tol, 
          because sometimes RDF peaks could be very close
    """

    rcut_buffer = tol
    io1 = 0
    io2 = 1
    io3 = 2
    delta = arr[io2] - arr[io1]

    while delta < rcut_buffer and io3 < len(arr):
        io1 = io1 + 1
        io2 = io2 + 1
        io3 = io3 + 1
        delta = arr[io2] - arr[io1]
    rcut = (arr[io2] + arr[io1]) / float(2.0)
    return rcut


def smooth_kde(x, y):
    """
    For making smooth distributions
    """
    denn = gaussian_kde(y)
    denn.covariance_factor = lambda: 0.25
    denn._compute_covariance()
    xs = np.linspace(0, max(x), 100)
    kde = denn(xs)
    return kde


def get_prdf(
    s=None, cutoff=10.0, intvl=0.1, plot_prdf=False, std_tol=0.25, filename="prdf.png"
):
    """
    Get partial radial distribution function

    Args:
         s: Structure object
         cutoff: maximum cutoff in Angstrom
         intvl: bin-size
         plot_prdf: whether to plot PRDF
         std_tol: when calculating the dihedral angles, the code got stuck for some N-containing compounds
         so we averaged all the first neghbors and added a fraction of the standard deviation as a tolerance
         this new tolerance is used for dihedral part only, while angular and RDF part remains intact
         filename: if plotted the name of the output file
    Returns:
          max-cutoff to ensure all the element-combinations are included
    """

    neighbors_lst = s.get_all_neighbors(cutoff)
    comb = el_combs(
        s=s
    )  # [str(i[0])+str('-')+str(i[1]) for i in list(itertools.product(sps, repeat=2))]
    info = {}
    for c in comb:
        for i, ii in enumerate(neighbors_lst):
            for j in ii:
                comb1 = str(s[i].specie) + str("-") + str(j[0].specie)
                comb2 = str(j[0].specie) + str("-") + str(s[i].specie)
                if comb1 == c or comb2 == c:
                    info.setdefault(c, []).append(j[1])
    plt.rcParams.update({"font.size": 22})
    for i in info.items():
        i[1].sort()
        dist_hist, dist_bins = np.histogram(
            i[1], bins=np.arange(0, cutoff + intvl, intvl), density=False
        )
        shell_vol = (
            4.0 / 3.0 * pi * (np.power(dist_bins[1:], 3) - np.power(dist_bins[:-1], 3))
        )
        number_density = s.num_sites / s.volume
        rdf = dist_hist / shell_vol / number_density / len(neighbors_lst)
        if plot_prdf == True:
            plt.plot(dist_bins[:-1], rdf, label=i[0], linewidth=2)

    if plot_prdf == True:

        plt.legend(prop={"size": 26})
        plt.xlabel("r ($\AA$)")
        plt.ylabel("g(r)")
        plt.ylim(ymin=0)
        plt.tight_layout()
        return plt
        # plt.savefig(filename)
        # plt.close()

    cut_off = {}

    for i, j in info.items():
        cut_off[i] = flatten_out(arr=j, tol=0.1)
    vals = np.array([float(i) for i in cut_off.values()])
    # print (cut_off.items(),' fwfwf',np.mean(vals)+0.25*np.std(vals))
    return (
        max(cut_off.items(), key=itemgetter(1))[1],
        np.mean(vals) + 0.25 * np.std(vals),
    )


def get_rdf(s=None, cutoff=10.0, intvl=0.1):
    """
    Get total radial distribution function

    Args:
         s: Structure object
         cutoff: maximum distance for binning
         intvl: bin-size
    Returns:
           bins, RDF, bond-order distribution
    """

    neighbors_lst = s.get_all_neighbors(cutoff)
    all_distances = np.concatenate(
        tuple(map(lambda x: [itemgetter(1)(e) for e in x], neighbors_lst))
    )
    rdf_dict = {}
    dist_hist, dist_bins = np.histogram(
        all_distances, bins=np.arange(0, cutoff + intvl, intvl), density=False
    )  # equivalent to bon-order
    shell_vol = (
        4.0 / 3.0 * pi * (np.power(dist_bins[1:], 3) - np.power(dist_bins[:-1], 3))
    )
    number_density = s.num_sites / s.volume
    rdf = dist_hist / shell_vol / number_density / len(neighbors_lst)
    return (
        dist_bins[:-1],
        [round(i, 4) for i in rdf],
        dist_hist / float(len(s)),
    )  # [{'distances': dist_bins[:-1], 'distribution': rdf}]


def get_dist_cutoffs(s="", max_cut=5.0):
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

    x, y, z = get_rdf(s)
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
    rcut, rcut_dihed = get_prdf(s=s)
    # rcut_dihed=min(rcut_dihed,max_dihed)
    return rcut, rcut1, rcut2, rcut_dihed


def get_structure_data(s="", c_size=10.0):
    """
    Convert Structure object to a represnetation where only unique sites are keps
    Used for angular distribution terms
    Args:
        s: Structure object
    Returns:
          info with coords,dim,nat,lat,new_symbs
    """
    info = {}
    coords = s.frac_coords
    box = s.lattice.matrix
    dim1 = int(float(c_size) / float(max(abs(box[0])))) + 1
    dim2 = int(float(c_size) / float(max(abs(box[1])))) + 1
    dim3 = int(float(c_size) / float(max(abs(box[2])))) + 1
    dim = [dim1, dim2, dim3]
    dim = np.array(dim)
    all_symbs = [i.symbol for i in s.species]
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


@timeout(timelimit)
def ang_dist1(struct_info={}, c_size=10.0, max_n=500, plot=True, max_cut=5.0, rcut=""):
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
    coords = struct_info["coords"]
    dim = struct_info["dim"]
    nat = struct_info["nat"]
    lat = struct_info["lat"]
    znm = 0
    bond_arr = []
    deg_arr = []
    nn = np.zeros((nat), dtype="int")
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
                if np.fabs(diff[v]) >= dim05[v]:
                    diff[v] = diff[v] - np.sign(diff[v])
            new_diff = np.dot(diff, lat)
            dd = np.linalg.norm(new_diff)
            if dd < rcut and dd >= 0.1:
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
    # print ('maxnn',max(nn))

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
    ang_hist1, ang_bins1 = np.histogram(
        angs, weights=norm, bins=np.arange(1, 181.0, 1), density=False
    )

    if plot == True:
        plt.bar(ang_bins1[:-1], ang_hist1)
        plt.savefig("ang1.png")
        plt.close()

    return ang_hist1, ang_bins1


@timeout(timelimit)
def ang_dist2(struct_info={}, c_size=10.0, max_n=500, max_cut=5.0, rcut2="", plot=True):
    """
    Get  angular distribution function upto second neighbor

    Args:
        struct_info: struct information
        max_n: maximum number of neigbors
        c_size: max. cell size
        plot: whether to plot distributions
        max_cut: max. bond cut-off for angular distribution
    Retruns:
         ang_hist2: Angular distribution upto secondt cut-off
         ang_bins2: angle bins
    """
    coords = struct_info["coords"]
    coords = struct_info["coords"]
    dim = struct_info["dim"]
    nat = struct_info["nat"]
    lat = struct_info["lat"]

    ### 2nd neighbors
    znm = 0
    bond_arr = []
    deg_arr = []
    nn = np.zeros((nat), dtype="int")
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
                if np.fabs(diff[v]) >= dim05[v]:
                    diff[v] = diff[v] - np.sign(diff[v])
            new_diff = np.dot(diff, lat)
            dd = np.linalg.norm(new_diff)
            if dd < rcut2 and dd >= 0.1:
                nn_index = nn[i]  # index of the neighbor
                nn[i] = nn[i] + 1
                dist[nn_index][i] = dd  # nn_index counter id
                nn_id[nn_index][i] = j  # exact id
                bondx[nn_index, i] = new_diff[0]
                bondy[nn_index, i] = new_diff[1]
                bondz[nn_index, i] = new_diff[2]
                nn_index1 = nn[j]  # index of the neighbor
                nn[j] = nn[j] + 1
                dist[nn_index1][j] = dd  # nn_index counter id
                nn_id[nn_index1][j] = i  # exact id
                bondx[nn_index1, j] = -new_diff[0]
                bondy[nn_index1, j] = -new_diff[1]
                bondz[nn_index1, j] = -new_diff[2]

    # print ('ang at2')
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
                    # ang_at.setdefault(deg, []).append(i)
                    ang_at.setdefault(round(deg, 3), []).append(i)
                else:
                    znm = znm + 1
    angs = np.array([float(i) for i in ang_at.keys()])
    norm = np.array([float(len(i)) / float(len(set(i))) for i in ang_at.values()])

    ang_hist2, ang_bins2 = np.histogram(
        angs, weights=norm, bins=np.arange(1, 181.0, 1), density=False
    )
    if plot == True:
        plt.bar(ang_bins2[:-1], ang_hist2)
        plt.savefig("ang2.png")
        plt.close()
    return ang_hist2, ang_bins2


@timeout(timelimit)
def dihedrals(
    struct_info={}, new_symbs=[], max_n=500, c_size=10.0, rcut_dihed="", plot=False
):
    """
    Get dihedral angle distribution
    """

    coords = struct_info["coords"]
    dim = struct_info["dim"]
    nat = struct_info["nat"]
    lat = struct_info["lat"]

    znm = 0
    bond_arr = []
    deg_arr = []
    nn = np.zeros((nat), dtype="int")
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
                if np.fabs(diff[v]) >= dim05[v]:
                    diff[v] = diff[v] - np.sign(diff[v])
            new_diff = np.dot(diff, lat)
            dd = np.linalg.norm(new_diff)
            if dd < rcut_dihed and dd >= 0.1:
                # if dd<rcut1 and dd>=0.1:
                nn_index = nn[i]  # index of the neighbor
                nn[i] = nn[i] + 1
                dist[nn_index][i] = dd  # nn_index counter id
                nn_id[nn_index][i] = j  # exact id
                bondx[nn_index, i] = new_diff[0]
                bondy[nn_index, i] = new_diff[1]
                bondz[nn_index, i] = new_diff[2]
                nn_index1 = nn[j]  # index of the neighbor
                nn[j] = nn[j] + 1
                dist[nn_index1][j] = dd  # nn_index counter id
                nn_id[nn_index1][j] = i  # exact id
                bondx[nn_index1, j] = -new_diff[0]
                bondy[nn_index1, j] = -new_diff[1]
                bondz[nn_index1, j] = -new_diff[2]

    dih_at = {}
    for i in range(nat):
        for in1 in range(nn[i]):
            j1 = nn_id[in1][i]
            if j1 > i:
                for in2 in range(nn[i]):  # all other nn of i that are not j
                    j2 = nn_id[in2][i]
                    if j2 != j1:
                        for in3 in range(nn[j1]):  # all other nn of j that are not i
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


def get_chgdescrp_arr(elm=""):
    """
      Get charge descriptors for an element

      Args:
           elm: element name
      Returns:
             arr: array value
      """
    arr = []

    try:
        f = open(el_chrg_json, "r")
        emdat = json.load(f)
        f.close()
        arr = emdat[elm][0][1]
    except:
        pass
    return arr


def get_descrp_arr_name(elm="Al"):
    """
      Get chemical descriptors for an element

      Args:
           elm: element name
      Returns:
             arr: array value
      """
    arr = []
    try:
        f = open(el_chem_json, "r")
        dat = json.load(f)
        f.close()

        d = dat[elm]
        arr = []
        for k, v in d.items():
            arr.append(k)
    except:
        pass
    return arr


def get_descrp_arr(elm=""):
    """
      Get chemical descriptors for an element

      Args:
           elm: element name
      Returns:
             arr: array value
      """
    arr = []
    try:
        f = open(el_chem_json, "r")
        dat = json.load(f)
        f.close()

        d = dat[elm]
        arr = []
        for k, v in d.items():
            arr.append(v)
        arr = np.array(arr).astype(float)
    except:
        pass
    return arr


def packing_fraction(s=None):
    """
    Get packing fraction

    Args:
         s: Structure object
    Returns:
           packing fraction
    """
    total_rad = 0
    for site in s:
        total_rad += site.specie.atomic_radius ** 3
    pf = np.array([4 * pi * total_rad / (3 * s.volume)])
    return pf[0]


def get_comp_descp(
    struct="",
    jcell=True,
    jmean_chem=True,
    jmean_chg=True,
    jrdf=False,
    jrdf_adf=True,
    print_names=False,
):
    """
        Get chemo-structural CFID decriptors

        Args:
            struct: Structure object
            jcell: whether to use cell-size descriptors
            jmean_chem: whether to use average chemical descriptors
            jmean_chg: whether to use average charge distribution descriptors
            jmean_rdf: whether to use radial distribution descriptors
            jrdf_adf: whether to use radial as well as angle distribution descriptors
            print_names: whether to print names of descriptors
        Returns:
              cat: catenated final descriptors
        """
    cat = []
    s = struct
    # s = get_effective_structure(struct)
    cell = []
    mean_chem = []
    rdf = []
    adf = []
    nn = []
    mean_chg = []
    adfa = []
    adfb = []
    ddf = []

    if jmean_chem == True:
        comp = s.composition
        el_dict = comp.get_el_amt_dict()
        arr = []
        for k, v in el_dict.items():
            des = get_descrp_arr(k)
            arr.append(des)
        mean_chem = np.mean(arr, axis=0)
        # print ('mean_chem',len(mean_chem))

    if jcell == True:
        v_pa = round(float(s.volume) / float(s.composition.num_atoms), 5)
        logv_pa = round(log(float(s.volume) / float(s.composition.num_atoms)), 5)
        pf = round(packing_fraction(s), 5)
        density = round(s.density, 5)
        cell = np.array([v_pa, logv_pa, pf, density])
        # print ('jcell',len(cell))

    if jrdf == True:
        distrdf, bins, bo = get_rdf(s=s)
        rdf = np.array(distrdf)
        # print ('rdf',len(rdf))

    if jrdf_adf == True:
        rcut, rcut1, rcut2, rcut_dihed = get_dist_cutoffs(s=s)
        # print ('rcut,rcut1,rcut2,rcut_dihed',rcut,rcut1,rcut2,rcut_dihed)
        struct_info = get_structure_data(s=s)
        # print ('struct_info',struct_info)
        bins, rdf, nn = get_rdf(s=s)
        try:
            adfa = np.zeros(179)
            adfa, _ = ang_dist1(struct_info=struct_info, plot=False, rcut=rcut_dihed)
        except:
            print(
                "Angular distribution1 part is taking too long a time,setting it to zero"
            )
            pass
        try:
            adfb = np.zeros(179)
            adfb, _ = ang_dist2(struct_info=struct_info, plot=False, rcut2=rcut2)
        except:
            print(
                "Angular distribution2 part is taking too long a time,setting it to zero"
            )
            pass
        try:
            ddf = np.zeros(179)
            ddf, _ = dihedrals(
                struct_info=struct_info, plot=False, rcut_dihed=rcut_dihed
            )
        except:
            print(
                "Dihedral distribution part is taking too long a time,setting it to zero"
            )
            pass
        adfa = np.array(adfa)
        adfb = np.array(adfb)
        rdf = np.array(rdf)
        ddf = np.array(ddf)
        nn = np.array(nn)
        # print ('adfa',len(adfa))
        # print ('ddf',len(ddf))
        # print ('adfb',len(adfb))
        # print ('rdf',len(rdf))
        # print ('nn',len(nn))

    if jmean_chg == True:
        chgarr = []
        for k, v in el_dict.items():
            chg = get_chgdescrp_arr(k)
            chgarr.append(chg)
        mean_chg = np.mean(chgarr, axis=0)
        # print ('mean_chg',len(mean_chg))

    if print_names == True:
        nmes = []
        chem_nm = get_descrp_arr_name()
        for d, nm in zip(
            [mean_chem, cell, mean_chg, rdf, adfa, adfb, ddf, nn],
            ["mean_chem", "cell", "mean_chg", "rdf", "adfa", "adfb", "ddf", "nn"],
        ):
            if d != []:
                for ff, dd in enumerate(d):
                    cat.append(dd)
                    if nm == "mean_chem":
                        tag = chem_nm[ff]
                    else:
                        tag = str(nm) + str("_") + str(ff)
                    nmes.append(str(tag))
        cat = np.array(cat).astype(float)
        # print (nmes,len(nmes))
        return nmes
    else:
        for d, nm in zip(
            [mean_chem, cell, mean_chg, rdf, adfa, adfb, ddf, nn],
            ["mean_chem", "cell", "mean_chg", "rdf", "adfa", "adfb", "ddf", "nn"],
        ):
            if d != []:
                for ff, dd in enumerate(d):
                    cat.append(dd)
        cat = np.array(cat).astype(float)
    ##except:
    # pass
    return cat


def get_chemonly(string=""):
    """
  Get only mean chemical descriptors for a chemical formula, say Al2O3

  Args:
       string: chemical formula, say Al2O3, NaCl etc.
  Returns:
         mean_chem: average chemical descriptors
  """
    comp = Composition(string)
    el_dict = comp.get_el_amt_dict()
    arr = []
    for k, v in el_dict.items():
        des = get_descrp_arr(k)
        arr.append(des)
    mean_chem = np.mean(arr, axis=0)
    return mean_chem


if __name__ == "__main__":

    # chemo-structural features
    s = Structure.from_file("POSCAR")
    ss = s.copy()
    # Making supercell
    ss.make_supercell([2, 3, 4])
    x = get_comp_descp(s)
    xx = get_comp_descp(ss)
    print("len", len(x))
    # print(len(x))
    count = 0
    for i, j in zip(x, xx):
        count = count + 1
        if i != j:
            print(count, i, j)
    # only chemical features for a structure
    y = get_comp_descp(
        struct=s,
        jcell=False,
        jmean_chem=True,
        jmean_chg=False,
        jrdf=False,
        jrdf_adf=False,
        print_names=False,
    )
    print(len(y))

    # chemical features based on composition only
    chm = get_chemonly("Al2O3")
    print(len(chm))

    # charge descriptors
    comp = Composition("Al2O3")
    el_dict = comp.get_el_amt_dict()
    chgarr = []
    for k, v in el_dict.items():
        chg = get_chgdescrp_arr(k)
        chgarr.append(chg)
    mean_chg = np.mean(chgarr, axis=0)
    print(len(mean_chg))
