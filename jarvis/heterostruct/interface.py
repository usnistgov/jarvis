from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from matplotlib.colors import LinearSegmentedColormap
from math import floor
from monty.serialization import loadfn, MontyDecoder
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.substrate_analyzer import ZSLGenerator, SubstrateAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from ase.lattice.surface import surface
from pymatgen.core.surface import Slab, SlabGenerator
from jarvis.db.static.explore_db import get_3d_dataset, get_2d_dataset

plt.switch_backend("agg")


def center_struct(s=""):
    coords = s.cart_coords
    avg_x = np.sum(coords[:, 0]) / len(s)
    avg_y = np.sum(coords[:, 1]) / len(s)
    avg_z = np.sum(coords[:, 2]) / len(s)
    coords = s.cart_coords - [avg_x, avg_y, avg_z]
    new_s = Structure(s.lattice, s.species, coords, coords_are_cartesian=True)
    return new_s


def rangexyz(s=""):
    coords = s.cart_coords
    range_x = abs(max(coords[:, 0]) - min(coords[:, 0]))
    range_y = abs(max(coords[:, 1]) - min(coords[:, 1]))
    range_z = abs(max(coords[:, 2]) - min(coords[:, 2]))
    return range_z


def get_ase_surf(s=[], miller=(0, 0, 1), layers=1):
    strt = center_struct(s)
    ase_atoms = AseAtomsAdaptor().get_atoms(strt)
    ase_slab = surface(ase_atoms, miller, layers)
    ase_slab.center(axis=2)
    tmp = AseAtomsAdaptor().get_structure(ase_slab)
    range_z = rangexyz(strt)
    # print ('rz',abs(strt.lattice.matrix[2][2]),range_z)
    if abs(strt.lattice.matrix[2][2] - range_z) <= 6.0:
        ase_slab.center(axis=2, vacuum=12.0)
    else:
        ase_slab.center(axis=2)
    slab_pymatgen = AseAtomsAdaptor().get_structure(ase_slab)
    slab_pymatgen.sort()
    print("Surface")
    print(Poscar(slab_pymatgen))
    return slab_pymatgen


def mismatch_strts(film=[], subs=[]):

    z = ZSLGenerator()
    matches = list(z(film.lattice.matrix[:2], subs.lattice.matrix[:2], lowest=True))
    info = {}
    info["mismatch_u"] = "na"
    info["mismatch_v"] = "na"
    info["mismatch_angle"] = "na"
    info["area1"] = "na"
    info["area2"] = "na"
    info["film_sl"] = film
    info["subs_sl"] = subs

    try:
        uv1 = matches[0]["sub_sl_vecs"]
        uv2 = matches[0]["film_sl_vecs"]
        # uv1=[[-8.52917200e+00,  9.12745800e+00,  3.66344517e-17],[1.27937580e+01, 9.12745800e+00, 1.34228735e-15]]
        # uv2=[[7.02403800e+00, 1.05360570e+01, 1.07524571e-15],[-1.40480760e+01,  7.02403800e+00, -4.30098283e-16]]
        u = np.array(uv1)
        v = np.array(uv2)
        a1 = u[0]
        a2 = u[1]
        b1 = v[0]
        b2 = v[1]
        mismatch_u = np.linalg.norm(b1) / np.linalg.norm(a1) - 1
        mismatch_v = np.linalg.norm(b2) / np.linalg.norm(a2) - 1
        angle1 = (
            np.arccos(np.dot(a1, a2) / np.linalg.norm(a1) / np.linalg.norm(a2))
            * 180
            / np.pi
        )
        angle2 = (
            np.arccos(np.dot(b1, b2) / np.linalg.norm(b1) / np.linalg.norm(b2))
            * 180
            / np.pi
        )
        mismatch_angle = abs(angle1 - angle2)
        area1 = np.linalg.norm(np.cross(a1, a2))
        area2 = np.linalg.norm(np.cross(b1, b2))
        uv_substrate = uv1
        uv_mat2d = uv2
        substrate_latt = Lattice(
            np.array(
                [uv_substrate[0][:], uv_substrate[1][:], subs.lattice.matrix[2, :]]
            )
        )
        # to avoid numerical issues with find_mapping
        mat2d_fake_c = (
            film.lattice.matrix[2, :] / np.linalg.norm(film.lattice.matrix[2, :]) * 5.0
        )
        mat2d_latt = Lattice(np.array([uv_mat2d[0][:], uv_mat2d[1][:], mat2d_fake_c]))
        mat2d_latt_fake = Lattice(
            np.array(
                [film.lattice.matrix[0, :], film.lattice.matrix[1, :], mat2d_fake_c]
            )
        )
        _, __, scell = subs.lattice.find_mapping(substrate_latt, ltol=0.05, atol=1)
        scell[2] = np.array([0, 0, 1])
        subs.make_supercell(scell)
        _, __, scell = mat2d_latt_fake.find_mapping(mat2d_latt, ltol=0.05, atol=1)
        scell[2] = np.array([0, 0, 1])
        film.make_supercell(scell)
        # modify the substrate lattice so that the 2d material can be
        # grafted on top of it
        lmap = Lattice(
            np.array(
                [
                    subs.lattice.matrix[0, :],
                    subs.lattice.matrix[1, :],
                    film.lattice.matrix[2, :],
                ]
            )
        )
        film.modify_lattice(lmap)
        # print ("film",film)
        # print ("subs",subs)

        info["mismatch_u"] = mismatch_u
        info["mismatch_v"] = mismatch_v
        info["mismatch_angle"] = mismatch_angle
        info["area1"] = area1
        info["area2"] = area2
        info["film_sl"] = film
        info["subs_sl"] = subs
    except:
        pass
    return info


def get_hetero_type(A={}, B={}):
    stack = "na"
    int_type = "na"
    try:
        # if A['phi']>B['phi']:
        if A["scf_vbm"] - A["avg_max"] < B["scf_vbm"] - B["avg_max"]:
            stack = "BA"
        else:
            C = A
            D = B
            A = D
            B = C
            stack = "AB"
            # tmp=B
            # B=A
            # A=tmp
        vbm_a = A["scf_vbm"] - A["avg_max"]
        vbm_b = B["scf_vbm"] - B["avg_max"]
        cbm_a = A["scf_cbm"] - A["avg_max"]
        cbm_b = B["scf_cbm"] - B["avg_max"]
        #  print ('vbm_a,vbm_b,cbm_b,cbm_a',vbm_a,vbm_b,cbm_b,cbm_a)
        if vbm_a < vbm_b and vbm_b < cbm_b and cbm_b < cbm_a:
            int_type = "I"
        elif vbm_a < vbm_b and vbm_b < cbm_a and cbm_a < cbm_b:
            int_type = "II"
        elif vbm_a < cbm_a and cbm_a < vbm_b and vbm_b < cbm_b:
            int_type = "III"
    except:
        pass
    return int_type, stack


def get_hetero(film, substrate, seperation=3.0):

    # unique site coordinates in the substrate top layers
    coords_uniq_sub = substrate.cart_coords

    # unique site coordinates in the 2D material bottom layers
    coords_uniq_film = film.cart_coords

    substrate_top_z = max(substrate.cart_coords[:, 2])

    film_bottom = min(film.cart_coords[:, 2])

    # shift normal to the surface by 'seperation'
    sub_z = substrate.lattice.matrix[2, :]
    origin = np.array([0, 0, substrate_top_z])
    shift_normal = sub_z / np.linalg.norm(sub_z) * seperation
    # generate all possible interfaces, one for each combination of
    # unique substrate and unique 2d materials site in the layers .i.e
    # an interface structure for each parallel shift
    # interface = 2D material + substrate
    interface = substrate.copy()
    hetero_interfaces = []

    for site in film:
        new_coords = site.coords
        new_coords[2] = site.coords[2] - film_bottom
        new_coords = new_coords + origin + shift_normal
        interface.append(site.specie, new_coords, coords_are_cartesian=True)

    return interface


def get_direct_hetero(
    bulk_film="",
    bulk_subs="",
    film_miller=(1, 1, 1),
    film_subs=(0, 0, 1),
    seperation=3.0,
    film_layers=2,
    bulk_layers=1,
):
    bulk_film_cvn = SpacegroupAnalyzer(bulk_film).get_conventional_standard_structure()
    bulk_subs_cvn = SpacegroupAnalyzer(bulk_subs).get_conventional_standard_structure()
    slab_film = get_ase_surf(s=bulk_film_cvn, miller=film_miller, layers=film_layers)
    slab_subs = get_ase_surf(s=bulk_subs_cvn, miller=film_subs, layers=bulk_layers)
    hetero = get_hetero(slab_film, slab_subs, seperation=seperation)

    return hetero


def get_strt(jid=""):
    dat_2d = get_2d_dataset()
    for i in dat_2d:
        if i["jid"] == jid:
            return i["final_str"]


def make_2d_jid_strts(jid1="", jid2="", tol=3.0):
    mat1 = get_strt(jid1)
    mat2 = get_strt(jid2)
    rz1 = rangexyz(s=mat1)
    rz2 = rangexyz(s=mat2)
    x = get_direct_hetero(
        bulk_film=mat1,
        bulk_subs=mat2,
        film_miller=(0, 0, 1),
        film_subs=(0, 0, 1),
        film_layers=1,
        bulk_layers=1,
        seperation=tol,
    )
    return x
    # return x


if __name__ == "__main__":

    s = Structure.from_file("POSCAR-Al.vasp")
    x = get_direct_hetero(bulk_film=s, bulk_subs=s)
    print(Poscar(x))

    # x.to(fmt="poscar", filename="POSCAR-Al.vasp")
    x = make_2d_jid_strts(jid1="JVASP-664", jid2="JVASP-652")
    print(Poscar(x))

    s1 = Structure.from_file("POSCAR-Al.vasp")
    s2 = Structure.from_file("POSCAR-Al2O3.vasp")
    x = get_direct_hetero(bulk_film=s1, bulk_subs=s2)
    print(Poscar(x))
