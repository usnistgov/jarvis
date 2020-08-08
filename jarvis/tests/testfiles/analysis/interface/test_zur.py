from jarvis.analysis.interface.zur import (
    ZSLGenerator,
    get_hetero_type,
    make_interface,
    add_atoms,
)
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os
from jarvis.db.figshare import get_jid_data


def get_2d_hetero_jids(jid1="JVASP-664", jid2="JVASP-52"):
    from jarvis.db.figshare import get_jid_data

    m1 = get_jid_data(jid1)["atoms"]
    m2 = get_jid_data(jid2)["atoms"]
    mat1 = Atoms.from_dict(m1)
    mat2 = Atoms.from_dict(m2)
    vac = max(mat1.lattice_mat[2][2], mat2.lattice_mat[2][2])
    print("jid1", jid1)
    print(mat1)
    print("jid2", jid2)
    print(mat2)
    combined = make_interface(
        film=mat1.center_around_origin(),
        # max_area=1000,
        max_area=400,
        max_area_ratio_tol=0.09,
        # ltol=0.5,
        ltol=0.05,
        # atol=1.1,
        atol=0.1,
        subs=mat2.center_around_origin(),
    )["interface"]
    return combined


# Good 2D examples
jids = [
    "JVASP-664",
    "JVASP-667",
    "JVASP-661",
    "JVASP-688",
    "JVASP-652",
    "JVASP-676",
    "JVASP-6841",
    "JVASP-771",
    "JVASP-780",
    "JVASP-5899",
    "JVASP-646",
    "JVASP-5956",
    "JVASP-649",
    "JVASP-658",
    "JVASP-744",
    "JVASP-5872",
    "JVASP-76196",
    "JVASP-5983",
    "JVASP-20002",
    "JVASP-27836",
    "JVASP-27862",
    "JVASP-6838",
    "JVASP-687",
    "JVASP-6613",
    "JVASP-6667",
    "JVASP-6028",
    "JVASP-6079",
    "JVASP-5944",
    "JVASP-680",
    "JVASP-6955",
    "JVASP-6901",
    "JVASP-19987",
    "JVASP-6238",
    "JVASP-13600",
    "JVASP-27940",
    "JVASP-27775",
    "JVASP-27780",
    "JVASP-730",
    "JVASP-5929",
    "JVASP-19989",
    "JVASP-792",
    "JVASP-6268",
    "JVASP-6181",
    "JVASP-12027",
    "JVASP-27785",
    "JVASP-28152",
    "JVASP-28268",
    "JVASP-6994",
    "JVASP-31373",
    "JVASP-762",
    "JVASP-756",
    "JVASP-27724",
    "JVASP-13514",
    "JVASP-19987",
    "JVASP-31368",
    "JVASP-6007",
    "JVASP-741",
    "JVASP-720",
    "JVASP-789",
    "JVASP-777",
    "JVASP-5932",
    "JVASP-728",
    "JVASP-673",
    "JVASP-60477",
    "JVASP-6742",
    "JVASP-12064",
    "JVASP-6034",
    "JVASP-5935",
    "JVASP-60525",
    "JVASP-14431",
    "JVASP-13632",
    "JVASP-31379",
    "JVASP-28106",
    "JVASP-27864",
    "JVASP-27855",
    "JVASP-5950",
    "JVASP-696",
    "JVASP-783",
    "JVASP-13526",
    "JVASP-13541",
    "JVASP-31353",
    "JVASP-5926",
    "JVASP-5959",
    "JVASP-726",
    "JVASP-60244",
    "JVASP-13536",
    "JVASP-27940",
    "JVASP-60497",
    "JVASP-705",
    "JVASP-27978",
    "JVASP-27865",
    "JVASP-5977",
    "JVASP-5974",
    "JVASP-750",
    "JVASP-738",
    "JVASP-655",
    "JVASP-31356",
    "JVASP-6862",
    "JVASP-765",
    "JVASP-13586",
]


def test_zur():
    m1 = get_jid_data("JVASP-664")["atoms"]
    m2 = get_jid_data("JVASP-652")["atoms"]
    s1 = Atoms.from_dict(m1)
    s2 = Atoms.from_dict(m2)
    info = make_interface(film=s1, subs=s2)
    combined = info["interface"]
    assert (
        round(info["mismatch_u"], 3),
        round(info["mismatch_angle"], 3),
    ) == (-0.041, 0.0,)


def test_mos2_bn():
    intf = get_2d_hetero_jids(jid1="JVASP-688", jid2="JVASP-664")
    print(intf)
    print(intf.density)
    assert round(intf.density, 4) == round(1.7346542883452571, 4)


def test_2d_interface():
    jids = [
        "JVASP-649",
        "JVASP-652",
        "JVASP-658",
        "JVASP-664",
        "JVASP-688",
        "JVASP-6841",
        "JVASP-5983",
        "JVASP-60244",
        "JVASP-60389",
        "JVASP-76195",
    ]
    # jids = ["JVASP-652", "JVASP-664"]
    jids = ["JVASP-76195", "JVASP-5983"]
    jids = ["JVASP-76195", "JVASP-60389"]
    # jids = ["JVASP-76195", "JVASP-60244"]
    jids = ["JVASP-688", "JVASP-664"]
    jids = ["JVASP-652", "JVASP-646"]
    jids = ["JVASP-652", "JVASP-658"]  # WS2,WSe2
    count = 0
    for i in jids:
        for j in jids:
            if count < 100 and i != j:
                # try:
                intf = get_2d_hetero_jids(jid1=i, jid2=j)
                if intf.num_atoms < 200:

                    count = count + 1
                    print(i, j)
                    p = Poscar(intf)
                    filename = "POSCAR-" + str(i) + "_" + str(j) + ".vasp"
                    p.comment = "Surf@" + str(i) + "_" + str(j)
                    print(p)
                    # p.write_file(filename)
                    print()
                    print()
                    print()
            # except:
            #  pass


def test_type():
    B = {}
    B["scf_vbm"] = -5
    B["scf_cbm"] = -4
    B["avg_max"] = -2
    A = {}
    A["scf_vbm"] = -7
    A["scf_cbm"] = -6
    A["avg_max"] = -2
    int_type, stack = get_hetero_type(A=A, B=B)
    assert (int_type, stack) == ("III", "BA")
    print(int_type, stack)
    int_type, stack = get_hetero_type(A=B, B=A)
    print(int_type, stack)

    B = {}
    B["scf_vbm"] = -5
    B["scf_cbm"] = -4
    B["avg_max"] = -1
    A = {}
    A["scf_vbm"] = -7
    A["scf_cbm"] = -3
    A["avg_max"] = -1
    int_type, stack = get_hetero_type(A=B, B=A)
    print(int_type, stack)

    B = {}
    B["scf_vbm"] = -4
    B["scf_cbm"] = -2
    B["avg_max"] = -1
    A = {}
    A["scf_vbm"] = -7
    A["scf_cbm"] = -3
    A["avg_max"] = -1
    int_type, stack = get_hetero_type(A=B, B=A)
    print(int_type, stack)


jids = [
    "JVASP-649",
    "JVASP-652",
    "JVASP-658",
    "JVASP-664",
    "JVASP-688",
    "JVASP-6841",
    "JVASP-5983",
    "JVASP-60244",
    "JVASP-60389",
    "JVASP-76195",
]


def test_2d_hetero():
    count = 0
    for i in jids:
        for j in jids:
            if count < 10 and i != j:
                try:
                    intf = get_2d_hetero_jids(jid1=i, jid2=j)
                    ats = intf.get_string(cart=False)
                    if intf.num_atoms < 200:
                        count = count + 1
                        print(i, j)
                        print(ats)
                        print()
                        print()
                        print()
                except:
                    pass


def test_metal_ceramic_interface():
    m1 = get_jid_data(jid="JVASP-816", dataset="dft_3d")["atoms"]
    m2 = get_jid_data(jid="JVASP-32", dataset="dft_3d")["atoms"]
    mat_Al = Atoms.from_dict(m1)
    mat_Al2O3 = Atoms.from_dict(m2)
    from jarvis.analysis.defects.surface import Surface

    mat1 = Surface(atoms=mat_Al, indices=[1, 1, 1], layers=3).make_surface()
    mat2 = Surface(atoms=mat_Al2O3, indices=[0, 0, 1], layers=1).make_surface()
    combined = make_interface(
        film=mat1,  # .center_around_origin(),
        max_area=500,
        max_area_ratio_tol=0.09,
        ltol=0.01,
        apply_strain=True,
        subs=mat2,  # .center_around_origin(),
    )["interface"]
    print(combined)
    # from ase.lattice.surface import surface
    # from pymatgen.io.ase import AseAtomsAdaptor
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    # from pymatgen.io.vasp.inputs import Poscar
    # mat_cvn = SpacegroupAnalyzer(mat_Al.pymatgen_converter()).get_conventional_standard_structure()
    # ase_atoms = AseAtomsAdaptor().get_atoms(mat_cvn)
    # ase_slab = surface(ase_atoms, [1,1,1], 3)
    # ase_slab.center(vacuum=18, axis=2)
    # slab_pymatgen = AseAtomsAdaptor().get_structure(ase_slab)
    # slab_pymatgen.sort()
    # print (Poscar(slab_pymatgen))
    # print ()
    # print ()
    print(mat1.center_around_origin().get_string(cart=False))


# test_mos2_bn()
# test_2d_interface()
# test_metal_ceramic_interface()
