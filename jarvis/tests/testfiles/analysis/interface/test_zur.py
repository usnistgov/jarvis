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
    combined = make_interface(film=mat1, subs=mat2)["interface"]
    return combined


s1 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-652"))
s2 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-664"))
s3 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-688"))
s4 = Poscar.from_file(os.path.join(os.path.dirname(__file__), "POSCAR-JVASP-75175"))

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
    info = make_interface(film=s1.atoms, subs=s2.atoms)
    combined = info["interface"]
    assert (round(info["mismatch_u"], 3), round(info["mismatch_angle"], 3),) == (
        0.043,
        0.0,
    )


def test_2d_interface():
    jids = [
        "JVASP-688",
        "JVASP-664",
        "JVASP-652",
        "JVASP-6841",
        "JVASP-5983",
        "JVASP-60244",
    ]
    jids = ["JVASP-652", "JVASP-664"]
    count = 0
    for i in jids:
        for j in jids:
            if count < 100 and i != j:
                # try:
                intf = get_2d_hetero_jids(jid1=i, jid2=j)
                ats = intf.get_string(cart=False)
                if intf.num_atoms < 200:

                    count = count + 1
                    print(i, j)
                    print(ats)
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
