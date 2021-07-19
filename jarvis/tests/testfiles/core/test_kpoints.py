from jarvis.core.kpoints import (
    Kpoints3D,
    generate_kpath,
    generate_kgrid,
    HighSymmetryKpoint3DFactory,
)
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.core.atoms import Atoms
from jarvis.io.vasp.inputs import Poscar
import os
import tempfile


s1 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR"
    )
).atoms
s2 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-Aem2"
    )
).atoms
s3 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-C2m"
    )
).atoms
s4 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-Cmcm"
    )
).atoms
s5 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-P-1"
    )
).atoms
s6 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "analysis",
        "structure",
        "POSCAR-tetragonal",
    )
).atoms
s7 = Poscar.from_file(
    os.path.join(
        os.path.dirname(__file__), "..", "analysis", "structure", "POSCAR-Pc"
    )
).atoms


def test_kp():
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    lattice_mat = Si.lattice_mat
    kp = Kpoints3D().automatic_length_mesh(lattice_mat=lattice_mat, length=20)
    td = kp.to_dict()
    fd = Kpoints3D.from_dict(td)
    sym = kp.high_symm_path(Si)._path
    x, y = kp.interpolated_points(Si)
    kpath = generate_kpath(kpath=[[0, 0, 0], [0, 0.5, 0.5]], num_k=5)
    kps = generate_kgrid(grid=[5, 5, 5])
    new_file, filename = tempfile.mkstemp()
    kp.write_file(filename)
    sym = kp.high_symm_path(s1)._path
    sym = kp.high_symm_path(s2)._path
    sym = kp.high_symm_path(s3)._path
    sym = kp.high_symm_path(s4)._path
    sym = kp.high_symm_path(s5)._path
    sym = kp.high_symm_path(s6)._path
    sym = kp.high_symm_path(s7)._path
    assert (len(x), x[0][0], y[0], kpath[0][0], kps[0][0]) == (
        166,
        0,
        "\Gamma",
        0.0,
        0.0,
    )
    kp = Kpoints3D().kpath(atoms=Si, unique_kp_only=True)


def test_extra_spgs():
    from jarvis.db.figshare import data

    #d = data("dft_3d")
    few_spgs = {
        "JVASP-4663": 119,
        "JVASP-4666": 194,
        "JVASP-588": 160,
        "JVASP-4669": 191,
        "JVASP-4672": 141,
        "JVASP-581": 164,
        "JVASP-4687": 139,
        "JVASP-4693": 62,
        "JVASP-4696": 187,
        "JVASP-4711": 225,
        "JVASP-4714": 2,
        "JVASP-4723": 60,
        "JVASP-32": 167,
        "JVASP-107": 186,
        "JVASP-4756": 63,
        "JVASP-96": 216,
        "JVASP-4334": 123,
        "JVASP-4222": 11,
        "JVASP-4804": 156,
        "JVASP-329": 129,
        "JVASP-4852": 163,
        "JVASP-155": 61,
        "JVASP-4340": 51,
        "JVASP-4343": 15,
        "JVASP-4361": 14,
        "JVASP-137": 72,
        "JVASP-4879": 166,
        "JVASP-4885": 19,
        "JVASP-4894": 31,
        "JVASP-4216": 12,
        "JVASP-4918": 121,
        "JVASP-4948": 25,
        "JVASP-4957": 221,
        "JVASP-4960": 65,
        "JVASP-5059": 189,
        "JVASP-5071": 66,
        "JVASP-5074": 150,
        "JVASP-5086": 74,
        "JVASP-4388": 64,
        "JVASP-4391": 136,
        "JVASP-5155": 127,
        "JVASP-5185": 43,
        "JVASP-5197": 59,
        "JVASP-5212": 29,
        "JVASP-164": 176,
        "JVASP-5224": 137,
        "JVASP-5227": 148,
        "JVASP-4234": 193,
        "JVASP-5257": 7,
        "JVASP-5266": 140,
        "JVASP-4397": 33,
        "JVASP-5317": 8,
        "JVASP-5332": 13,
        "JVASP-5353": 16,
        "JVASP-5371": 4,
        "JVASP-5407": 229,
        "JVASP-5416": 147,
        "JVASP-5509": 52,
        "JVASP-5536": 70,
        "JVASP-5560": 71,
        "JVASP-5680": 28,
        "JVASP-5740": 6,
        "JVASP-5839": 162,
        "JVASP-5863": 79,
        "JVASP-110": 99,
        "JVASP-579": 38,
        "JVASP-4501": 5,
        "JVASP-91": 227,
        "JVASP-41": 154,
        "JVASP-4516": 96,
        "JVASP-4564": 82,
        "JVASP-4645": 130,
        "JVASP-152": 55,
        "JVASP-4792": 88,
        "JVASP-5041": 107,
        "JVASP-5425": 87,
        "JVASP-5464": 115,
        "JVASP-5650": 157,
        "JVASP-4450": 36,
        "JVASP-22520": 152,
        "JVASP-22523": 205,
        "JVASP-11998": 53,
        "JVASP-22528": 58,
        "JVASP-22533": 41,
        "JVASP-12091": 32,
        "JVASP-22541": 215,
        "JVASP-22543": 97,
        "JVASP-22549": 9,
        "JVASP-12103": 138,
        "JVASP-22575": 40,
        "JVASP-22602": 180,
        "JVASP-22611": 182,
        "JVASP-12022": 85,
        "JVASP-22637": 111,
        "JVASP-12139": 10,
        "JVASP-22709": 92,
        "JVASP-12060": 20,
        "JVASP-12064": 1,
        "JVASP-12194": 185,
        "JVASP-13885": 128,
        "JVASP-14096": 56,
        "JVASP-14020": 122,
        "JVASP-13904": 54,
        "JVASP-13783": 135,
        "JVASP-14213": 26,
        "JVASP-14034": 23,
        "JVASP-14158": 113,
        "JVASP-14256": 21,
        "JVASP-32150": 217,
        "JVASP-28392": 224,
        "JVASP-32180": 146,
        "JVASP-32197": 57,
        "JVASP-31813": 67,
        "JVASP-31819": 219,
        "JVASP-29262": 39,
        "JVASP-29281": 102,
        "JVASP-31825": 226,
        "JVASP-29425": 143,
        "JVASP-29441": 42,
        "JVASP-29520": 73,
        "JVASP-29526": 18,
        "JVASP-29555": 149,
        "JVASP-29597": 151,
        "JVASP-29705": 174,
        "JVASP-31921": 199,
        "JVASP-33147": 161,
        "JVASP-33220": 46,
        "JVASP-33344": 114,
        "JVASP-30263": 44,
        "JVASP-33832": 155,
        "JVASP-30458": 69,
        "JVASP-30461": 144,
        "JVASP-30518": 132,
        "JVASP-36540": 198,
        "JVASP-36548": 220,
        "JVASP-36568": 3,
        "JVASP-36573": 145,
        "JVASP-35061": 131,
        "JVASP-35137": 125,
        "JVASP-35222": 204,
        "JVASP-35344": 109,
        "JVASP-42053": 173,
        "JVASP-40216": 91,
        "JVASP-38594": 165,
        "JVASP-37330": 86,
        "JVASP-37573": 84,
        "JVASP-36714": 47,
        "JVASP-36754": 98,
        "JVASP-59313": 230,
        "JVASP-46893": 95,
        "JVASP-46896": 76,
        "JVASP-45779": 30,
        "JVASP-45831": 158,
        "JVASP-46446": 35,
        "JVASP-44393": 159,
        "JVASP-44773": 22,
        "JVASP-47741": 78,
        "JVASP-47811": 181,
        "JVASP-48055": 94,
        "JVASP-48916": 34,
        "JVASP-49907": 190,
        "JVASP-50342": 223,
        "JVASP-50360": 68,
        "JVASP-50431": 37,
        "JVASP-52902": 142,
        "JVASP-52377": 24,
        "JVASP-50791": 214,
        "JVASP-54512": 108,
        "JVASP-56567": 213,
        "JVASP-54867": 126,
        "JVASP-55180": 81,
        "JVASP-57780": 212,
        "JVASP-57807": 118,
        "JVASP-57816": 50,
        "JVASP-57100": 197,
        "JVASP-57138": 116,
        "JVASP-58233": 124,
        "JVASP-59682": 200,
        "JVASP-20180": 206,
        "JVASP-21448": 203,
        "JVASP-40468": 90,
        "JVASP-42914": 17,
        "JVASP-21594": 100,
        "JVASP-21706": 188,
        "JVASP-22783": 218,
        "JVASP-24006": 202,
        "JVASP-30879": 83,
        "JVASP-31186": 49,
        "JVASP-21031": 110,
        "JVASP-21116": 192,
        "JVASP-25245": 134,
        "JVASP-7066": 169,
        "JVASP-13017": 112,
        "JVASP-58953": 105,
        "JVASP-8682": 183,
        "JVASP-9902": 80,
        "JVASP-34882": 208,
        "JVASP-34330": 179,
        "JVASP-34502": 48,
        "JVASP-62982": 178,
    }
    """
 mem=[]
 spgs=[]
 pos=[]
 info=defaultdict()
 for i in d:
     a=Atoms.from_dict(i['atoms'])
     spg=Spacegroup3D(a).space_group_number
     if spg not in spgs:
        spgs.append(spg)
        pos.append(i['jid'])
        info[i['jid']]=spg
        print (info)
        if len(spgs)==230:
            break
 """
    for i, j in few_spgs.items():
                name='POSCAR-'+str(i)
                a = Atoms.from_poscar(name)
                spg = Spacegroup3D(a).space_group_number
                assert j == spg

                lattice_mat = a.lattice_mat
                kp = Kpoints3D().automatic_length_mesh(
                    lattice_mat=lattice_mat, length=40
                )
                sym = kp.high_symm_path(a)._path
                # print (ii['jid'],a)
                # TODO: following line failing for JVASP-4222
                # x, y = kp.interpolated_points(a)
                break


def test_highsym():
    kpts = HighSymmetryKpoint3DFactory().cubic()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().fcc()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().bcc()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().tet()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().bctet1(3, 2)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().bctet2(3, 2)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orc()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcf1(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcf2(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcf3(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orci(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().orcc(1, 2, 3)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().hex()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().rhl1(47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().rhl2(47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mcl(2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc1(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc2(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc3(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc4(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().mclc5(1, 2, 3, 47)._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().tria()._kpoints["\\Gamma"]
    assert kpts[0] == 0
    kpts = HighSymmetryKpoint3DFactory().trib()._kpoints["\\Gamma"]
    assert kpts[0] == 0


test_extra_spgs()
# test_highsym()
# test_kp()
