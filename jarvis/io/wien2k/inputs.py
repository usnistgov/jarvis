"""Module to prepare input files for WIEN2k."""
from scipy import array, zeros, sqrt, dot
import os
from fractions import gcd
import numpy as np
from functools import reduce
# from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D


def get_wien_kpoints(
    atoms=None, write_file=False, line_density=1, filename="MyKpoints"
):
    """Get WIEN2k style kpoints for bandstructure calculation."""
    uniqe_lbl = []
    uniqe_k = []
    kpoints, labels = Kpoints3D().interpolated_points(
        atoms, line_density=line_density
    )
    kp = Kpoints3D().high_kpath(atoms)["kpoints"]
    for i, j in kp.items():
        uniqe_lbl.append(i)
        uniqe_k.append(j)
    legend = uniqe_lbl
    BS = uniqe_k
    data = []
    BS = array(BS, dtype=float)
    Nt = 500
    dl = zeros(len(BS) - 1)
    for i in range(len(BS) - 1):
        dr = BS[i + 1] - BS[i]
        dl[i] = sqrt(dot(dr, dr))
    dN = dl / sum(dl) * Nt
    Ni = [int(round(dn)) for dn in dN]

    Mi = zeros(len(Ni), dtype=int)
    for i in range(len(Ni)):
        tmp = np.concatenate((BS[i + 1], BS[i]))
        fracts = 1 // array(list(filter(lambda x: x > 1e-6, tmp)))
        fact = int(reduce(lambda x, y: x * y // gcd(x, y), fracts))
        if Ni[i] % fact == 0:
            Mi[i] = Ni[i]
        else:
            Mi[i] = Ni[i] * fact

    for p in range(len(Ni)):
        NAME = legend[p]
        for i in range(Ni[p]):
            kint = Mi[p] * BS[p] + (BS[p + 1] - BS[p]) * i * Mi[p] // Ni[p]
            if i > 0:
                NAME = "   "
            data.append([NAME, kint[0], kint[1], kint[2], Mi[p], 1.0])
    NAME = legend[-1]
    kint = BS[-1] * Mi[-1]
    data.append([NAME, kint[0], kint[1], kint[2], Mi[-1], 1.0])

    if write_file:
        f = open(filename, "w")
        for i in data:
            line = (
                str(i[0])
                + str(" ")
                + str(i[1])
                + str(" ")
                + str(i[2])
                + str(" ")
                + str(i[3])
                + str(" ")
                + str(i[4])
                + str(" ")
                + str(i[5])
                + str("\n")
            )
            f.write(line)
        line = str("END \n")
        f.write(line)
        f.close()
    return data


def prepare_wien_input(atoms=None, filename="case.cif"):
    """Prepare basic WIEN2k input files."""
    atoms.write_cif(filename=filename)
    cmd = str("cif2struct ") + str(filename)
    os.system(cmd)
    cmd = str(
        "init_lapw -b -red 0 -vxc 13 -ecut -7.0  -numk 100 -sp >init_w2k_out"
    )
    os.system(cmd)
    cmd = str("runsp_lapw -cc 0.0001 -ec 0.0001 -p -i 500  >run_w2k_out")
    os.system(cmd)


"""
if __name__ == "__main__":
    legend = ["Gamma", "H", "P", "N", "Gamma", "P"]
    BS = [
        [0, 0, 0],
        [0, 1, 0],
        [1 / 2.0, 1 / 2.0, 1 / 2.0],
        [1 / 2.0, 1 / 2.0, 0],
        [0, 0, 0],
        [1 / 2.0, 1 / 2.0, 1 / 2.0],
    ]

    s = Atoms.from_poscar("POSCAR")

    data = get_wien_kpoints(atoms=s, write_file=True)

    for i in data:
        print("%-10s%5d%5d%5d%5d%5.1f" % (i[0], i[1], i[2], i[3], i[4], i[5]))
"""
