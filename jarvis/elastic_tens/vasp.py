from pymatgen.core.structure import Structure
import numpy as np
from pymatgen.analysis.elasticity.elastic import ElasticTensor
import os


def elastic_props(outcar="", vacuum=False):
    """
    Obtain elastic tensor and calculate related properties

    Args:
        outcar: OUTCAR file path
        vacuum: whether the structure has vaccum such as 2D materials
        for vacuum structures bulk and shear mod. needs extra attenstion
        and elastic tensor are in Nm^-1 rather than GPa
    Returns:
          info: data for elastic tensor (in string and object representation), bulk, shear modulus, and phonon modes
    """
    if vacuum == True:
        contcar = Structure.from_file(str(outcar).replace("OUTCAR", "POSCAR"))
        ratio_c = 0.1 * float(
            abs(contcar.lattice.matrix[2][2])
        )  # *(10**9)*(10**-10) #N/m unit
    el_tens = "na"
    KV = "na"
    GV = "na"
    spin = "na"
    info = {}
    ratio_c = 1.0
    v = open(outcar, "r")
    lines = v.read().splitlines()
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            c11 = lines[i + 3].split()[1]
            c12 = lines[i + 3].split()[2]
            c13 = lines[i + 3].split()[3]
            c14 = lines[i + 3].split()[4]
            c15 = lines[i + 3].split()[5]
            c16 = lines[i + 3].split()[6]
            c21 = lines[i + 4].split()[1]
            c22 = lines[i + 4].split()[2]
            c23 = lines[i + 4].split()[3]
            c24 = lines[i + 4].split()[4]
            c25 = lines[i + 4].split()[5]
            c26 = lines[i + 4].split()[6]
            c31 = lines[i + 5].split()[1]
            c32 = lines[i + 5].split()[2]
            c33 = lines[i + 5].split()[3]
            c34 = lines[i + 5].split()[4]
            c35 = lines[i + 5].split()[5]
            c36 = lines[i + 5].split()[6]
            c41 = lines[i + 6].split()[1]
            c42 = lines[i + 6].split()[2]
            c43 = lines[i + 6].split()[3]
            c44 = lines[i + 6].split()[4]
            c45 = lines[i + 6].split()[5]
            c46 = lines[i + 6].split()[6]
            c51 = lines[i + 7].split()[1]
            c52 = lines[i + 7].split()[2]
            c53 = lines[i + 7].split()[3]
            c54 = lines[i + 7].split()[4]
            c55 = lines[i + 7].split()[5]
            c56 = lines[i + 7].split()[6]
            c61 = lines[i + 8].split()[1]
            c62 = lines[i + 8].split()[2]
            c63 = lines[i + 8].split()[3]
            c64 = lines[i + 8].split()[4]
            c65 = lines[i + 8].split()[5]
            c66 = lines[i + 8].split()[6]
            c11 = round(ratio_c * float(c11) / float(10), 1)
            c12 = round(ratio_c * float(c12) / float(10), 1)
            c13 = round(ratio_c * float(c13) / float(10), 1)
            c14 = round(ratio_c * float(c14) / float(10), 1)
            c15 = round(ratio_c * float(c15) / float(10), 1)
            c16 = round(ratio_c * float(c16) / float(10), 1)
            c21 = round(ratio_c * float(c21) / float(10), 1)
            c22 = round(ratio_c * float(c22) / float(10), 1)
            c23 = round(ratio_c * float(c23) / float(10), 1)
            c24 = round(ratio_c * float(c24) / float(10), 1)
            c25 = round(ratio_c * float(c25) / float(10), 1)
            c26 = round(ratio_c * float(c26) / float(10), 1)
            c31 = round(float(c31) / float(10), 1)
            c32 = round(float(c32) / float(10), 1)
            c33 = round(float(c33) / float(10), 1)
            c34 = round(float(c34) / float(10), 1)
            c35 = round(float(c35) / float(10), 1)
            c36 = round(float(c36) / float(10), 1)
            c41 = round(float(c41) / float(10), 1)
            c42 = round(float(c42) / float(10), 1)
            c43 = round(float(c43) / float(10), 1)
            c44 = round(float(c44) / float(10), 1)
            c45 = round(float(c45) / float(10), 1)
            c46 = round(float(c46) / float(10), 1)
            c51 = round(float(c51) / float(10), 1)
            c52 = round(float(c52) / float(10), 1)
            c53 = round(float(c53) / float(10), 1)
            c54 = round(float(c54) / float(10), 1)
            c55 = round(float(c55) / float(10), 1)
            c56 = round(float(c56) / float(10), 1)
            c61 = round(float(c61) / float(10), 1)
            c62 = round(float(c62) / float(10), 1)
            c63 = round(float(c63) / float(10), 1)
            c64 = round(float(c64) / float(10), 1)
            c65 = round(float(c65) / float(10), 1)
            c66 = round(float(c66) / float(10), 1)
            KV = float((c11 + c22 + c33) + 2 * (c12 + c23 + c31)) / float(9)
            GV = float(
                (c11 + c22 + c33) - (c12 + c23 + c31) + 3 * (c44 + c55 + c66)
            ) / float(15)
            KV = round(KV, 3)
            GV = round(GV, 3)
            break
    v.close()

    # Convenient string representation for storage

    el_tens = (
        str(c11)
        + str(",")
        + str(c12)
        + str(",")
        + str(c13)
        + str(",")
        + str(c14)
        + str(",")
        + str(c15)
        + str(",")
        + str(c16)
        + str(",")
        + str(c21)
        + str(",")
        + str(c22)
        + str(",")
        + str(c23)
        + str(",")
        + str(c24)
        + str(",")
        + str(c25)
        + str(",")
        + str(c26)
        + str(",")
        + str(c31)
        + str(",")
        + str(c32)
        + str(",")
        + str(c33)
        + str(",")
        + str(c34)
        + str(",")
        + str(c35)
        + str(",")
        + str(c36)
        + str(",")
        + str(c41)
        + str(",")
        + str(c42)
        + str(",")
        + str(c43)
        + str(",")
        + str(c44)
        + str(",")
        + str(c45)
        + str(",")
        + str(c46)
        + str(",")
        + str(c51)
        + str(",")
        + str(c52)
        + str(",")
        + str(c53)
        + str(",")
        + str(c54)
        + str(",")
        + str(c55)
        + str(",")
        + str(c56)
        + str(",")
        + str(c61)
        + str(",")
        + str(c62)
        + str(",")
        + str(c63)
        + str(",")
        + str(c64)
        + str(",")
        + str(c65)
        + str(",")
        + str(c66)
    )

    cij = np.empty((6, 6), dtype=float)
    elast = np.array(el_tens.split(","), dtype="float")
    count = 0
    for ii in range(6):
        for jj in range(6):
            cij[ii][jj] = elast[count]
    el_tens2 = ElasticTensor.from_voigt(cij)

    modes = []
    try:
        for i in lines:
            if "cm-1" in i and "meV" in i:

                mod = float(i.split()[7])
                if "f/i" in i:
                    mod = mod * -1
                if mod not in modes:
                    modes.append(float(mod))
    except:
        pass

    info["el_tens_str"] = el_tens
    info["el_tens_obj"] = el_tens2
    info["KV"] = KV
    info["GV"] = GV
    info["modes"] = modes

    return info


def get_et(elast_str=""):
    cij = np.empty((6, 6), dtype=float)
    elast = np.array(elast_str.split(","), dtype="float")
    count = 0
    for ii in range(6):
        for jj in range(6):
            cij[ii][jj] = elast[count]
            count = count + 1
    et = ElasticTensor.from_voigt(cij)
    return et


if __name__ == "__main__":
    out = path = str(
        os.path.join(
            os.path.dirname(__file__),
            "../vasp/examples/SiOptb88/MAIN-ELASTIC-bulk@mp_149/OUTCAR",
        )
    )
    ep = elastic_props(out)
    print(ep)
