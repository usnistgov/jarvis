"""Module to get Wannier TB format for Force-constants."""
import numpy as np
from datetime import datetime


def get_phonon_hr(
    fcmat,
    smallest_vectors,
    mass,
    multi,
    super_pos,
    p2s_map,
    s2p_map,
    num_satom,
    num_patom,
    out_file,
):
    """
    Get phonon-TB Hamiltonian similar to wannier90_hr.dat.

    Its build by FORCE_CONSTANTS
    num_satom: number of atoms in super cell
    num_patom: number of atoms in primitive cell
    fcmat: force constansts num_satom*num_satom*3*3
    Forked from WannierTools package, please cite it
    along with jarvis-tools f you use this part of the code.
    """
    # maximum dimension for hr matrix
    ndim = 51
    sdim = 10
    nrpt_max = 51 ** 3
    # hr matrix
    norbs = num_patom * 3
    hr_mat = np.zeros((ndim, ndim, ndim, norbs, norbs), dtype=np.complex128)
    hr_mat0 = np.zeros((nrpt_max, norbs, norbs), dtype=np.complex128)
    # WS points
    rpts = np.zeros((nrpt_max, 3), dtype=np.int32)
    # degeneracy
    dege = np.zeros((nrpt_max), dtype=np.int32)

    for iatom in range(num_patom):  # atoms in primitive cell
        for jatom in range(num_patom):  # atoms in primitive cell
            mass_sqrt = np.sqrt(mass[iatom] * mass[jatom])
            for katom in range(num_satom):  # atoms in supercell
                if s2p_map[katom] != p2s_map[jatom]:
                    continue
                for rr in range(np.int32(multi[katom, iatom])):
                    # find which rpt
                    rvec = (
                        smallest_vectors[katom, iatom, rr]
                        + super_pos[p2s_map[iatom]]
                        - super_pos[s2p_map[katom]]
                    )
                    for ii in range(3):
                        if abs(rvec[ii]) < 1.0e-6:
                            rvec[ii] = 0.0
                    rx = np.int32(np.round(rvec[0]))
                    ry = np.int32(np.round(rvec[1]))
                    rz = np.int32(np.round(rvec[2]))
                    idx = iatom * 3
                    idy = jatom * 3
                    rx = rx + sdim
                    ry = ry + sdim
                    rz = rz + sdim
                    hr_mat[rx, ry, rz, idx : idx + 3, idy : idy + 3] = fcmat[
                        p2s_map[iatom], katom
                    ] / (multi[katom, iatom] * mass_sqrt)

    # collect all hr at rpts with none-zero data
    irpt = 0  # count for useful rpts
    for rx in range(-sdim, sdim + 1):
        ixo = rx + sdim
        for ry in range(-sdim, sdim + 1):
            iyo = ry + sdim
            for rz in range(-sdim, sdim + 1):
                izo = rz + sdim
                if (
                    abs(hr_mat[ixo, iyo, izo, :, :]).sum() < 1.0e-9
                ):  # ommit too small
                    continue
                dege[irpt] = 1
                rpts[irpt, 0] = rx
                rpts[irpt, 1] = ry
                rpts[irpt, 2] = rz
                hr_mat0[irpt, :, :] = hr_mat[ixo, iyo, izo, :, :]
                irpt = irpt + 1
    nrpt = irpt
    dege_rpts = dege[0:nrpt]
    norbs = num_patom * 3
    with open(out_file, "w") as f:
        line = (
            " Writen on "
            + str(datetime.now())
            + "\n"
            + "          "
            + str(norbs)
            + "\n"
            + "        "
            + str(nrpt)
            + "\n"
        )
        f.write(line)
        nl = np.int32(np.ceil(nrpt / 15.0))
        for n in range(nl):
            line = "    " + "    ".join(
                [str(np.int32(i)) for i in dege_rpts[n * 15 : (n + 1) * 15]]
            )
            f.write(line)
            f.write("\n")
        for irpt in range(nrpt):
            rx = rpts[irpt, 0]
            ry = rpts[irpt, 1]
            rz = rpts[irpt, 2]
            for jatomorb in range(norbs):
                for iatomorb in range(norbs):
                    rp = hr_mat0[irpt, iatomorb, jatomorb].real
                    ip = hr_mat0[irpt, iatomorb, jatomorb].imag
                    tmp = "{:8d}{:8d}{:8d}{:8d}{:8d}{:20.10f}{:20.10f}\n"
                    line = tmp.format(
                        rx, ry, rz, iatomorb + 1, jatomorb + 1, rp, ip
                    )
                    f.write(line)
