"""
Class for reading wannier outouts.

Such as wannier90.wout and wannier90_hr.dat
"""

import os
import json
import matplotlib.pyplot as plt
import cmath
import numpy as np
import math
import time
from jarvis.core.kpoints import Kpoints3D
from jarvis.io.vasp.outputs import Vasprun
import scipy
import scipy.optimize as opt
import copy as copy


def get_projectors_for_formula(
    semi_core_states=None, formula_dict={"Bi": 2, "Se": 3}
):
    """Get semi core states from formula dict."""
    if semi_core_states is None:
        path_semi_core = str(
            os.path.join(os.path.dirname(__file__), "default_semicore.json")
        )
        f = open(path_semi_core, "r")
        semi_core_states = json.load(f)
        f.close()
    arr = []
    for i, j in formula_dict.items():
        dat = semi_core_states[i]
        arr.append([i, j, [str(k) for k in dat[2].split(",")]])
    return arr


def get_orbitals(
    projection_info=[["Bi", 2, ["s", "p"]], ["Se", 3, ["s", "p"]]],
    desired_orbitals=[["Bi", "p"]],
    soc=True,
    ncells=1,
    supercell=[1, 1, 6],
    surfaceonly=False,
):
    """Get spdf orbitals."""
    # projection_info example for Bi2Se3 with s and p orbital projections
    # [["Bi", 2, ["s","p"]], ["Se", 3, ["s","p"]]]

    # orbitals wanted example
    # [["Bi"]]   all Bi orbitals

    # [["Bi", "p"]] all Bi p orbitals

    # [["Bi", "px"], ["Bi" ,"py"]] all Bi px, py oribitals

    # [["Bi", "s"], ["Se", "s"]]  Bi s and Se s

    # so for spin-orbit

    c = 0

    projection_dict = {}

    # print "get_orbs ", so
    # print projection_info

    for proj in projection_info:

        atom = proj[0]
        natom = proj[1]
        orbitals = proj[2]

        for n in range(natom):
            for o in orbitals:
                if o == "s":
                    if (atom, "s") not in projection_dict:
                        projection_dict[(atom, "s")] = []
                    projection_dict[(atom, "s")].append(c)
                    c += 1

                elif o == "p":
                    if (atom, "p") not in projection_dict:
                        projection_dict[(atom, "p")] = []
                        projection_dict[(atom, "pz")] = []
                        projection_dict[(atom, "py")] = []
                        projection_dict[(atom, "px")] = []

                    projection_dict[(atom, "p")].append(c)
                    projection_dict[(atom, "pz")].append(c)
                    c += 1

                    projection_dict[(atom, "p")].append(c)
                    projection_dict[(atom, "px")].append(c)
                    c += 1

                    projection_dict[(atom, "p")].append(c)
                    projection_dict[(atom, "py")].append(c)
                    c += 1

                elif o == "d":
                    if (atom, "p") not in projection_dict:
                        projection_dict[(atom, "d")] = []
                        projection_dict[(atom, "dz2")] = []
                        projection_dict[(atom, "dxz")] = []
                        projection_dict[(atom, "dyz")] = []
                        projection_dict[(atom, "dx2y2")] = []
                        projection_dict[(atom, "dxy")] = []

                        projection_dict[(atom, "d")].append(c)
                        projection_dict[(atom, "dz2")].append(c)
                        c += 1

                        projection_dict[(atom, "d")].append(c)
                        projection_dict[(atom, "dxz")].append(c)
                        c += 1

                        projection_dict[(atom, "d")].append(c)
                        projection_dict[(atom, "dyz")].append(c)
                        c += 1

                        projection_dict[(atom, "d")].append(c)
                        projection_dict[(atom, "dx2y2")].append(c)
                        c += 1

                        projection_dict[(atom, "d")].append(c)
                        projection_dict[(atom, "dxy")].append(c)
                        c += 1

    nwan = c
    if soc:
        for (atom, orb) in projection_dict.keys():
            new_ind = []
            for i in projection_dict[(atom, orb)]:
                new_ind.append(i + nwan)
            projection_dict[(atom, orb)] += new_ind
        nwan = nwan * 2
    print("nwan = ", nwan)

    inds = []
    for d in desired_orbitals:
        if len(d) == 1:
            for orb in ["s", "p", "d"]:
                if (d[0], orb) in projection_dict:
                    new_orbs = projection_dict[(d[0], orb)]
                    inds += new_orbs
        else:
            new_orbs = projection_dict[tuple(d)]
            inds += new_orbs

        if ncells > 1 and not surfaceonly:

            inds_super = []
            for n in range(ncells):
                for i in inds:
                    inds_super.append(i + n * nwan)

            return inds_super

        elif ncells > 1 and surfaceonly:
            inds_super = []
            for n in [0, ncells - 1]:
                for i in inds:
                    inds_super.append(i + n * nwan)

            return inds_super

    return inds


class WannierHam(object):
    """Construct WannierHamltonian object."""

    def __init__(
        self,
        filename="wannier90_hr.dat",
        nwan=None,
        nr=None,
        sym_r=None,
        H_int=None,
        H_val=None,
        H=None,
        HR=None,
    ):
        """Initialize the class."""
        self.filename = filename
        self.nr = nr
        self.nwan = nwan
        self.sym_r = sym_r
        self.H_int = H_int
        self.H_val = H_val
        self.H = H
        self.HR = HR

        if self.nr is None:
            self.read_ham()

    def to_dict(self):
        """Convert to dictionary."""
        info = {}
        info["nr"] = self.nr
        info["filename"] = self.filename
        info["nwan"] = self.nwan
        info["sym_r"] = self.sym_r
        info["H_int"] = self.H_int
        info["H_val"] = self.H_val
        info["H"] = self.H
        info["HR"] = self.HR
        return info

    def find_nodes(
        self,
        num_occ=8,
        origin=[-0.5, -0.5, -0.5],
        k1=[1, 0, 0],
        k2=[0, 1, 0],
        k3=[0, 0, 1],
        nk1=10,
        nk2=10,
        nk3=10,
        sig=[0.01, 0.02, 0.05, 0.1, 0.2, 0.3],
        thresh=-10,
        node_tol=0.001,
        use_min=True,
    ):
        """Find nodes, Dirac, Wyel points."""
        K = np.zeros((nk1 * 4, nk2 * 4, nk3 * 4, 3), dtype=float)
        k1 = np.array(k1, dtype=float)
        k2 = np.array(k2, dtype=float)
        k3 = np.array(k3, dtype=float)

        for c1 in range(4 * nk1):
            for c2 in range(4 * nk2):
                for c3 in range(4 * nk3):
                    K[c1, c2, c3, :] = (
                        origin
                        + k1 * float(c1) / float(nk1 * 4)
                        + k2 * float(c2) / float(nk2 * 4)
                        + k3 * float(c3) / float(nk3 * 4)
                    )

                    #        VALS = np.zeros((nk1,nk2, ham.nwan))

        IMAGE = np.zeros((nk1 * 4, nk2 * 4, nk3 * 4, len(sig)))
        DIRECTGAP = np.zeros((nk1 * 4, nk2 * 4, nk3 * 4))
        DIRECTGAP[:, :, :] = 100.0

        VAL = np.ones((nk1 * 4, nk2 * 4, nk3 * 4)) * 100.0

        #        sig_0 = 0.2

        # initial pass

        for c1 in range(0, 4 * nk1, 4):
            print("c1 ", c1)
            for c2 in range(0, 4 * nk2, 4):
                for c3 in range(0, 4 * nk3, 4):
                    k = K[c1, c2, c3, :]
                    val, vect, p = self.solve_ham(k, proj=None)

                    d = val[num_occ] - val[num_occ - 1]
                    DIRECTGAP[c1, c2, c3] = d

                    # IMAGE[c1,c2,c3,-1] = np.exp( -(d)**2 / sig[-1]**2)

        if thresh == -10:
            thresh = max(np.min(DIRECTGAP) * 1.5, 0.03)
            print("set thresh ", thresh)

            # second pass

            #        maxval = np.max(IMAGE[:,:,:,-1])
            #        thresh = maxval * 0.1
            #        print("thresh ", maxval, thresh)

        def fun_to_min(k):
            val, vect, p = self.solve_ham(k, proj=None)
            return val[num_occ] - val[num_occ - 1]

        num_thresh = 0

        min_gap = 1000000000.0
        min_cond = 1000000000.0
        max_val = -1000000000.0

        dk1 = float(1) / float(nk1 * 4) / 2.0
        dk2 = float(1) / float(nk2 * 4) / 2.0
        dk3 = float(1) / float(nk3 * 4) / 2.0

        weyl_points = []
        dirac_points = []
        higher_order_points = []

        for c1 in range(0, 4 * nk1, 4):
            print("c1 ", c1)
            for c2 in range(0, 4 * nk2, 4):
                for c3 in range(0, 4 * nk3, 4):
                    for c1a in range(4):
                        for c2a in range(4):
                            for c3a in range(4):

                                if c1a > 0:
                                    c1p = (c1 + 4) % (4 * nk1)
                                else:
                                    c1p = c1
                                if c2a > 0:
                                    c2p = (c2 + 4) % (4 * nk2)
                                else:
                                    c2p = c2
                                if c3a > 0:
                                    c3p = (c3 + 4) % (4 * nk3)
                                else:
                                    c3p = c3

                                if (
                                    DIRECTGAP[c1, c2, c3] < thresh
                                    or DIRECTGAP[c1, c2, c3p] < thresh
                                    or DIRECTGAP[c1, c2p, c3] < thresh
                                    or DIRECTGAP[c1, c2p, c3p] < thresh
                                    or DIRECTGAP[c1p, c2, c3] < thresh
                                    or DIRECTGAP[c1p, c2, c3p] < thresh
                                    or DIRECTGAP[c1p, c2p, c3] < thresh
                                    or DIRECTGAP[c1p, c2p, c3p] < thresh
                                ):
                                    num_thresh += 1
                                    k = K[c1 + c1a, c2 + c2a, c3 + c3a, :]

                                    if use_min:
                                        B = opt.Bounds(
                                            k - [dk1, dk2, dk3],
                                            k + [dk1, dk2, dk3],
                                            keep_feasible=True,
                                        )
                                        val, vect, p = self.solve_ham(
                                            k, proj=None
                                        )
                                        if (
                                            val[num_occ] - val[num_occ - 1]
                                            < 0.15
                                        ):
                                            res = opt.minimize(
                                                fun_to_min,
                                                k,
                                                method="Powell",
                                                bounds=B,
                                                tol=0.5e-3,
                                                options={"maxfev": 5},
                                            )
                                            if res.fun < 0.005:
                                                res = opt.minimize(
                                                    fun_to_min,
                                                    res.x,
                                                    method="Powell",
                                                    bounds=B,
                                                    tol=1e-5,
                                                    options={"maxfev": 30},
                                                )
                                            val, vect, p = self.solve_ham(
                                                res.x, proj=None
                                            )

                                            k = res.x
                                    else:
                                        val, vect, p = self.solve_ham(
                                            k, proj=None
                                        )

                                    d = val[num_occ] - val[num_occ - 1]
                                    if d < node_tol:
                                        # degen = 0
                                        gap_mean = (
                                            val[num_occ] + val[num_occ - 1]
                                        ) / 2.0
                                        count = np.sum(
                                            np.abs(val - gap_mean)
                                            < node_tol * 2
                                        )
                                        if count == 2:
                                            weyl_points.append(copy.copy(k))
                                            print("weyl point ", k, d, count)
                                        elif count == 4 or count == 4:
                                            dirac_points.append(copy.copy(k))
                                            print("dirac point ", k, d, count)
                                        else:
                                            higher_order_points.append(
                                                copy.copy(k)
                                            )
                                            print(
                                                "nontrivial point ",
                                                k,
                                                d,
                                                count,
                                            )

                                    DIRECTGAP[c1 + c1a, c2 + c2a, c3 + c3a] = d
                                    VAL[c1 + c1a, c2 + c2a, c3 + c3a] = 0.5 * (
                                        val[num_occ] + val[num_occ - 1]
                                    )
                                    min_gap = min(d, min_gap)
                                    min_cond = min(min_cond, val[num_occ])
                                    max_val = max(max_val, val[num_occ - 1])

        for (cs, s) in enumerate(sig):
            IMAGE[:, :, :, cs] = np.exp(-((DIRECTGAP) ** 2) / s ** 2)

        print(
            "num thresh ", num_thresh, " out of ", nk1 * nk2 * nk3 * 4 * 4 * 4
        )
        print(
            "min direct gap ", min_gap, " indirect gap  ", min_cond - max_val
        )

        print("weyl ", weyl_points)
        print("dirac ", dirac_points)
        print("higher order ", higher_order_points)

        return (
            IMAGE,
            DIRECTGAP,
            VAL,
            [min_gap, min_cond - max_val],
            weyl_points,
            dirac_points,
            higher_order_points,
        )

    def chern_number_simple(
        self,
        nocc=8,
        k1=[0, 0, 0],
        k2=[0, 0, 0],
        nk1=20,
        nk2=20,
        Kmat=None,
        usemod=True,
    ):
        """Calculate Chern number."""
        if Kmat is None:

            K = np.zeros((nk1, nk2, 3), dtype=float)

            k1 = np.array(k1, dtype=float)
            k2 = np.array(k2, dtype=float)

            #        print( 'K')

            for c1 in range(nk1):
                for c2 in range(nk2):
                    K[c1, c2, :] = k1 * float(c1) / float(nk1) + k2 * float(
                        c2
                    ) / float(nk2)

        else:
            K = Kmat

        if usemod:

            lim1 = nk1
            lim2 = nk2

        else:

            lim1 = nk1 + 1
            lim2 = nk2 + 1

        VECT = np.zeros((lim1, lim2, self.nwan, self.nwan), dtype=complex)
        gap_min = 100000000000.0
        val_max = -10000000000000.0
        cond_min = 10000000000.0

        # cond_min_p1 = 10000000000.0

        # gap_min2 = 100000000000.0

        # nwan = self.nwan

        for c1 in range(lim1):
            for c2 in range(lim2):
                k = K[c1, c2, :]
                val, vect, p = self.solve_ham(k, proj=None)
                VECT[c1, c2, :, :] = vect
                if (val[nocc] - val[nocc - 1]) < gap_min:
                    gap_min = val[nocc] - val[nocc - 1]

                if val[nocc - 1] > val_max:
                    val_max = val[nocc - 1]

                if val[nocc] < cond_min:
                    cond_min = val[nocc]

        print(
            "minimum_direct_gap", gap_min, "indirect_gap", cond_min - val_max
        )
        direct_gap = gap_min
        indirect_gap = cond_min - val_max

        Chern = 0.0
        for c1 in range(nk1):
            for c2 in range(nk2):

                if usemod:
                    c1p = (c1 + 1) % nk1
                    c2p = (c2 + 1) % nk2
                else:
                    c1p = c1 + 1
                    c2p = c2 + 1

                klist = [[c1, c2], [c1p, c2], [c1p, c2p], [c1, c2p], [c1, c2]]

                phi = 0.0
                for i in range(4):
                    M = np.dot(
                        VECT[klist[i][0], klist[i][1], :, 0:nocc].T.conj(),
                        VECT[klist[i + 1][0], klist[i + 1][1], :, 0:nocc],
                    )
                    phi = phi - (cmath.log(np.linalg.det(M))).imag

                if phi > np.pi:
                    phi = phi - 2 * np.pi
                elif phi < -np.pi:
                    phi = phi + 2 * np.pi

                if abs(phi) > 0.75:
                    print(
                        "WARNING in chern calculation, phi",
                        c1,
                        c2,
                        phi,
                        " try more k-points",
                    )

                Chern += phi / (2.0 * np.pi)

        if usemod:
            return Chern, direct_gap, indirect_gap
        else:
            return Chern, direct_gap, val_max, cond_min

    @classmethod
    def from_dict(self, info):
        """Load from dictionary."""
        w = WannierHam(
            nr=info["nr"],
            nwan=info["nwan"],
            sym_r=info["sym_r"],
            H_int=info["H_int"],
            H_val=info["H_val"],
            H=info["H"],
            HR=info["HR"],
        )
        return w

    def read_ham(self):
        """Read _hr.dat file.."""
        f = open(self.filename, "r")
        lines = f.read().splitlines()
        # allines = f.readlines()
        f.close()
        self.nwan = int(lines[1])
        self.nr = int(lines[2])
        if self.nr % 15 == 0:
            lines_r = int(math.floor(self.nr // 15))
        else:
            lines_r = int(math.floor(self.nr // 15) + 1)

        # print ('self.nwan,self.nr',self.nwan,self.nr,lines_r)
        self.sym_r = np.zeros(self.nr, dtype=float)
        # load sym ops
        for i in range(lines_r):

            num = 3 + i
            start = i * 15
            end = (i + 1) * 15
            if end > self.nr:
                end = self.nr
            self.sym_r[start:end] = [float(j) for j in lines[num].split()]
        # print (self.sym_r)

        tot = self.nwan ** 2 * self.nr
        self.H_int = np.zeros((tot, 5), dtype=int)
        self.H_val = np.zeros(tot, dtype=complex)

        c = 0
        rnum = 0
        for i in range(lines_r + 3, lines_r + 3 + tot):

            rnum = (c) // self.nwan ** 2
            # print ('lines[i].split()[0:5]',c,self.nwan**2)#,self.sym_r[rnum])
            self.H_int[c, :] = [int(j) for j in lines[i].split()[0:5]]
            # print ('lines[i][5]',lines[i])
            self.H_val[c] = (
                float(lines[i].split()[5]) + 1j * float(lines[i].split()[6])
            ) / float(self.sym_r[rnum])
            c += 1

        # print (self.H_int[0,:])
        # print (self.H_val[0])

        # print (self.H_int[-1,:])
        # print (self.H_val[-1])

        # print ('loaded ', self.filename)
        # print ('nwan: ', self.nwan)
        # print ('nr:   ', self.nr)
        # print ()

        # reshape

        nx1 = np.min(self.H_int[:, 0])
        nx2 = np.max(self.H_int[:, 0])
        ny1 = np.min(self.H_int[:, 1])
        ny2 = np.max(self.H_int[:, 1])
        nz1 = np.min(self.H_int[:, 2])
        nz2 = np.max(self.H_int[:, 2])

        self.ind = [[nx1, nx2], [ny1, ny2], [nz1, nz2]]

        ix = nx2 - nx1 + 1
        iy = ny2 - ny1 + 1
        iz = nz2 - nz1 + 1

        print("H size", ix, iy, iz, self.nwan, self.nwan)

        self.H = np.zeros((ix, iy, iz, self.nwan, self.nwan), dtype=complex)

        self.ind_dict = {}
        for i in range(self.H_val.shape[0]):
            ind = self.get_ind(self.H_int[i, 0:3])

            self.ind_dict[(ind[0], ind[1], ind[2])] = i

            nw1 = self.H_int[i, 3]
            nw2 = self.H_int[i, 4]

            self.H[ind[0], ind[1], ind[2], nw1 - 1, nw2 - 1] = self.H_val[i]

        # print ('done reshaping1')
        nr = ix * iy * iz
        self.R = np.zeros((nr, 3), dtype=float)
        self.HR = np.zeros((nr, self.nwan ** 2), dtype=complex)

        c = 0
        for x in range(nx1, nx2 + 1):
            for y in range(ny1, ny2 + 1):
                for z in range(nz1, nz2 + 1):
                    #                    ind = self.ind_dict[(x,y,z)]
                    ind = self.get_ind([x, y, z])
                    self.R[c, :] = [x, y, z]
                    self.HR[c, :] = np.reshape(
                        self.H[ind[0], ind[1], ind[2], :, :], self.nwan ** 2
                    )
                    c += 1

        if c != nr:
            ValueError("c is not equal to r", c, nr)

        return self.R, self.H, self.HR

    def get_ind(self, nxyz):
        """Get index."""
        return [
            nxyz[0] - self.ind[0][0],
            nxyz[1] - self.ind[1][0],
            nxyz[2] - self.ind[2][0],
        ]

    def solve_ham(self, k=[0, 0, 0], proj=None):
        """Solve Wannier Hamiltonian at a k-point."""
        nr = self.R.shape[0]
        # print ('nr==',nr,self.nr)
        hk = np.zeros((self.nwan, self.nwan), dtype=complex)

        kmat = np.tile(k, (nr, 1))

        exp_ikr = np.exp(1.0j * 2 * np.pi * np.sum(kmat * self.R, 1))

        temp = np.zeros(self.nwan ** 2, dtype=complex)
        for i in range(nr):
            temp += exp_ikr[i] * self.HR[i, :]

        hk = np.reshape(temp, (self.nwan, self.nwan))

        hk = (hk + hk.T.conj()) / 2.0

        val, vect = np.linalg.eigh(hk)

        if proj is not None:
            p = np.real(np.sum(vect[proj, :] * np.conj(vect[proj, :]), 0))
        else:
            p = np.ones(val.shape)

        #        print 'proj', np.sum(np.sum(p))

        return val.real, vect, p

    def band_structure_eigs(self, kpath=None, proj=None, efermi=0.0):
        """Get eigenvalues for band eigenvalues."""
        eigs = []
        for i in kpath:
            # print (i)
            val, vect, p = self.solve_ham(k=i, proj=proj)
            eigs.append(val - efermi)
        return np.array(eigs)

    def get_bandstructure_plot(
        self, atoms=None, efermi=0.0, filename="bs.png", yrange=[-4, 4]
    ):
        """Get bandstructure plot.."""
        kpoints, labels = Kpoints3D().interpolated_points(atoms)  # _path

        # print ('kpath',kpoints)

        eigs = self.band_structure_eigs(kpath=kpoints, efermi=efermi).T

        for i, ii in enumerate(eigs):
            plt.plot(ii, color="b")

        plt.ylim(yrange)

        plt.savefig(filename)

        plt.close()

    def compare_dft_wann(
        self,
        vasprun_path="",
        energy_tol=0.75,
        plot=True,
        kp_labels_points=[],
        kp_labels=[],
        filename="compare.png",
    ):
        """Compare DFT and Wannier bands to check accuracy."""
        vrun = Vasprun(vasprun_path)
        kpoints = vrun.kpoints._kpoints
        fermi = vrun.efermi
        eigs_wan = self.band_structure_eigs(
            kpath=kpoints, efermi=vrun.efermi
        )  # -fermi#[::-1]#.T
        eigs_vrun = vrun.eigenvalues[0][:, :, 0] - fermi  # [::-1]#.T
        nbands = eigs_vrun.shape[1]
        nwann = eigs_wan.shape[1]
        info = {}
        print(
            "eigs.shape,eigs_vrun.shape",
            eigs_wan.shape,
            eigs_vrun.shape,
            nbands,
            nwann,
        )
        min_arr = []
        erange = [-energy_tol, energy_tol]
        for k in range(len(kpoints)):
            for n in eigs_wan[k]:
                diff_arr = []
                if n > erange[0] and n < erange[1]:
                    for v in eigs_vrun[k]:
                        diff = abs(n - v)
                        diff_arr.append(diff)
                if diff_arr != []:
                    tmp = np.min(diff_arr)
                    min_arr.append(tmp)
        maxdiff = "na"
        if min_arr != []:
            # print ('min_arr',min_arr)
            print("MAX diff", max(min_arr))
            maxdiff = max(min_arr)
        print("maxdiff", maxdiff)
        if plot:
            for i, ii in enumerate(eigs_wan.T):
                plt.plot(ii, color="b")
            for i, ii in enumerate(eigs_vrun.T):
                plt.plot(ii, color="r")
            plt.ylim([-energy_tol, energy_tol])
            if kp_labels_points != [] and kp_labels != []:
                plt.xticks(kp_labels_points, kp_labels)

            plt.savefig(filename)
            plt.close()
        info["eigs_wan"] = list(eigs_wan.tolist())
        info["eigs_vrun"] = list(eigs_vrun.tolist())
        info["kp_labels_points"] = list(kp_labels_points)
        info["kp_labels"] = kp_labels
        info["maxdiff"] = maxdiff
        info["efermi"] = fermi
        # print (info)
        return info  # ,eigs_wan.T,eigs_vrun.T

    def dos(
        self,
        kpoints=[],
        proj=None,
        efermi=0.0,
        xrange=None,
        nenergy=100,
        sig=0.02,
        pdf="dos.pdf",
        show=True,
    ):
        """Get density of states."""
        # plt.clf()

        # kgrid = generate_kgrid(grid)
        nk = len(kpoints)
        nwan = self.nwan

        vals = np.zeros((nk, nwan), dtype=float)
        pvals = np.zeros((nk, nwan), dtype=float)

        for i, k in enumerate(kpoints):
            val, vect, p = self.solve_ham(k, proj)
            vals[i, :] = val - efermi
            pvals[i, :] = p

        # print vals
        # print "pvals"
        # print pvals
        # print "np.sum pvals ", np.sum(np.sum(pvals))

        if xrange is None:
            vmin = np.min(vals[:])
            vmax = np.max(vals[:])
            vmin2 = vmin - (vmax - vmin) * 0.05
            vmax2 = vmax + (vmax - vmin) * 0.05
            xrange = [vmin2, vmax2]
            # plt.xlim(xrange)

        energies = np.arange(
            xrange[0],
            xrange[1] + 1e-5,
            (xrange[1] - xrange[0]) / float(nenergy),
        )
        dos = np.zeros(np.size(energies))
        pdos = np.zeros(np.size(energies))

        v = vals

        condmin = np.min(v[v > 0.0])
        valmax = np.max(v[v < 0.0])

        print("DOS BAND GAP ", condmin - valmax, "    ", valmax, " ", condmin)

        c = -0.5 / sig ** 2
        for i in range(np.size(energies)):
            arg = c * (v - energies[i]) ** 2
            dos[i] = np.sum(np.exp(arg))
            if proj is not None:
                pdos[i] = np.sum(np.exp(arg) * pvals)

        de = energies[1] - energies[0]
        dos = dos / sig / (2.0 * np.pi) ** 0.5 / float(nk)
        if proj is not None:
            pdos = pdos / sig / (2.0 * np.pi) ** 0.5 / float(nk)
        print("np.sum(dos) ", np.sum(dos * de))
        if proj is not None:
            print("np.sum(pdos) ", np.sum(pdos * de))

        return energies, dos, pdos

    def generate_supercell(
        self, supercell=[2, 2, 2], index=[0, 0, 0], sparse=False
    ):
        """Generate supercell."""
        t0 = time.time()

        nw = self.nwan

        factor = np.prod(supercell)
        NWAN = factor * self.nwan

        def plus_r(rold, subcell):
            rnew = subcell + rold
            # print ('cell',rnew,supercell)
            cellnew = rnew // supercell  # this is integer division
            subnew = rnew % supercell

            return cellnew, subnew

        def subcell_index(ss):
            t = (
                ss[0] * supercell[1] * supercell[2]
                + ss[1] * supercell[2]
                + ss[2]
            )
            return range(t * nw, (t + 1) * nw)

        RH_new = {}
        h_temp = np.zeros((self.nwan, self.nwan), dtype=complex)
        subcell = np.zeros(3, dtype=int)

        t1 = time.time()

        for ii in range(self.R.shape[0]):

            rold = np.array(self.R[ii, :], dtype=int)

            for i in range(supercell[0]):
                for j in range(supercell[1]):
                    for k in range(supercell[2]):
                        subcell[:] = [i, j, k]

                        cellnew, subnew = plus_r(rold, subcell)

                        if (
                            (index[0] > 0 and cellnew[0] != 0)
                            or (index[1] > 0 and cellnew[1] != 0)
                            or (index[2] > 0 and cellnew[2] != 0)
                        ):
                            continue

                        if tuple(cellnew) not in RH_new:
                            if sparse:
                                sps = scipy.sparse.csr_matrix  # TO CHECK
                                RH_new[tuple(cellnew)] = [
                                    cellnew,
                                    sps.lil_matrix(
                                        (NWAN, NWAN), dtype=complex
                                    ),
                                ]
                            else:
                                RH_new[tuple(cellnew)] = [
                                    cellnew,
                                    np.zeros((NWAN, NWAN), dtype=complex),
                                ]

                        h_temp[:, :] = np.reshape(
                            self.HR[ii, :], (self.nwan, self.nwan)
                        )

                        r1 = subcell_index(subcell)
                        r2 = subcell_index(subnew)

                        for c1, c2 in enumerate(r1):
                            for d1, d2 in enumerate(r2):
                                RH_new[tuple(cellnew)][1][c2, d2] += h_temp[
                                    c1, d1
                                ]

        t2 = time.time()

        nr = len(RH_new)
        hbig = WannierHam(nr=nr)
        # hbig = wan_ham()
        # if sparse:
        #    hbig.sparse = True

        # nwan = NWAN
        hbig.nwan = NWAN
        # hbig.nr = rn

        # R = np.zeros((nr, 3), dtype=float)
        hbig.R = np.zeros((nr, 3), dtype=float)
        if sparse:
            # HR = sps.lil_matrix((nr, NWAN ** 2), dtype=complex)
            hbig.HR = sps.lil_matrix((nr, NWAN ** 2), dtype=complex)
        else:
            # HR = np.zeros((nr, NWAN ** 2), dtype=complex)
            hbig.HR = np.zeros((nr, NWAN ** 2), dtype=complex)

        for c, i in enumerate(RH_new):
            h = RH_new[i][1]
            r = RH_new[i][0]
            if sparse:
                # HR[c, :] = sps.lil_matrix.reshape(h, NWAN * NWAN)
                hbig.HR[c, :] = sps.lil_matrix.reshape(h, NWAN * NWAN)
            else:
                # HR[c, :] = np.reshape(h, NWAN * NWAN)
                hbig.HR[c, :] = np.reshape(h, NWAN * NWAN)

            # R[c, :] = r
            hbig.R[c, :] = r

        t3 = time.time()
        print("TIME SUPERCELL", t1 - t0, t2 - t1, t3 - t2)
        return hbig

    def fermi_surf_2d(
        self,
        fermi=7.8163,
        origin=np.array([-1.0, -1.0, 0]),
        k1=np.array([4.0, 0.0, 0.0]),
        k2=np.array([0.0, 4.0, 0.0]),
        nk1=50,
        nk2=50,
        sig=0.3,
    ):
        """Generate 2D projected Fermi-surface."""
        K = np.zeros((nk1, nk2, 3), dtype=float)
        for c1 in range(nk1):
            for c2 in range(nk2):
                K[c1, c2, :] = (
                    origin
                    + k1 * float(c1) / float(nk1)
                    + k2 * float(c2) / float(nk2)
                )

        image = np.zeros((nk1, nk2))

        for c1 in range(nk1):
            for c2 in range(nk2):
                k = K[c1, c2, :]
                val, vect, p = self.solve_ham(k=k)
                image[c1, c2] = np.sum(
                    np.exp(-((val - fermi) ** 2) / sig ** 2)
                )

        return image


class Wannier90wout(object):
    """Construct wannier90.out related object."""

    def __init__(self, wout_path="wannier90.wout"):
        """Initialize with file path."""
        self.wout = wout_path

    def give_wannier_centers(self):
        """Get wannier charge centers."""
        f = open(self.wout, "r")
        lines = f.read().splitlines()
        f.close()
        final = False
        wan_cnts = []
        for ii, i in enumerate(lines):
            if "Final State" in i:
                final = True
            if final:
                if "WF centre and spread" in i:
                    tmp = [
                        float(j)
                        for j in i.split("(")[1].split(")")[0].split(",")
                    ]
                    wan_cnts.append(tmp)
        return wan_cnts


class Wannier90eig(object):
    """Construct wannier90.eig related object."""

    def __init__(self, weig_path="wannier90.eig"):
        """Initialize with file path."""
        self.weig = weig_path

    def give_wannier_eigs(self):
        """Get wannier eigs."""
        eigs = np.loadtxt(self.weig)
        return eigs

    def neigs(self):
        """Get wannier eigs."""
        ne = int(self.give_wannier_eigs()[-1][0])
        return ne

    def nk(self):
        """Get wannier eigs."""
        nk = int(len(self.give_wannier_eigs()) / self.neigs())
        return nk
