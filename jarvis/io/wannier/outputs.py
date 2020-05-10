"""
Class for reading wannier outouts such as
wannier90.wout and wannier90_hr.dat
"""

import sys
import matplotlib
import os
import json
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import math
import time
import copy as copy
from jarvis.core.kpoints import Kpoints3D
from jarvis.core.kpoints import generate_kgrid
from jarvis.io.vasp.inputs import Poscar
from jarvis.io.vasp.outputs import Vasprun


def get_projectors_for_formula(semi_core_states=None, formula_dict={"Bi": 2, "Se": 3}):
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
    soc=True, ncells = 1, supercell = [1,1,6],surfaceonly = False
):
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


        if ncells > 1 and surfaceonly == False:

            inds_super = []
            for n in range(ncells):
                for i in inds:
                    inds_super.append(i+ n * nwan)

            return inds_super

        elif ncells > 1 and surfaceonly == True:
            inds_super = []
            for n in [0, ncells-1]:
                for i in inds:
                    inds_super.append(i+ n * nwan)

            return inds_super



    return inds


class WannierHam(object):
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

    @classmethod
    def from_dict(self, info):
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
        f = open(self.filename, "r")
        lines = f.read().splitlines()
        allines = f.readlines()
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
            self.H_val[c] = (float(lines[i].split()[5]) + 1j * float(
                lines[i].split()[6]
            )) / float(self.sym_r[rnum])
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

        return [
            nxyz[0] - self.ind[0][0],
            nxyz[1] - self.ind[1][0],
            nxyz[2] - self.ind[2][0],
        ]

    def solve_ham(self, k=[0, 0, 0], proj=None):

        # print ('solve', self.nwan, self.R.shape, self.HR.shape)

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
        eigs = []
        for i in kpath:
            # print (i)
            val, vect, p = self.solve_ham(k=i, proj=proj)
            eigs.append(val - efermi)
        return np.array(eigs)

    def get_bandstructure_plot(
        self, atoms=None, efermi=0.0, filename="bs.png", yrange=[-4, 4]
    ):
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
            "eigs.shape,eigs_vrun.shape", eigs_wan.shape, eigs_vrun.shape, nbands, nwann
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
        if plot == True:
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

    # def generate_kgrid(self, grid):
    #
    #    t = []
    #    for i in range(grid[0]):
    #        for j in range(grid[1]):
    #            for k in range(grid[2]):
    #                t.append([float(i)/(float(grid[0])) , float(j)/(float(grid[1])), float(k)/(float(grid[2]))])
    #    return t

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
            xrange[0], xrange[1] + 1e-5, (xrange[1] - xrange[0]) / float(nenergy)
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
            if not proj is None:
                pdos[i] = np.sum(np.exp(arg) * pvals)

        de = energies[1] - energies[0]
        dos = dos / sig / (2.0 * np.pi) ** 0.5 / float(nk)
        if not proj is None:
            pdos = pdos / sig / (2.0 * np.pi) ** 0.5 / float(nk)
        print("np.sum(dos) ", np.sum(dos * de))
        if not proj is None:
            print("np.sum(pdos) ", np.sum(pdos * de))

        return energies, dos, pdos

    def generate_supercell(self, supercell=[2, 2, 2], index=[0, 0, 1], sparse=False):

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
            t = ss[0] * supercell[1] * supercell[2] + ss[1] * supercell[2] + ss[2]
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

                        #                        print 'rs', rold, subcell, 'new', cellnew, subnew
                        if tuple(cellnew) not in RH_new:
                            if sparse:
                                RH_new[tuple(cellnew)] = [
                                    cellnew,
                                    sps.lil_matrix((NWAN, NWAN), dtype=complex),
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

                        #                        print r1,r2,h_temp.shape, RH_new[tuple(cellnew)][1][r1,r2].shape
                        for c1, c2 in enumerate(r1):
                            for d1, d2 in enumerate(r2):
                                RH_new[tuple(cellnew)][1][c2, d2] += h_temp[c1, d1]

        t2 = time.time()

        nr = len(RH_new)
        hbig = WannierHam(nr=nr)
        # hbig = wan_ham()
        # if sparse:
        #    hbig.sparse = True

        nwan = NWAN
        hbig.nwan = NWAN
        # hbig.nr = rn

        R = np.zeros((nr, 3), dtype=float)
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


class Wannier90wout(object):
    def __init__(self, wout_path="wannier90.wout"):
        self.wout = wout_path

    def give_wannier_centers(self):
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
                    tmp = [float(j) for j in i.split("(")[1].split(")")[0].split(",")]
                    wan_cnts.append(tmp)
        return wan_cnts

"""
if __name__ == "__main__":
    hr = "../../tests/testfiles/io/wannier/wannier90_hr.dat"
    w = WannierHam(filename=hr)  # get_bandstructure_plot(atoms=p)
    big = w.generate_supercell()

    pp = get_projectors_for_formula()
    x = get_orbitals()
    print(w.nwan, big.nwan)
    sys.exit()

    # info = w.to_dict()
    # print (info['nwan'])
    # ww = WannierHam()
    # dd = WannierHam.from_dict(info)
    # print (dd.to_dict()['nwan'])
    energies, dos, pdos = w.dos([5, 5, 5])
    print(round(dos[75], 3))
    # w.compare_dft_wann(vasprun_path=run)
    import sys

    vrun = Vasprun(run)
    eigs_wan = w.band_structure_eigs(
        kpath=vrun.kpoints._kpoints, efermi=vrun.efermi
    )  # .T
    eigs_vrun = vrun.eigenvalues[0][:, :, 0]  # .T
    nbands = eigs_vrun.shape[1]
    nwann = eigs_wan.shape[1]
    fermi = vrun.efermi
    print("eigs.shape,eigs_vrun.shape", eigs_wan.shape, eigs_vrun.shape, nbands, nwann)
    min_arr = []
    erange = [-2 - fermi, 2 + fermi]
    for k in range(len(eigs_vrun)):
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

    wout = "/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/wannier90.wout"
    centers = Wannier90wout(wout_path=wout).give_wannier_centers()
    # print (centers)

    # val,vect,p = WannierHam(filename=hr).solve_ham()
    # print ('val=',val)
    # print ('vect=',vect)
    # print ('p=',p)

    p = Poscar.from_file(
        "/rk2/knc6/Chern3D/JVASP-1067_mp-541837_PBEBO/MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/POSCAR"
    ).atoms
    print(p)
    WannierHam(filename=hr).get_bandstructure_plot(atoms=p)

    # R,H, HR = WannierHam(filename=hr).read_ham()
    # print ('R',R)
    # print ()
    # print ()
    # print ('HR',HR)
"""
