#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from math import sqrt
from vasp_constant import *
from scipy.fftpack import fftfreq, fftn, ifftn

############################################################
"""
This program is based on the code written by Ren Hao <renh@upc.edu.cn>.
"""


class vaspwfc:
    """
    Class for VASP Pseudowavefunction stored in WAVECAR

    The format of VASP WAVECAR, as shown in
        http://www.andrew.cmu.edu/user/feenstra/wavetrans/
    is:
        Record-length #spin components RTAG(a value specifying the precision)
        #k-points #bands ENCUT(maximum energy for plane waves)
        LatVec-A
        LatVec-B
        LatVec-C
        Loop over spin
           Loop over k-points
              #plane waves, k vector
              Loop over bands
                 band energy, band occupation
              End loop over bands
              Loop over bands
                 Loop over plane waves
                    Plane-wave coefficient
                 End loop over plane waves
              End loop over bands
           End loop over k-points
        End loop over spin
    """

    def __init__(self, fnm="WAVECAR", lsorbit=False, lgamma=False):
        """
        Initialization.
        """

        self._fname = fnm
        self._lsoc = lsorbit
        self._lgam = lgamma

        assert not (lsorbit and lgamma), "The two settings conflict!"

        try:
            self._wfc = open(self._fname, "rb")
        except:
            raise IOError("Failed to open %s" % self._fname)

        # read the basic information
        self.readWFHeader()
        # read the band information
        self.readWFBand()

        if self._lsoc:
            assert self._nspin == 1, "NSPIN = 1 for noncollinear version WAVECAR!"

    def isSocWfc(self):
        """
        Is the WAVECAR from an SOC calculation?
        """
        return True if self._lsoc else False

    def printWF(self):
        """
        print 'KFG modification'
        print
        print 
        print 'nkpts', self._nkpts 
        print 'nbands', self._nbands
        print 'encut', self._encut 
        print 'Acell', self._Acell 
        print 'Omega', self._Omega 
        print 'Bcell', self._Bcell 
        print 'grid', self._ngrid
        print 'spin orbit', self._lsoc
#        print 'kpoints'
#        print  self._kvecs
        """
        print

    def readWFHeader(self):
        """
        Read the system information from WAVECAR, which is written in the first
        two record.

        rec1: recl, nspin, rtag
        rec2: nkpts, nbands, encut, ((cell(i,j) i=1, 3), j=1, 3)
        """

        # goto the start of the file and read the first record
        self._wfc.seek(0)
        self._recl, self._nspin, self._rtag = np.array(
            np.fromfile(self._wfc, dtype=np.float, count=3), dtype=int
        )
        self._WFPrec = self.setWFPrec()
        # the second record
        self._wfc.seek(self._recl)
        dump = np.fromfile(self._wfc, dtype=np.float, count=12)

        self._nkpts = int(dump[0])  # No. of k-points
        self._nbands = int(dump[1])  # No. of bands
        self._encut = dump[2]  # Energy cutoff
        self._Acell = dump[3:].reshape((3, 3))  # real space supercell basis
        self._Omega = np.linalg.det(self._Acell)  # real space supercell volume
        self._Bcell = np.linalg.inv(self._Acell).T  # reciprocal space supercell volume

        # Minimum FFT grid size
        Anorm = np.linalg.norm(self._Acell, axis=1)
        CUTOF = np.ceil(sqrt(self._encut / RYTOEV) / (TPI / (Anorm / AUTOA)))
        self._ngrid = np.array(2 * CUTOF + 1, dtype=int)

    def setWFPrec(self):
        """
        Set wavefunction coefficients precision:
            TAG = 45200: single precision complex, np.complex64, or complex(qs)
            TAG = 45210: double precision complex, np.complex128, or complex(q)
        """
        if self._rtag == 45200:
            return np.complex64
        elif self._rtag == 45210:
            return np.complex128
        elif self._rtag == 53300:
            raise ValueError("VASP5 WAVECAR format, not implemented yet")
        elif self._rtag == 53310:
            raise ValueError(
                "VASP5 WAVECAR format with double precision "
                + "coefficients, not implemented yet"
            )
        else:
            raise ValueError("Invalid TAG values: {}".format(self._rtag))

    def readWFBand(self, ispin=1, ikpt=1, iband=1):
        """
        Extract KS energies and Fermi occupations from WAVECAR.
        """

        self._nplws = np.zeros(self._nkpts, dtype=int)
        self._kvecs = np.zeros((self._nkpts, 3), dtype=float)
        self._bands = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)
        self._occs = np.zeros((self._nspin, self._nkpts, self._nbands), dtype=float)

        for ii in range(self._nspin):
            for jj in range(self._nkpts):
                rec = self.whereRec(ii + 1, jj + 1, 1) - 1
                self._wfc.seek(rec * self._recl)
                dump = np.fromfile(
                    self._wfc, dtype=np.float, count=4 + 3 * self._nbands
                )
                if ii == 0:
                    self._nplws[jj] = int(dump[0])
                    self._kvecs[jj] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self._bands[ii, jj, :] = dump[:, 0]
                self._occs[ii, jj, :] = dump[:, 2]

        if self._nkpts > 1:
            tmp = np.linalg.norm(
                np.dot(np.diff(self._kvecs, axis=0), self._Bcell), axis=1
            )
            self._kpath = np.concatenate(([0], np.cumsum(tmp)))
        else:
            self._kpath = None
        return self._kpath, self._bands

    def gvectors(self, ikpt=1):
        """
        Generate the G-vectors that satisfies the following relation
            (G + k)**2 / 2 < ENCUT
        """
        assert 1 <= ikpt <= self._nkpts, "Invalid kpoint index!"

        kvec = self._kvecs[ikpt - 1]
        # fx, fy, fz = [fftfreq(n) * n for n in self._ngrid]
        # fftfreq in scipy.fftpack is a little different with VASP frequencies
        fx = [
            ii if ii < self._ngrid[0] / 2 + 1 else ii - self._ngrid[0]
            for ii in range(self._ngrid[0])
        ]
        fy = [
            jj if jj < self._ngrid[1] / 2 + 1 else jj - self._ngrid[1]
            for jj in range(self._ngrid[1])
        ]
        fz = [
            kk if kk < self._ngrid[2] / 2 + 1 else kk - self._ngrid[2]
            for kk in range(self._ngrid[2])
        ]
        if self._lgam:
            # parallel gamma version of VASP WAVECAR exclude some planewave
            # components, -DwNGZHalf
            kgrid = np.array(
                [
                    (fx[ii], fy[jj], fz[kk])
                    for kk in range(self._ngrid[2])
                    for jj in range(self._ngrid[1])
                    for ii in range(self._ngrid[0])
                    if (
                        (fz[kk] > 0)
                        or (fz[kk] == 0 and fy[jj] > 0)
                        or (fz[kk] == 0 and fy[jj] == 0 and fx[ii] >= 0)
                    )
                ],
                dtype=float,
            )
        else:
            kgrid = np.array(
                [
                    (fx[ii], fy[jj], fz[kk])
                    for kk in range(self._ngrid[2])
                    for jj in range(self._ngrid[1])
                    for ii in range(self._ngrid[0])
                ],
                dtype=float,
            )

        # Kinetic_Energy = (G + k)**2 / 2
        # HSQDTM    =  hbar**2/(2*ELECTRON MASS)
        KENERGY = (
            HSQDTM
            * np.linalg.norm(
                np.dot(kgrid + kvec[np.newaxis, :], TPI * self._Bcell), axis=1
            )
            ** 2
        )
        # find Gvectors where (G + k)**2 / 2 < ENCUT
        Gvec = kgrid[np.where(KENERGY < self._encut)[0]]

        if self._lsoc:
            assert Gvec.shape[0] == self._nplws[ikpt - 1] / 2, (
                "No. of planewaves not consistent for an SOC WAVECAR! %d %d %d"
                % (Gvec.shape[0], self._nplws[ikpt - 1], np.prod(self._ngrid))
            )
        else:
            assert Gvec.shape[0] == self._nplws[ikpt - 1], (
                "No. of planewaves not consistent! %d %d %d"
                % (Gvec.shape[0], self._nplws[ikpt - 1], np.prod(self._ngrid))
            )

        return np.asarray(Gvec, dtype=int)

    def save2vesta(self, phi=None, poscar="POSCAR", prefix="wfc"):
        """
        Save the real space pseudo-wavefunction as vesta format.
        """
        nx, ny, nz = phi.shape
        try:
            pos = open(poscar, "r")
            head = ""
            for line in pos:
                if line.strip():
                    head += line
                else:
                    break
            head += "\n%5d%5d%5d\n" % (nx, ny, nz)
        except:
            raise IOError("Failed to open %s" % poscar)

        with open(prefix + "_r.vasp", "w") as out:
            out.write(head)
            nwrite = 0
            for kk in range(nz):
                for jj in range(ny):
                    for ii in range(nx):
                        nwrite += 1
                        out.write("%16.8E " % phi.real[ii, jj, kk])
                        if nwrite % 10 == 0:
                            out.write("\n")
        if not self._lgam:
            with open(prefix + "_i.vasp", "w") as out:
                out.write(head)
                nwrite = 0
                for kk in range(nz):
                    for jj in range(ny):
                        for ii in range(nx):
                            nwrite += 1
                            out.write("%16.8E " % phi.imag[ii, jj, kk])
                            if nwrite % 10 == 0:
                                out.write("\n")

    def wfc_r(self, ispin=1, ikpt=1, iband=1, gvec=None, ngrid=None, norm=False):
        """
        Obtain the pseudo-wavefunction of the specified KS states in real space
        by performing FT transform on the reciprocal space planewave
        coefficients.  The 3D FT grid size is determined by ngrid, which
        defaults to self._ngrid if not given.  Gvectors of the KS states is used
        to put 1D planewave coefficients back to 3D grid.
        """
        self.checkIndex(ispin, ikpt, iband)

        if ngrid is None:
            ngrid = self._ngrid.copy()
        else:
            ngrid = np.array(ngrid, dtype=int)
            assert ngrid.shape == (3,)
            assert np.alltrue(ngrid >= self._ngrid), (
                "Minium FT grid size: (%d, %d, %d)"
                % (self._ngrid[0], self._ngrid[1], self._ngrid[2])
            )
        if gvec is None:
            gvec = self.gvectors(ikpt)

        if self._lgam:
            phi_k = np.zeros(
                (ngrid[0], ngrid[1], ngrid[2] / 2 + 1), dtype=np.complex128
            )
        else:
            phi_k = np.zeros(ngrid, dtype=np.complex128)

        gvec %= ngrid[np.newaxis, :]

        if self._lsoc:
            wfc_spinor = []
            dump = self.readBandCoeff(ispin, ikpt, iband, norm)
            nplw = dump.shape[0] / 2

            # spinor up
            phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = dump[:nplw]
            wfc_spinor.append(ifftn(phi_k))
            # spinor down
            phi_k[:, :, :] = 0.0j
            phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = dump[nplw:]
            wfc_spinor.append(ifftn(phi_k))

            del dump
            return wfc_spinor

        else:
            phi_k[gvec[:, 0], gvec[:, 1], gvec[:, 2]] = self.readBandCoeff(
                ispin, ikpt, iband, norm
            )

            if self._lgam:
                # add some components that are excluded and perform c2r FFT
                for ii in range(ngrid[0]):
                    for jj in range(ngrid[1]):
                        fx = ii if ii < ngrid[0] / 2 + 1 else ii - ngrid[0]
                        fy = ii if ii < ngrid[1] / 2 + 1 else ii - ngrid[1]
                        if (fy > 0) or (fy == 0 and fx >= 0):
                            continue
                        phi_k[ii, jj, 0] = phi_k[-ii, -jj, 0].conjugate()
                phi_k /= np.sqrt(2.0)
                phi_k[0, 0, 0] *= np.sqrt(2.0)
                return np.fft.irfftn(phi_k, s=ngrid)
            else:
                # perform complex2complex FFT
                return ifftn(phi_k)

    def readBandCoeff(self, ispin=1, ikpt=1, iband=1, norm=False):
        """
        Read the planewave coefficients of specified KS states.
        """

        self.checkIndex(ispin, ikpt, iband)

        rec = self.whereRec(ispin, ikpt, iband)
        self._wfc.seek(rec * self._recl)

        nplw = self._nplws[ikpt - 1]
        dump = np.fromfile(self._wfc, dtype=self._WFPrec, count=nplw)

        cg = np.asarray(dump, dtype=np.complex128)
        if norm:
            cg /= np.linalg.norm(cg)
        return cg

    def whereRec(self, ispin=1, ikpt=1, iband=1):
        """
        Return the rec position for specified KS state.
        """

        self.checkIndex(ispin, ikpt, iband)

        rec = (
            2
            + (ispin - 1) * self._nkpts * (self._nbands + 1)
            + (ikpt - 1) * (self._nbands + 1)
            + iband
        )
        return rec

    def checkIndex(self, ispin, ikpt, iband):
        """
        Check if the index is valid!
        """
        assert 1 <= ispin <= self._nspin, "Invalid spin index!"
        assert 1 <= ikpt <= self._nkpts, "Invalid kpoint index!"
        assert 1 <= iband <= self._nbands, "Invalid band index!"

    def TransitionDipoleMoment(self, ks_i, ks_j, norm=False):
        """
        calculate Transition Dipole Moment between two KS states.
        TDM in momentum representation
                                             ___              
                                  i⋅h        ╲                
        <psi_a| r | psi_b> =    --------- ⋅   ╲   Cai⋅Cbi⋅Gi
                                 Eb - Ea      ╱               
                                             ╱                
                                             ‾‾‾              
                                              i               
        Note: |psi_a> and |psi_b> should be bloch function with 
              the same k vector.

        The KS states ks_i (ks_j) is specified by list of index (ispin, ikpt, iband).
        """

        ks_i = list(ks_i)
        ks_j = list(ks_j)
        assert len(ks_i) == len(ks_j) == 3, "Must be there indexes!"
        assert ks_i[1] == ks_j[1], "k-point of the two states differ!"
        self.checkIndex(*ks_i)
        self.checkIndex(*ks_j)

        # according to the above equation, G = 0 does NOT contribute to TDM.
        gvec = np.dot(self.gvectors(ikpt=ks_i[1]), self._Bcell * TPI)
        # planewave coefficients of the two states
        phi_i = self.readBandCoeff(*ks_i, norm=norm)
        phi_j = self.readBandCoeff(*ks_j, norm=norm)
        # energy differences between the two states
        dE = (
            self._bands[ks_j[0] - 1, ks_j[1] - 1, ks_j[2] - 1]
            - self._bands[ks_i[0] - 1, ks_i[1] - 1, ks_i[2] - 1]
        )

        tmp1 = phi_i.conjugate() * phi_j
        ovlap = np.sum(tmp1)
        if self._lgam:
            tmp2 = phi_i * phi_j.conjugate()
            tdm = (
                np.sum(tmp1[:, np.newaxis] * gvec, axis=0)
                - np.sum(tmp2[:, np.newaxis] * gvec, axis=0)
            ) / 2.0
        else:
            tdm = np.sum(tmp1[:, np.newaxis] * gvec, axis=0)

        tdm = 1j / (dE / (2 * RYTOEV)) * tdm * AUTOA * AUTDEBYE

        return dE, ovlap, tdm

    def inverse_participation_ratio(self, norm=True):
        """
        Calculate Inverse Paticipation Ratio (IPR) from the wavefunction. IPR is
        a measure of the localization of Kohn-Sham states. For a particular KS
        state \phi_j, it is defined as
            
                            \sum_n |\phi_j(n)|^4 
            IPR(\phi_j) = -------------------------
                          |\sum_n |\phi_j(n)|^2||^2

        where n iters over the number of grid points.
        """

        self.ipr = np.zeros((self._nspin, self._nkpts, self._nbands, 3))

        for ispin in range(self._nspin):
            for ikpt in range(self._nkpts):
                for iband in range(self._nbands):
                    phi_j = self.wfc_r(ispin + 1, ikpt + 1, iband + 1, norm=norm)
                    phi_j_abs = np.abs(phi_j)

                    # print 'Calculating IPR of #spin %4d, #kpt %4d, #band %4d' % (ispin+1, ikpt+1, iband+1)
                    self.ipr[ispin, ikpt, iband, 0] = (
                        self._kpath[ikpt] if self._kpath is None else 0
                    )
                    self.ipr[ispin, ikpt, iband, 1] = self._bands[ispin, ikpt, iband]
                    self.ipr[ispin, ikpt, iband, 2] = (
                        np.sum(phi_j_abs ** 4) / np.sum(phi_j_abs ** 2) ** 2
                    )

        np.save("ipr.npy", self.ipr)
        return self.ipr


############################################################

if __name__ == "__main__":
    # xx = vaspwfc('wavecar')
    # phi = xx.wfc_r(1, 30, 17, ngrid=(28, 28, 252))
    # xx.save2vesta(phi, poscar='POSCAR')

    # xx = vaspwfc('./gamma/WAVECAR')
    # phi = xx.wfc_r(1, 1, 317, ngrid=(60, 108, 160),
    #                gamma=True)
    # xx.save2vesta(phi, poscar='./gamma/POSCAR',gamma=True)

    # xx = vaspwfc('WAVECAR')
    # dE, ovlap, tdm = xx.TransitionDipoleMoment([1,30,17], [1,30,18], norm=True)
    # print dE, ovlap.real, np.abs(tdm)**2

    # print xx._recl, xx._nspin, xx._rtag
    # print xx._nkpts, xx._nbands, xx._encut
    # print xx._Acell, xx._Bcell
    # # print np.linalg.norm(xx._Acell, axis=1)
    # print xx._ngrid
    # print xx._bands[0,0,:]
    # print xx._kvecs
    # print xx._kpath
    # b = xx.readBandCoeff(1,1,1)
    # xx = np.savetxt('kaka.dat', xx.gvectors(2), fmt='%5d')
    # gvec = xx.gvectors(1)
    # gvec %= xx._ngrid[np.newaxis, :]
    # print gvec

    # ngrid=(28, 28, 252)
    # phi = xx.wfc_r(1, 30, 17, ngrid=(28, 28, 252))
    # header = open('POSCAR').read()
    # with open('wave_real.vasp', 'w') as out:
    #     out.write(header)
    #     out.write('%5d%5d%5d\n' % (ngrid[0], ngrid[1], ngrid[2]))
    #     nwrite=0
    #     for kk in range(ngrid[2]):
    #         for jj in range(ngrid[1]):
    #             for ii in range(ngrid[0]):
    #                 nwrite += 1
    #                 out.write('%22.16f ' % phi.real[ii,jj,kk])
    #                 if nwrite % 10 == 0:
    #                     out.write('\n')
    # with open('wave_imag.vasp', 'w') as out:
    #     out.write(header)
    #     out.write('%5d%5d%5d\n' % (ngrid[0], ngrid[1], ngrid[2]))
    #     nwrite=0
    #     for kk in range(ngrid[2]):
    #         for jj in range(ngrid[1]):
    #             for ii in range(ngrid[0]):
    #                 nwrite += 1
    #                 out.write('%22.16f ' % phi.imag[ii,jj,kk])
    #                 if nwrite % 10 == 0:
    #                     out.write('\n')

    # xx = vaspwfc('wave_tyz')
    # ipr = xx.inverse_participation_ratio()
    # print xx._nbands, xx._nkpts
    #
    # import matplotlib as mpl
    # import matplotlib.pyplot as plt
    #
    # fig = plt.figure()
    # ax = plt.subplot()
    #
    # ax.scatter(ipr[...,0], ipr[..., 1], s=ipr[..., 2] / ipr[..., 2].max() * 10, c=ipr[..., 2],
    #            cmap='jet_r')
    #
    # plt.show()

    xx = vaspwfc("examples/wfc_r/wavecar_mose2-wse2", lsorbit=True)
    phi_spinor = xx.wfc_r(1, 1, 36, ngrid=xx._ngrid * 2)
    for ii in range(2):
        phi = phi_spinor[ii]
        prefix = "spinor_{:02d}".format(ii)
        xx.save2vesta(phi, prefix=prefix, poscar="examples/wfc_r/poscar_mose2-wse2")

    pass
