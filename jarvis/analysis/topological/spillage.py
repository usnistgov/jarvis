"""
Code to calculate spin-orbit spillage.

Please find details in:
https://www.nature.com/articles/s41598-019-45028-y
https://www.nature.com/articles/s41524-020-0319-4
"""

import numpy as np
from jarvis.io.vasp.outputs import Wavecar  # vaspwfc


class Spillage(object):
    """
    Spin-orbit spillage criteria.

    Predict whether a material is topologically non-trival.
    The spillage criteria physically signifies
    number of band-inverted electrons.
    A non-zero, high value (generally >0.5) suggests non-trivial behavior
    """

    def __init__(self, wf_noso="", wf_so=""):
        """
        Require path to WAVECAR files with and without LSORBIT = .TRUE.

        Args:
            wf_noso : WAVECAR without spin-orbit coupling

            wf_so : WAVECAR with spin-orbit coupling
        """
        self.wf_noso = wf_noso
        self.wf_so = wf_so

    def isclose(self, n1, n2, rel_tol=1e-7):
        """Check if n1 and n2 are close."""
        if abs(n1 - n2) < rel_tol:
            return True
        else:
            return False

    def orth(self, A):
        """Orthogonalize a vector."""
        # As we do not store overlap matrix, we orthogonalize the
        # Wavefunctions.
        u, s, vh = np.linalg.svd(A, full_matrices=False)
        M, N = A.shape
        eps = np.finfo(float).eps
        tol = max(M, N) * np.amax(s) * eps
        num = np.sum(s > tol, dtype=int)
        Q = u[:, :num]
        return Q, num

    def overlap_so_spinpol(self):
        """Calculate spillage."""
        noso = Wavecar(filename=self.wf_noso, lsorbit=False)
        so = Wavecar(filename=self.wf_so, lsorbit=True)
        # band gap stuff. not needed per se, just a useful sanity check
        noso_k, noso_bands = noso.readWFBand()
        so_k, so_bands = so.readWFBand()

        # Calculate the number of occupied bands in non-soc calculation
        # at each k-point

        # Start from deep energy levels and then finds
        # first unoccupied bands with occupation less than 50%
        # noso._kvecs[nk1 - 1, :]: number of kpoints are rows,
        # number of columns three, kpoint units

        nelec_list = []
        for nk1 in range(1, noso._nkpts + 1):  # no spin orbit kpoints loop
            knoso = noso._kvecs[nk1 - 1, :]
            for nk2 in range(1, so._nkpts + 1):  # spin orbit
                kso = so._kvecs[nk2 - 1, :]
                if (
                    self.isclose(kso[0], knoso[0])
                    and self.isclose(kso[1], knoso[1])
                    and self.isclose(kso[2], knoso[2])
                ):  # do kpoints match?
                    if nk1 == nk2:
                        for c, e in enumerate(noso._occs[0, nk1 - 1, :]):
                            if e < 0.5:
                                cup = c
                                break
                        for c, e in enumerate(noso._occs[1, nk1 - 1, :]):
                            if e < 0.5:
                                cdn = c
                                break

                        nelec_list.append([cup, cdn, cup + cdn])
                        # total number of occupied bands at k point nk1

        n_arr = np.array(nelec_list)

        n_up = int(round(np.mean(n_arr[:, 0])))
        n_dn = int(round(np.mean(n_arr[:, 1])))
        n_tot = int(round(np.mean(n_arr[:, 2])))

        nelec = int(n_tot)

        # noso_homo_up = np.max(noso_bands[0, :, n_up - 1])
        # noso_lumo_up = np.min(noso_bands[0, :, n_up])
        # noso_homo_dn = np.max(noso_bands[1, :, n_dn - 1])
        # noso_lumo_dn = np.min(noso_bands[1, :, n_dn])

        so_homo = np.max(so_bands[0, :, nelec - 1])
        so_lumo = np.min(so_bands[0, :, nelec])

        # noso_direct_up = np.min(noso_bands[0, :, n_up] -
        # noso_bands[0, :, n_up - 1])
        # noso_direct_dn = np.min(noso_bands[1, :, n_dn] -
        # noso_bands[1, :, n_dn - 1])

        so_direct = np.min(so_bands[0, :, nelec] - so_bands[0, :, nelec - 1])

        noso_direct = 1000000.0
        noso_homo = -10000000.0
        noso_lumo = 100000000.0
        for i in range(noso_bands.shape[1]):
            homo_k = max(
                noso_bands[0, i, n_up - 1], noso_bands[1, i, n_dn - 1]
            )
            lumo_k = min(noso_bands[0, i, n_up], noso_bands[1, i, n_dn])
            noso_direct = min(noso_direct, lumo_k - homo_k)

            noso_homo = max(noso_homo, homo_k)
            noso_lumo = min(noso_lumo, lumo_k)

        gamma_k = []
        kpoints = []
        np.set_printoptions(precision=4)

        x = []
        y = []

        nelec_tot = 0.0
        for nk1 in range(1, noso._nkpts + 1):  # no spin orbit kpoints loop
            knoso = noso._kvecs[nk1 - 1, :]
            for nk2 in range(1, so._nkpts + 1):  # spin orbit
                kso = so._kvecs[nk2 - 1, :]
                if (
                    self.isclose(kso[0], knoso[0])
                    and self.isclose(kso[1], knoso[1])
                    and self.isclose(kso[2], knoso[2])
                ):  # do kpoints match?

                    # changes section 2
                    if nk1 == nk2:
                        nelec_up = n_arr[nk1 - 1, 0]
                        nelec_dn = n_arr[nk1 - 1, 1]
                        nelec_tot = n_arr[nk1 - 1, 2]

                        gamma_k.append(nelec_tot)
                        kpoints.append(kso)
                        Mmn = 0.0
                        vnoso = noso.readBandCoeff(
                            ispin=1, ikpt=nk1, iband=1, norm=False
                        )
                        # Getting size of matrices
                        n_noso1 = vnoso.shape[0]
                        # number of plane waves in noso
                        vnoso = noso.readBandCoeff(
                            ispin=2, ikpt=nk1, iband=1, norm=False
                        )
                        # spin, kpt, bands, number of plane wave
                        # n_noso2 = vnoso.shape[0]
                        vso = so.readBandCoeff(
                            ispin=1, ikpt=nk2, iband=1, norm=False
                        )
                        n_so = vso.shape[0]  # number of plane waves in so
                        vs = min(n_noso1 * 2, n_so)
                        # Vnono,and so holds wavefunctions
                        Vnoso = np.zeros((vs, nelec_tot), dtype=complex)
                        Vso = np.zeros((vs, nelec_tot), dtype=complex)

                        # prepare matricies,
                        # putting the wavefunction coeffients
                        # align so that they have same structure as SOC case
                        # which has both spin and down together

                        for n1 in range(1, nelec_up + 1):
                            Vnoso[0 : vs // 2, n1 - 1] = noso.readBandCoeff(
                                ispin=1, ikpt=nk1, iband=n1, norm=False
                            )[0 : vs // 2]
                        for n1 in range(1, nelec_dn + 1):
                            Vnoso[
                                vs // 2 : vs, n1 - 1 + nelec_up
                            ] = noso.readBandCoeff(
                                ispin=2, ikpt=nk1, iband=n1, norm=False
                            )[
                                0 : vs // 2
                            ]

                        for n1 in range(1, nelec_tot + 1):
                            t = so.readBandCoeff(
                                ispin=1, ikpt=nk2, iband=n1, norm=False
                            )
                            Vso[0 : vs // 2, n1 - 1] = t[0 : vs // 2]
                            Vso[vs // 2 : vs, n1 - 1] = t[
                                n_so // 2 : n_so // 2 + vs // 2
                            ]
                        # make orthonormal basis?
                        # Occupied bands, both up and down
                        # Generally to get wavefunctions, we have
                        # <psi><S><psi*>
                        # But as we do't have S from DFT saved
                        # we orthogalize the wavefunctions hoping
                        # they span the same space
                        Qnoso, num_noso = self.orth(Vnoso)

                        Qso, num_so = self.orth(Vso)

                        # evalute main expression

                        a = []
                        for n1 in range(0, nelec_tot):  # noso occupied bands
                            v1 = Qnoso[:, n1]  # number of plane waves, bands
                            aa = 0.0
                            # represents the non band-inverted electrons
                            for n2 in range(0, nelec_tot):  # so occupied bands
                                v2 = Qso[:, n2]
                                # Inner product of a nonsoc
                                # and soc wavefunction
                                # index n1 and n2

                                # The inner product
                                # of bands with and without SOC
                                # will be 1 for trivial bands

                                # Delta function of n1 and n2 if SOC is weak

                                t = np.dot(np.conj(v1), v2)
                                Mmn += t * t.conj()
                                aa += (t * t.conj()).real
                            a.append(aa)
                        gamma_k[-1] -= Mmn  # eq 4 in prb 90 125133
                        x.append(kso)
                        y.append(np.real(gamma_k[-1]))
                        if gamma_k[-1] > 0.5:
                            print(
                                "nk1 nk2 kpoint gamma_k ",
                                nk1,
                                nk2,
                                kso,
                                knoso,
                                np.real(gamma_k[-1]),
                                "!!!!!!!!!!",
                            )

        gmax = max(np.real(gamma_k))
        nkmax = np.argmax(np.real(gamma_k))
        kmax = kpoints[nkmax]

        print("------------------------------------")
        print()
        print("                   INDIRECT DIRECT      HOMO/LUMO (eV)")
        print(
            "no spin-orbit gaps",
            "{:+.3f}".format(float(noso_lumo - noso_homo)),
            "{:+.3f}".format(noso_direct),
            "   ",
            [noso_homo, noso_lumo],
        )
        print(
            "spin-orbit gaps   ",
            "{:+.3f}".format(float(so_lumo - so_homo)),
            "{:+.3f}".format(so_direct),
            "   ",
            [so_homo, so_lumo],
        )
        print()
        info = {}
        print("gamma max", np.real(gmax), " at k =  ", kmax)
        info["spillage"] = np.real(gmax)
        info["spillage_k"] = np.real(gamma_k).tolist()
        info["kpoints"] = [kk.tolist() for kk in kpoints]
        info["kmax"] = kmax.tolist()
        info["noso_direct"] = noso_direct
        info["so_direct"] = so_direct
        info["x_axis"] = [xx.tolist() for xx in x]
        info["y_axis"] = [yy.tolist() for yy in y]
        info["so_lumo"] = so_lumo
        info["so_homo"] = so_homo
        info["noso_lumo"] = noso_lumo
        info["noso_homo"] = noso_homo
        return info


# wf_so='/rk2/knc6/JARVIS-DFT/Bulk9-at30/mp-22260_PBEBO/MAIN-SOCSCFBAND-JVASP-59757_mp-22260/WAVECAR'
# wf_noso='/rk2/knc6/JARVIS-DFT/Bulk9-at30/mp-22260_PBEBO/MAIN-MAGSCFBAND-JVASP-59757_mp-22260/WAVECAR'

# spl = Spillage(wf_noso=wf_noso, wf_so=wf_so)
# info = spl.overlap_so_spinpol()
"""
if __name__ == "__main__":
    # JVASP-1044
    wf_noso = "../../tests/testfiles/analysis/spillage/WAVECAR.nosoc"
    wf_so = "../../tests/testfiles/analysis/spillage/WAVECAR.soc"
    spl = Spillage(wf_noso=wf_noso, wf_so=wf_so)
    info = spl.overlap_so_spinpol()
    print(info["spillage"])
"""
