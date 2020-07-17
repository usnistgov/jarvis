"""Class for writing wt.in for wanniertools."""

import os
from jarvis.core.kpoints import Kpoints3D
from jarvis.io.wannier.outputs import Wannier90wout
import json
from collections import OrderedDict
from jarvis.core.atoms import Atoms


class WTin(object):
    """Construct wt.in object."""

    def __init__(
        self,
        atoms=None,
        nelect=8,
        miller=[0, 0, 1],
        wannierout="wannier90.wout",
        wtin="wt.in",
        efermi=None,
        semi_core_states=None,
        soc=True,
        exclude=0,
        nwan=10,
    ):
        """Initialize class."""
        self.atoms = atoms
        self.nelect = nelect
        self.wannierout = wannierout
        self.efermi = efermi
        self.soc = soc
        self.exclude = exclude
        self.nwan = nwan
        self.wtin = wtin
        self.miller = " ".join(map(str, miller))
        self.semi_core_states = semi_core_states
        if self.semi_core_states is None:
            path_semi_core = str(
                os.path.join(
                    os.path.dirname(__file__),
                    "..",
                    "wannier",
                    "default_semicore.json",
                )
            )
            f = open(path_semi_core, "r")
            semi_core_states = json.load(f)
            f.close()
            self.semi_core_states = semi_core_states

    def to_dict(self):
        """Convert the class to a dictionary."""
        d = OrderedDict()
        d["atoms"] = self.atoms.to_dict()
        d["nelect"] = self.nelect
        d["wannierout"] = self.wannierout
        d["efermi"] = self.efermi
        d["soc"] = self.soc
        d["exclude"] = self.exclude
        d["nwan"] = self.nwan
        d["wtin"] = self.wtin
        d["miller"] = self.miller
        d["semi_core_states"] = self.semi_core_states
        return d

    @classmethod
    def from_dict(self, d={}):
        """Construct class from a dictionary."""
        return WTin(
            atoms=Atoms.from_dict(d["atoms"]),
            nelect=d["nelect"],
            wannierout=d["wannierout"],
            efermi=d["efermi"],
            soc=d["soc"],
            exclude=d["exclude"],
            nwan=d["nwan"],
            wtin=d["wtin"],
            miller=d["miller"],
            semi_core_states=d["semi_core_states"],
        )

    def get_ibz_kp(self):
        """Get high-symmetry k-points."""
        frac_k_points, k_points_labels = Kpoints3D().interpolated_points(
            self.atoms
        )
        lines = []
        for i, j in zip(frac_k_points, k_points_labels):
            if j != "":
                line = (
                    str(j)
                    + str(" ")
                    + str(i[0])
                    + str(" ")
                    + str(i[1])
                    + str(" ")
                    + str(i[2])
                )
                line = line.replace("\\", "")
                lines.append(line)
        lines1 = lines
        lines2 = lines[1:]
        pairs = zip(lines1, lines2)
        lines = []
        for i in pairs:
            line = i[0] + str(" ") + str(i[1])  # +'\n'
            if i[0] != i[1]:
                lines.append(line)
        return lines

    def semi_core_wt(self, string=""):
        """Get emi-core states."""
        try:
            string = string.replace("p", "px,py,pz")
        except Exception:
            pass
        try:
            string = string.replace("d", "dxy,dxz,dyz,dx2-y2,dz2")
        except Exception:
            pass
        try:
            string = string.replace(
                "f", "fz3,fxz2,fyz2,fxyz,fzx2y2,fxx23y2,fy3x2y2"
            )
        except Exception:
            pass
        return string

    def write_wt_in(self):
        """Write et.in."""
        strt = self.atoms
        # nwan = self.nwan
        exclude = self.exclude

        f = open(self.wtin, "w")
        line = str("&TB_FILE \n")
        f.write(line)
        line = str("Hrfile = 'wannier90_hr.dat' \n")
        f.write(line)
        line = str("Package = 'VASP' \n")
        f.write(line)
        line = str("/ \n")
        f.write(line)

        wan_cts = Wannier90wout(
            wout_path=self.wannierout
        ).give_wannier_centers()

        control = {
            "BulkBand_calc": "T",
            "Z2_3D_calc": "T",
            "Chern_3D_calc": "T",
            "SlabSS_calc": "T",
            "FindNodes_calc": "T",
            "SlabArc_calc": "F",
            "BerryCurvature_calc": "F",
        }
        nele = self.nelect
        if self.soc:
            soc = 1
        else:
            soc = 0
        noc = int(float(nele)) - int(float(exclude))
        system = {
            "NSLAB": 10,
            "NumOccupied": noc,
            "SOC": soc,
            "E_FERMI": self.efermi,
            "Bx": 0,
            "By": 0,
            "Bz": 0,
            "surf_onsite": 0,
        }
        parameters = {
            "Eta_Arc": 0.001,
            "E_arc": 0.0,
            "OmegaNum": 500,
            "OmegaMin": -0.6,
            "OmegaMax": 0.5,
            "Nk1": 21,
            "Nk2": 21,
            "Nk3": 21,
            "NP": 2,
            "Gap_threshold": 0.0001,
        }
        kp_lines = self.get_ibz_kp()
        kp_dat = []
        kp_dat.append(len(kp_lines))
        for ii in kp_lines:
            kp_dat.append(ii)
        kp_dict = {
            "KPATH_BULK": [
                4,
                "G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000",
                "Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000",
                "F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000",
                "G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000",
            ],
            "KPATH_SLAB": [2, "K 0.33 0.67 G 0.0 0.0", "G 0.0 0.0 M 0.5 0.5"],
            "KPLANE_SLAB": ["-0.1 -0.1", "0.2  0.0", "0.0  0.2"],
            "KCUBE_BULK": [
                "-0.50 -0.50 -0.50",
                "1.00  0.00  0.00",
                "0.00  1.00  0.00",
                "0.00  0.00  1.00",
            ],
        }
        kp_dict = {
            "KPATH_BULK": kp_dat,
            "KPATH_SLAB": [2, "K 0.33 0.67 G 0.0 0.0", "G 0.0 0.0 M 0.5 0.5"],
            "KPLANE_SLAB": ["-0.1 -0.1", "0.2  0.0", "0.0  0.2"],
            "KPLANE_BULK": [
                "0.0 0.0  0.00",
                "1.00  0.00  0.00 ",
                " 0.00  1.00  0.00",
            ],
            "KCUBE_BULK": [
                "-0.50 -0.50 -0.50",
                "1.00  0.00  0.00",
                "0.00  1.00  0.00",
                "0.00  0.00  1.00",
            ],
        }
        print(kp_dat)
        wt_dict = {
            "&CONTROL": control,
            "&SYSTEM": system,
            "&PARAMETERS": parameters,
        }

        for i, j in wt_dict.items():
            line = str(i) + "\n"
            f.write(line)
            for k, l in j.items():
                line = str(k) + str(" = ") + str(l) + "\n"
                f.write(line)
            line = str("/") + "\n"
            f.write(line)

        line = str("LATTICE \n")
        f.write(line)
        line = str("Angstrom \n")
        f.write(line)
        line = (
            str(strt.lattice_mat[0][0])
            + str(" ")
            + str(strt.lattice_mat[0][1])
            + str(" ")
            + str(strt.lattice_mat[0][2])
            + str("\n")
        )
        f.write(line)
        line = (
            str(strt.lattice_mat[1][0])
            + str(" ")
            + str(strt.lattice_mat[1][1])
            + str(" ")
            + str(strt.lattice_mat[1][2])
            + str("\n")
        )
        f.write(line)
        line = (
            str(strt.lattice_mat[2][0])
            + str(" ")
            + str(strt.lattice_mat[2][1])
            + str(" ")
            + str(strt.lattice_mat[2][2])
            + str("\n")
        )
        f.write(line)
        line = str("ATOM_POSITIONS \n")
        f.write(line)
        line = str(strt.num_atoms) + str(" \n")
        f.write(line)
        line = str("Direct \n")
        f.write(line)

        for i, j in zip(strt.elements, strt.frac_coords):
            line = (
                str(i)
                + str(" ")
                + str(j[0])
                + str(" ")
                + str(j[1])
                + str(" ")
                + str(j[2])
                + "\n"
            )
            f.write(line)
        line = str("PROJECTORS \n")
        f.write(line)

        # semicore=ZVAL in POTCAR -  actual valence elec.
        semi_core = self.semi_core_states

        prjs = ""
        for i in strt.elements:
            symb = str(i)
            projs = semi_core[symb][2]
            projs = self.semi_core_wt(projs)
            num_projs = len(projs.split(","))
            prjs = prjs + str(" ") + str(num_projs)
        line = str(prjs) + "\n"
        f.write(line)

        for i in strt.elements:
            symb = str(i)
            projs = semi_core[symb][2]
            projs = self.semi_core_wt(projs)
            num_projs = len(projs.split(","))
            line = str(symb) + str(" ") + str(projs) + "\n"
            f.write(line)
        line = str("MILLER_INDICES \n")
        f.write(line)
        line = str(self.miller) + "\n"  # "0 0 1 \n")
        f.write(line)

        for i, j in kp_dict.items():
            line = str(i) + "\n"
            f.write(line)
            for k in j:
                line = str(k) + "\n"
                f.write(line)
        line = str("WANNIER_CENTRES \n")
        f.write(line)
        line = str("Cartesian \n")
        f.write(line)
        for i in wan_cts:
            print("wan cts", i)

            line = (
                str(i[0])
                + str(" ")
                + str(i[1])
                + str(" ")
                + str(i[2])
                + str(" ")
                + str("\n")
            )
            f.write(line)
        f.close()


"""
if __name__ == "__main__":
    from jarvis.io.vasp.inputs import Poscar
    hr = "wannier90_hr.dat"
    wout = "MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/wannier90.wout"
    centers = Wannier90wout(wout_path=wout).give_wannier_centers()

    p = Poscar.from_file(
        "MAIN-WANN-SOC-bulk@JVASP-1067_mp-541837/POSCAR"
    ).atoms
    lines = WTin(atoms=p).get_ibz_kp()
    # print(lines)
    wtin = WTin(atoms=p, wannierout=wout).write_wt_in()
"""
