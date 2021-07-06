"""Class for writing Wannier90 input files."""

import os
import json
import shutil
from collections import OrderedDict
from jarvis.core.atoms import Atoms


class Wannier90win(object):
    """Class for writing wannier90.win."""

    def __init__(
        self,
        struct="",
        dis_num_iter=100,
        kmesh_tol=0.000000001,
        dis_mix_ratio=0.5,
        num_iter=100,
        num_print_cycles=50,
        frozen_tol=2,
        efermi=None,
        soc=True,
        semi_core_states=None,
    ):
        """

        Provide information ndeeded for writing wannier90.win.

        At least Atoms object needed.

        Args:
            struct : Atoms object

            efermi : Fermi-energy

            soc : spin-orbit coupling True/False

            dis_num_iter : number of disentanglement iteration steps
                           generally 100-200 is enough

            num_iter : number of wannierisation iterations

            kmesh_tol : important for wannier convergence,
                        try to increase/decrease if
                        wannierization fails

            dis_mix_ratio : mixing ratio for the disentanglement routine

            frozen_tol :  disentanglement inner (frozen) window with +/- values

            semi_core_states : important to specify according to
                              pseudopotentials. A generic one is
                              given for GGA-PBE, change if you use
                              different POTCARs
        """
        self.struct = struct
        self.efermi = efermi
        self.soc = soc
        self.dis_num_iter = dis_num_iter
        self.dis_mix_ratio = dis_mix_ratio
        self.num_iter = num_iter
        self.num_print_cycles = num_print_cycles
        self.frozen_tol = frozen_tol
        self.semi_core_states = semi_core_states
        self.kmesh_tol = kmesh_tol
        # Check default potcars after Er in the json
        # Not tested completely
        if self.semi_core_states is None:
            path_semi_core = str(
                os.path.join(
                    os.path.dirname(__file__), "default_semicore.json"
                )
            )
            f = open(path_semi_core, "r")
            semi_core_states = json.load(f)
            f.close()
            self.semi_core_states = semi_core_states

    def to_dict(self):
        """Convert to a dictionary."""
        d = OrderedDict()
        d["struct"] = self.struct.to_dict()
        d["efermi"] = self.efermi
        d["soc"] = self.soc
        d["dis_num_iter"] = self.dis_num_iter
        d["dis_mix_ratio"] = self.dis_mix_ratio
        d["num_iter"] = self.num_iter
        d["num_print_cycles"] = self.num_print_cycles
        d["frozen_tol"] = self.frozen_tol
        d["semi_core_states"] = self.semi_core_states
        d["kmesh_tol"] = self.kmesh_tol
        return d

    @classmethod
    def from_dict(self, d={}):
        """Convert class from a dictionary."""
        return Wannier90win(
            struct=Atoms.from_dict(d["struct"]),
            efermi=d["efermi"],
            soc=d["soc"],
            dis_num_iter=d["dis_num_iter"],
            dis_mix_ratio=d["dis_mix_ratio"],
            num_iter=d["num_iter"],
            num_print_cycles=d["num_print_cycles"],
            frozen_tol=d["frozen_tol"],
            semi_core_states=d["semi_core_states"],
            kmesh_tol=d["kmesh_tol"],
        )

    def write_win(self, name="win.input"):
        """Write .win file."""
        if self.soc:
            mult = 2
        else:
            mult = 1
        semi_core = self.semi_core_states

        xx = self.struct.composition.to_dict()
        elements = []
        ELEMENTCOUNT = []
        for i, j in xx.items():
            elements.append(i)
            ELEMENTCOUNT.append(int(j))
        nwan = 0
        nelectrons = 0
        exclude = 0
        projections = ""
        c = 0
        for e in elements:
            s = semi_core[e]
            exclude += s[0] * ELEMENTCOUNT[c]
            nwan += mult * s[3] * ELEMENTCOUNT[c]
            nelectrons += s[1] * ELEMENTCOUNT[c]
            projections += e + ":" + s[2] + "\n"
            c += 1
        print("excluded", exclude)
        nwan_tot = nwan + exclude
        print("nwan_tot", nwan_tot)
        print("projections", projections)
        if exclude == 0:
            exclude_st = ""
        else:
            exclude_st = "exclude_bands     = 1-" + str(int(exclude)) + "\n"

        f = open(name, "w")
        line = "guiding_centres = TRUE \n"
        f.write(line)
        line = "kmesh_tol=" + str(self.kmesh_tol) + "\n"
        f.write(line)
        # line=str('num_wann =')+str(nwan_tot)+'\n'
        line = str("num_wann =") + str(nwan) + "\n"
        f.write(line)
        f.write(exclude_st)

        tol = self.frozen_tol
        if self.efermi is not None:
            dmax = self.efermi + tol
            # dmin = self.efermi - tol
            line = str("dis_froz_max =") + str(dmax) + str("\n")
            f.write(line)
            line = str("dis_froz_min =") + str(-1000) + str("\n")
            # line = str("dis_froz_min =") + str(dmin) + str("\n")
            f.write(line)

        line = str("dis_num_iter =") + str(self.dis_num_iter) + str("\n")
        f.write(line)
        line = str("dis_mix_ratio =") + str(self.dis_mix_ratio) + str("\n")
        f.write(line)
        line = str("num_iter =") + str(self.num_iter) + str("\n")
        f.write(line)
        line = (
            str("num_print_cycles =") + str(self.num_print_cycles) + str("\n")
        )
        f.write(line)
        line = str("begin projections") + str("\n")
        f.write(line)
        f.write(projections)
        line = str("end projections") + str("\n")
        f.write(line)
        f.close()
        return nwan, exclude

    def write_hr_win(
        self,
        hr_tag="hr_plot",
        prev_win="wannier90.win",
        hr="hr_wannier.win",
        nbands=18,
        soc="",
    ):
        """Write hr_plot or write_hr as .true."""
        f = open(prev_win, "r")
        lines = f.read().splitlines()
        f.close()
        bak = str(prev_win) + str(".bak")
        shutil.copy2(prev_win, bak)

        # line=str('write_hr=.true. \n')
        line = str(hr_tag) + str("=.true. \n")

        f = open(hr, "w")
        f.write(line)
        line = str("num_bands = ") + str(nbands) + ("\n")
        f.write(line)
        for i in lines:
            line = str(i) + str("\n")
            f.write(line)
        f.close()


"""
if __name__ == "__main__":
    from jarvis.core.atoms import Atoms
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Wannier90win(struct=Si, efermi=0.0).write_win()
    Wannier90win().write_hr_win(prev_win="win.input")
"""
