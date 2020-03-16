from jarvis.core.atoms import Atoms
import os
import json
import shutil


class Wannier90win(object):
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
        if self.semi_core_states is None:
            path_semi_core = str(
                os.path.join(os.path.dirname(__file__), "default_semicore.json")
            )
            f = open(path_semi_core, "r")
            semi_core_states = json.load(f)
            f.close()
            self.semi_core_states = semi_core_states

    def write_win(self, name="win.input"):

        if self.soc == True:
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
            dmax = self.efermi - tol
            dmin = self.efermi + tol
            line = str("dis_froz_max =") + str(dmax) + str("\n")
            f.write(line)
            line = str("dis_froz_min =") + str(dmin) + str("\n")
            f.write(line)

        line = str("dis_num_iter =") + str(self.dis_num_iter) + str("\n")
        f.write(line)
        line = str("dis_mix_ratio =") + str(self.dis_mix_ratio) + str("\n")
        f.write(line)
        line = str("num_iter =") + str(self.num_iter) + str("\n")
        f.write(line)
        line = str("num_print_cycles =") + str(self.num_print_cycles) + str("\n")
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
        f = open(prev_win, "r")
        lines = f.read().splitlines()
        f.close()
        bak = str(prev_win) + str(".bak")
        shutil.copy2(prev_win, bak)

        # line=str('write_hr=.true. \n')
        line = str(hr_tag) + str(".true. \n")

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
    box = [[2.715, 2.715, 0], [0, 2.715, 2.715], [2.715, 0, 2.715]]
    coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
    elements = ["Si", "Si"]
    Si = Atoms(lattice_mat=box, coords=coords, elements=elements)
    Wannier90win(struct=Si, efermi=0.0).write_win()
    Wannier90win().write_hr_win(prev_win="win.input")
"""
