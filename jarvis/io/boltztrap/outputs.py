"""Class for analyzing BoltzTrap outputs."""

import os
import numpy as np
from collections import OrderedDict, defaultdict

Ry_to_ev = 13.6056980659


class BoltzTrapOutput(object):
    """Analyze BoltzTrap output."""

    def __init__(
        self,
        path="",
        outtrans_data=[],
        intrans_data=[],
        condtens_fixdoping=[],
        halltens_fixdoping=[],
    ):
        """Specify boltztrap folder  as path to analze data."""
        self.path = path
        if outtrans_data == []:
            outtrans_data = self.read_outputtrans()
        self.outtrans_data = outtrans_data
        if intrans_data == []:
            intrans_data = self.read_intrans()
        self.intrans_data = intrans_data
        if condtens_fixdoping == []:
            condtens_fixdoping = self.read_condtens_fixdoping()
        self.condtens_fixdoping = condtens_fixdoping
        if halltens_fixdoping == []:
            halltens_fixdoping = self.read_halltens_fixdoping()
        self.halltens_fixdoping = halltens_fixdoping

    def to_dict(self):
        """Return output as a dictionary."""
        d = OrderedDict()
        d["path"] = self.path
        d["outtrans_data"] = self.outtrans_data
        d["intrans_data"] = self.intrans_data
        d["halltens_fixdoping"] = self.halltens_fixdoping
        d["condtens_fixdoping"] = self.condtens_fixdoping
        return d

    def read_intrans(self, filename=""):
        """Read intrans file."""
        if filename == "":
            filename = os.path.join(self.path, "boltztrap.intrans")
        f = open(filename, "r")
        lines = f.read().splitlines()
        f.close()
        self.intrans_data = lines
        return lines

    def read_outputtrans(self, filename=""):
        """Read outtrans file."""
        gap = "na"
        Ef = "na"
        vbm = "na"
        cbm = "na"
        warning = "na"
        excessN_doping = "na"
        if filename == "":
            filename = os.path.join(self.path, "boltztrap.outputtrans")
        info = {}
        excessN_doping = {}
        f = open(filename, "r")
        lines = f.read().splitlines()
        f.close()
        for ii, i in enumerate(lines):
            if "WARNING" in i:
                warning = i
            elif "Doping level number" in i:
                tmp_dop = float(i.split("=")[1].split()[0])
                tmp_excess = float(lines[ii + 1].split()[3])
                excessN_doping.setdefault(tmp_excess, tmp_dop)
            elif "Egap" in i:
                gap = float(i.split()[1]) * Ry_to_ev
            elif "VBM" in i:
                vbm = float(i.split()[1]) * Ry_to_ev
                cbm = float(i.split()[3]) * Ry_to_ev
                Ef = float(i.split()[5]) * Ry_to_ev

        info["gap"] = gap
        info["Ef"] = Ef
        info["vbm"] = vbm
        info["cbm"] = cbm
        info["warning"] = warning
        info["excessN_doping"] = excessN_doping
        # self.outtrans_data=info
        return info

    def dopinglevel_for_excessN(self, excessN):
        """Return doping level for excees concentration."""
        excessN_doping = self.read_outputtrans()["excessN_doping"]
        for i, j in excessN_doping.items():
            if np.isclose(i, excessN):
                return j

    def read_condtens_fixdoping(self, filename=""):
        """Read condtens_fixdoping file."""
        if filename == "":
            filename = os.path.join(self.path, "boltztrap.condtens_fixdoping")
        f = open(filename, "r")
        lines = f.read().splitlines()
        f.close()
        full_doping_data = []
        for i in lines:
            if "#" not in i and len(i) > 2:
                full_doping_data.append([float(j) for j in i.split()])
        d = np.array(full_doping_data)
        all_data = {}
        p_dict = defaultdict(dict)
        n_dict = defaultdict(dict)
        for i in d:
            T = i[0]
            ef = i[29] * Ry_to_ev
            N = i[1]
            N_cm3 = self.dopinglevel_for_excessN(N)
            # print ('N,N_cm3',N,N_cm3)
            cond = i[2:11]
            seeb = i[11:20]
            kappa = i[20:29]
            if N > 0:
                info = {}
                info["cond"] = cond
                info["seeb"] = seeb
                info["kappa"] = kappa
                info["Ef"] = ef
                info["N_cm3"] = N_cm3
                p_dict[T][N] = info
            else:
                info = {}
                info["cond"] = cond
                info["seeb"] = seeb
                info["kappa"] = kappa
                info["Ef"] = ef
                info["N_cm3"] = N_cm3
                n_dict[T][N] = info
        all_data["p"] = p_dict
        all_data["n"] = n_dict
        # self.condtens_fixdoping=all_data
        return all_data

    def read_halltens_fixdoping(self, filename=""):
        """Read halltens file."""
        if filename == "":
            filename = os.path.join(self.path, "boltztrap.halltens_fixdoping")
        f = open(filename, "r")
        lines = f.read().splitlines()
        f.close()
        full_doping_data = []
        for i in lines:
            if "#" not in i and len(i) > 2:
                full_doping_data.append([float(j) for j in i.split()])
        d = np.array(full_doping_data)
        all_data = {}
        p_dict = defaultdict(dict)
        n_dict = defaultdict(dict)
        for i in d:
            T = i[0]
            ef = i[29] * Ry_to_ev
            N = i[1]
            N_cm3 = self.dopinglevel_for_excessN(N)
            # print ('N,N_cm3',N,N_cm3)
            hall = i[2:29]
            if N > 0:
                info = {}
                info["hall"] = hall
                info["Ef"] = ef
                info["N_cm3"] = N_cm3
                p_dict[T][N] = info
            else:
                info = {}
                info["hall"] = hall
                info["Ef"] = ef
                info["N_cm3"] = N_cm3
                n_dict[T][N] = info
        all_data["p"] = p_dict
        all_data["n"] = n_dict
        # self.halltens_fixdoping=all_data
        return all_data

    @classmethod
    def from_dict(self, d={}):
        """Load from a dictionary."""
        return BoltzTrapOutput(
            path=d["path"],
            outtrans_data=d["outtrans_data"],
            intrans_data=d["intrans_data"],
            halltens_fixdoping=d["halltens_fixdoping"],
            condtens_fixdoping=d["condtens_fixdoping"],
        )


"""
if __name__ == "__main__":
    condtens_fix = "boltztrap/boltztrap.condtens_fixdoping"
    b = BoltzTrapOutput(
        path="boltztrap"
    ).read_condtens_fixdoping()
    b = BoltzTrapOutput(
        path="boltztrap"
    ).read_outputtrans()
    b = BoltzTrapOutput(
        path="boltztrap"
    ).to_dict()
    print(b)
"""
