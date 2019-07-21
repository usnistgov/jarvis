from __future__ import unicode_literals

# source ~/anaconda2/envs/my_boltz2/bin/activate my_boltz2
import matplotlib.pyplot as plt

plt.switch_backend("Agg")
import unittest
import json, os
import warnings
from monty.serialization import loadfn, dumpfn, MontyEncoder, MontyDecoder
from pymatgen.io.vasp import Vasprun
from interruptingcow import timeout
from pymatgen.electronic_structure.boltztrap2 import (
    BandstructureLoader,
    VasprunLoader,
    BztInterpolator,
    BztTransportProperties,
    BztPlotter,
    merge_up_down_doses,
)


from pymatgen.electronic_structure.core import Spin, OrbitalType

import numpy as np
from monty.serialization import loadfn


vrunfile = "pymatgen-master/test_files/boltztrap2/vasprun.xml"
vrunfile1 = "/users/knc6/Software/JARVIS-TOOLS/jarvis/jarvis/vasp/examples/SiOptb88/MAIN-RELAX-bulk@mp_149/vasprun.xml"
vrunfile2 = "/users/knc6/Software/JARVIS-TOOLS/jarvis/jarvis/vasp/examples/SiOptb88/MAIN-MBJ-bulk@mp_149/vasprun.xml"
vrunfile3 = "/users/knc6/Software/JARVIS-TOOLS/jarvis/jarvis/vasp/examples/SiOptb88/MAIN-OPTICS-bulk@mp_149/vasprun.xml"


@timeout(3600)
def nonspin_boltz(vrunfile="", k_latt=1.0, write_json=False):
    fname = vrunfile.replace("vasprun.xml", "boltz2data.json")
    if not os.path.isfile(fname):
        vrun = Vasprun(vrunfile, parse_projected_eigen=True)
        loader = VasprunLoader().from_file(vrunfile)
        temp_r = np.array([300, 400, 500, 600, 700, 800])
        doping = np.array([0, 10 ** 18, 10 ** 19, 10 ** 20, 10 ** 21, 10 ** 22])
        bztInterp = BztInterpolator(loader, lpfac=2)
        bztTransp = BztTransportProperties(bztInterp, temp_r=temp_r, doping=doping)
        xx = bztTransp.compute_properties_doping(doping=doping)
        # 4 temps, 2 doping
        Conductivity_doping = bztTransp.Conductivity_doping
        Seebeck_doping = bztTransp.Seebeck_doping
        Kappa_doping = bztTransp.Kappa_doping
        Effective_mass_doping = bztTransp.Effective_mass_doping
        Power_Factor_doping = bztTransp.Power_Factor_doping
        mu_r_eV = bztTransp.mu_r_eV
        info = {}
        info["mu_r_eV"] = mu_r_eV
        info["temp_r"] = temp_r
        info["doping"] = doping
        info["Conductivity_doping"] = Conductivity_doping
        info["Seebeck_doping"] = Seebeck_doping
        info["Kappa_doping"] = Kappa_doping
        info["Effective_mass_doping"] = Effective_mass_doping
        info["Power_Factor_doping"] = Power_Factor_doping

        info["Conductivity_mu"] = bztTransp.Conductivity_mu
        info["Seebeck_mu"] = bztTransp.Seebeck_mu
        info["Kappa_mu"] = bztTransp.Kappa_mu
        info["Power_Factor_mu"] = bztTransp.Power_Factor_mu
        info["Effective_mass_mu"] = bztTransp.Effective_mass_mu
        info["Hall_carrier_conc_trace_mu"] = bztTransp.Hall_carrier_conc_trace_mu

        zt = []
        temp_zt = []
        for i, ii in enumerate(info["temp_r"]):
            for j, k in zip(info["Power_Factor_mu"][i], info["Kappa_mu"][i]):
                temp_zt.append(0.001 * j * ii / (k + k_latt))
            zt.append(temp_zt)
            temp_zt = []
        zt = np.array(zt)
        info["zt_mu"] = zt
        if write_json == True:
            f = open(fname, "w")
            f.write(json.dumps(info, cls=MontyEncoder))
            f.close()

        return info
    else:
        print("File exists")
    # print (bztTransp.Effective_mass_mu[0][0])
    # bztPlotter = BztPlotter(bztTransp,bztInterp)
    # fig = bztPlotter.plot_props('S', 'mu', 'temp', temps=[300, 500])
    # fig.savefig('first_boltz2.png')
    # fig.close()


@timeout(3600)
def spin_boltz(vrunfile="", spin=1, k_latt=1.0, write_json=True):
    fname = vrunfile.replace("vasprun.xml", "boltz2data.json")
    if not os.path.isfile(fname):
        kp = vrunfile.replace("vasprun.xml", "KPOINTS")
        v = Vasprun(vrunfile)
        nelect = v.parameters["NELECT"]
        bs = v.get_band_structure(kp, line_mode=False)
        # doping=10.**np.arange(20,22)
        temp_r = np.array([300, 400, 500, 600, 700, 800])
        doping = np.array([0, 10 ** 18, 10 ** 19, 10 ** 20, 10 ** 21, 10 ** 22])
        loader = BandstructureLoader(bs, v.structures[-1], spin=spin, nelect=nelect)
        bztInterp = BztInterpolator(loader, lpfac=2, curvature=True)
        bztTransp = BztTransportProperties(bztInterp, doping=doping, temp_r=temp_r)
        xx = bztTransp.compute_properties_doping(doping=doping)
        # 4 temps, 2 doping
        Conductivity_doping = bztTransp.Conductivity_doping
        Seebeck_doping = bztTransp.Seebeck_doping
        Kappa_doping = bztTransp.Kappa_doping
        Effective_mass_doping = bztTransp.Effective_mass_doping
        Power_Factor_doping = bztTransp.Power_Factor_doping
        mu_r_eV = bztTransp.mu_r_eV

        info = {}
        info["mu_r_eV"] = mu_r_eV
        info["temp_r"] = temp_r
        info["doping"] = doping
        info["Conductivity_doping"] = Conductivity_doping
        info["Seebeck_doping"] = Seebeck_doping
        info["Kappa_doping"] = Kappa_doping
        info["Effective_mass_doping"] = Effective_mass_doping
        info["Power_Factor_doping"] = Power_Factor_doping

        info["Conductivity_mu"] = bztTransp.Conductivity_mu
        info["Seebeck_mu"] = bztTransp.Seebeck_mu
        info["Kappa_mu"] = bztTransp.Kappa_mu
        info["Power_Factor_mu"] = bztTransp.Power_Factor_mu
        info["Effective_mass_mu"] = bztTransp.Effective_mass_mu
        info["Hall_carrier_conc_trace_mu"] = bztTransp.Hall_carrier_conc_trace_mu
        zt = []
        temp_zt = []
        for i, ii in enumerate(info["temp_r"]):
            for j, k in zip(info["Power_Factor_mu"][i], info["Kappa_mu"][i]):
                temp_zt.append(0.001 * j * ii / (k + k_latt))
            zt.append(temp_zt)
            temp_zt = []
        zt = np.array(zt)
        info["zt_mu"] = zt
        if write_json == True:
            f = open(fname, "w")
            f.write(json.dumps(info, cls=MontyEncoder))
            f.close()
        return info
    else:
        print("File exists")


def run_dir(
    dir="/users/knc6/Software/JARVIS-TOOLS/jarvis/jarvis/vasp/examples/SiOptb88"
):
    for i in os.listdir(dir):
        if "json" in i:
            if "MAIN-RELAX" in i:
                try:
                    mvrun = (
                        str(dir)
                        + str("/")
                        + str(i.split(".json")[0])
                        + str("/vasprun.xml")
                    )
                    info = spin_boltz(mvrun)
                    print("mvrun", mvrun, info)
                except:
                    print("skpiing")
                    pass
            if "MAIN-MBJ" in i:
                try:
                    mbjvrun = (
                        str(dir)
                        + str("/")
                        + str(i.split(".json")[0])
                        + str("/vasprun.xml")
                    )
                    info = spin_boltz(mbjvrun)
                    print("mbjvrun", mbjvrun, info)
                except:
                    print("skpiing")
                    pass
            if "MAIN-OPTICS" in i:
                try:
                    mopvrun = (
                        str(dir)
                        + str("/")
                        + str(i.split(".json")[0])
                        + str("/vasprun.xml")
                    )
                    info = spin_boltz(mopvrun)
                    print("mopvrun", mopvrun, info)
                except:
                    print("skpiing")
                    pass


def get_specific_data(my_type="p", my_dop=10 ** 20, my_temp=600, bfile=""):
    dop = "na"
    temp = "na"
    d = loadfn(bfile, cls=MontyDecoder)
    for i, ii in enumerate(d["temp_r"]):
        if ii == my_temp:
            temp = i
    for j, jj in enumerate(d["doping"]):
        if jj == my_dop:
            dop = j
    if dop != "na" and temp != "na":
        seeb = max(d["Seebeck_doping"][my_type][temp][dop].flatten())
        cond = max(d["Conductivity_doping"][my_type][temp][dop].flatten())
        kapp = max(d["Kappa_doping"][my_type][temp][dop].flatten())
        pf = max(d["Power_Factor_doping"][my_type][temp][dop].flatten())
        eff_mass = max(d["Effective_mass_doping"][my_type][temp][dop].flatten())
        return seeb, cond, kapp, pf, eff_mass
    print("Selection does not exist")


def plot_prop_mu(prop="Seebeck_mu", bfile="", temp=300, filename="seebplot.png"):

    props = [
        "Conductivity_mu",
        "Seebeck_mu",
        "Kappa_mu",
        "Effective_mass_mu",
        "Power_Factor_mu",
        "Carrier_conc_mu",
        "Hall_carrier_conc_trace_mu",
        "zt_mu",
    ]
    props_lbl = [
        "Conductivity",
        "Seebeck",
        "$K_{el}$",
        "Effective mass",
        "Power Factor",
        "Carrier concentration",
        "Hall carrier conc.",
        "ZT",
    ]
    props_unit = (
        r"$(\mathrm{kS\,m^{-1}})$",
        r"($\mu$V/K)",
        r"$(W / (m \cdot K))$",
        r"$(m_e)$",
        r"$( mW / (m\cdot K^2)$",
        r"$(cm^{-3})$",
        r"$(cm^{-3})$",
        "",
    )
    ylbl = "na"
    for i, j, k in zip(props, props_lbl, props_unit):
        if i == prop:
            ylbl = j + " " + k
    print("ylbl", ylbl)
    d = loadfn(bfile, cls=MontyDecoder)
    itemp = "na"
    for i, ii in enumerate(d["temp_r"]):
        if ii == temp:
            itemp = i
    data = d[prop][itemp]
    dat_x = []
    dat_y = []
    dat_z = []
    for i in data:
        dat, vecs = np.linalg.eig(i)
        dat_x.append(dat[0])
        dat_y.append(dat[1])
        dat_z.append(dat[2])
    plt.close()
    plt.rcParams.update({"font.size": 22})
    fig = plt.figure(figsize=(10, 8))
    plt.plot(d["mu_r_eV"], dat_x, label="x")
    plt.plot(d["mu_r_eV"], dat_y, label="y")
    plt.plot(d["mu_r_eV"], dat_z, label="z")
    # plt.axhline(y=0.0, color='k', linestyle='-')
    plt.ylabel(ylbl)
    plt.xlim([-3, 3])
    plt.xlabel(r"$\mu$ (eV)")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


if __name__ == "__main__":
    path = path = str(
        os.path.join(os.path.dirname(__file__), "../vasp/examples/SiOptb88")
    )
    run_dir(path)
