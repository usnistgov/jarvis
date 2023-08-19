"""Module to parse WIEN2k outputs."""
import numpy as np
import matplotlib.pyplot as plt

plt.switch_backend("agg")


def read_band_energy(
    energy_file="FeSe.energy", plot=True, band_plot="band.png"
):
    """Read case.energy file."""
    f = open(energy_file, "r")
    lines = f.read().splitlines()
    f.close()
    eigs = []
    eig_bands = []
    start = True
    for i in lines:
        sp = i.split()
        if len(sp) == 2 and start:
            tmp = float(sp[1])
            eig_bands.append(tmp)
        if len(sp) == 7:
            start = True
            if eig_bands != []:
                eigs.append(np.array(eig_bands))
            eig_bands = []
    # print('eigs',eigs)
    # eigs = np.array(eigs)
    if plot:
        for ii, i in enumerate(eigs):
            plt.plot([ii for j in range(len(i))], i, ".", c="b")
        plt.savefig(band_plot)
        plt.close()
    return eigs


def read_scf(scf_file="FeSe.scf"):
    """Read case.scf file."""
    f = open(scf_file, "r")
    lines = f.read().splitlines()
    f.close()
    info = {}
    efermies = []
    nelects = []
    nkpts = []
    totens = []
    for i in lines:
        if ":FER  : F E R M I - ENERGY(TETRAH.M.)" in i:
            efermies.append(
                float(i.split(":FER  : F E R M I - ENERGY(TETRAH.M.)=")[1])
            )
        if ":ENE  : ********** TOTAL ENERGY IN Ry" in i:
            totens.append(
                float(i.split(":ENE  : ********** TOTAL ENERGY IN Ry =")[1])
            )
        if ":NOE  : NUMBER OF ELECTRONS          =" in i:
            nelects.append(
                float(i.split(":NOE  : NUMBER OF ELECTRONS          =")[1])
            )
        if ":KPT   :      NUMBER OF K-POINTS:" in i:
            nkpts.append(
                float(i.split(":KPT   :      NUMBER OF K-POINTS:")[1])
            )
    info["efermi"] = efermies[-1]
    info["nelect"] = nelects[-1]
    info["nkpt"] = nkpts[-1]
    info["tot_en"] = totens[-1]
    return info


def band_eigvals(energy_file="FeSe.energy", plot=False, band_file="band.png"):
    """Get bandstructure eigenvalues."""
    f = open(energy_file, "r")
    lines = f.read().splitlines()
    f.close()
    eigs = []
    eig_bands = []
    start = True
    for i in lines:
        sp = i.split()
        if len(sp) == 2 and start:
            tmp = float(sp[1])
            eig_bands.append(tmp)
        if len(sp) > 2:
            # if len(sp) == 7:
            start = True
            if eig_bands != []:
                eigs.append(eig_bands)
            eig_bands = []
    # eigs = np.array(eigs)
    print("TODO:Fix bug in np.array())")
    if plot:
        import matplotlib.pyplot as plt

        for i in eigs:
            # for i in eigs.T:
            plt.plot(i)
        plt.savefig(band_file)
        plt.close()
    return eigs


def read_spaghetti_ene(
    filename="ICSD-76748.spaghetti_ene",
):
    """Obtain data for plotting bandstructure."""
    f = open(filename, "r")
    lines = f.read().splitlines()
    f.close()
    k = []
    energy = []
    for i, ii in enumerate(lines):
        if "bandindex" not in ii:
            tmp = [float(j) for j in ii.split()]
            energy.append(tmp[-1])
            k.append(tmp[-2])
    k = np.array(k)
    energy = np.array(energy)
    return k, energy


"""
read_band_energy()
x = read_scf()
print(x)
"""
