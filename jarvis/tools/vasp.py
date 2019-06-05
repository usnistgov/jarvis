from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.io.vasp.outputs import Oszicar, Outcar, Vasprun
from pymatgen.util.plotting import pretty_plot as get_publication_quality_plot

# from pymatgen.util.plotting_utils import get_publication_quality_plot
import matplotlib, yaml, os
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Potcar, Incar, Kpoints
from pymatgen.electronic_structure.plotter import (
    BSPlotterProjected,
    BSPlotter,
    DosPlotter,
)
import matplotlib

matplotlib.use("Agg")
import math, glob
from pymatgen.core.structure import Structure
import numpy as np
from pymatgen.analysis.elasticity.elastic import ElasticTensor


def bandgap(vrun="", occu_tol=0.001):
    """
    Parses vasprun.xml and gives bandgap 

    Args:
        vrun: path to vasprun.xml 
        occu_tol: occupation tolerance
    Returns:
           eg: bandgap
           direct: direct/indirect bandgap
    """
    v = Vasprun(vrun, occu_tol=occu_tol)
    eg = float(v.eigenvalue_band_properties[0])
    direct = v.eigenvalue_band_properties[3]

    return eg, direct


def bandstr(vrun="", kpfile="", filename=".", plot=False):
    """
    Plot electronic bandstructure
 
    Args:
        vrun: path to vasprun.xml
        kpfile:path to line mode KPOINTS file 
    Returns:
           matplotlib object
    """

    run = Vasprun(vrun, parse_projected_eigen=True)
    bands = run.get_band_structure(kpfile, line_mode=True, efermi=run.efermi)
    bsp = BSPlotter(bands)
    zero_to_efermi = True
    bandgap = str(round(bands.get_band_gap()["energy"], 3))
    # print "bg=",bandgap
    data = bsp.bs_plot_data(zero_to_efermi)
    plt = get_publication_quality_plot(12, 8)
    plt.close()
    plt.clf()
    band_linewidth = 3
    x_max = data["distances"][-1][-1]
    # print (x_max)
    for d in range(len(data["distances"])):
        for i in range(bsp._nb_bands):
            plt.plot(
                data["distances"][d],
                [
                    data["energy"][d]["1"][i][j]
                    for j in range(len(data["distances"][d]))
                ],
                "b-",
                linewidth=band_linewidth,
            )
            if bsp._bs.is_spin_polarized:
                plt.plot(
                    data["distances"][d],
                    [
                        data["energy"][d]["-1"][i][j]
                        for j in range(len(data["distances"][d]))
                    ],
                    "r--",
                    linewidth=band_linewidth,
                )
    bsp._maketicks(plt)
    if bsp._bs.is_metal():
        e_min = -10
        e_max = 10
        band_linewidth = 3

    for cbm in data["cbm"]:
        plt.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)

        for vbm in data["vbm"]:
            plt.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

    plt.xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
    ylabel = (
        r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_to_efermi else r"$\mathrm{Energy\ (eV)}$"
    )
    plt.ylabel(ylabel, fontsize=30)
    plt.ylim(-4, 4)
    plt.xlim(0, x_max)
    plt.tight_layout()
    if plot == True:
        plt.savefig(filename, img_format="png")
        plt.close()

    return plt


def plot_dos(vrun="", low_lim=-5, up_lim=10, filename=".", extra="", plot=False):
    """
    Plot total density of states

    Args:
       vrun: vasprun.xml path
       low_lim: lower energy limit in eV
       up_lim: upper energy limit in eV
    Returns:
          matplotlib objects for
          total, orbital and element project density of states
    """

    run = Vasprun(vrun, occu_tol=1e-2)
    complete_dos = run.complete_dos
    totup = complete_dos.get_densities(spin=Spin.up)
    totdn = complete_dos.get_densities(spin=Spin.down)
    en = complete_dos.energies
    bandgap = str(run.eigenvalue_band_properties[0])
    ef = complete_dos.efermi
    en[:] = [x - ef for x in en]

    plt1 = get_publication_quality_plot(14, 10)
    # plt1.close()
    # plt1.clf()
    plt1.xlabel("Energies (eV)")
    plt1.ylabel("Tot. DOS (arb. unit.)")
    plt1.plot(en, totup, "b", linewidth=3)
    plt1.plot(en, -totdn, "r", linewidth=3)
    plt1.xlim(low_lim, up_lim)
    plt1.tight_layout()
    if plot == True:
        tmp = str(filename) + str("/TDos") + str(extra) + str(".png")
        plt1.savefig(tmp)

    spd_dos = complete_dos.get_spd_dos()
    sdosup = spd_dos[OrbitalType.s].densities[Spin.up]
    pdosup = spd_dos[OrbitalType.p].densities[Spin.up]
    ddosup = spd_dos[OrbitalType.d].densities[Spin.up]
    sdosdn = spd_dos[OrbitalType.s].densities[Spin.down]
    pdosdn = spd_dos[OrbitalType.p].densities[Spin.down]
    ddosdn = spd_dos[OrbitalType.d].densities[Spin.down]

    plt2 = get_publication_quality_plot(14, 10)
    # plt2.close()
    # plt2.clf()
    plt2.plot(en, sdosup, "r", linewidth=3, label="s")
    plt2.plot(en, -sdosdn, "r", linewidth=3)
    plt2.plot(en, pdosup, "g", linewidth=3, label="p")
    plt2.plot(en, -pdosdn, "g", linewidth=3)
    plt2.plot(en, ddosup, "b", linewidth=3, label="d")
    plt2.plot(en, -ddosdn, "b", linewidth=3)
    plt2.xlabel("Energies (eV)")
    plt2.ylabel("Orb. DOS (arb. unit.)")
    plt2.legend(prop={"size": 26})
    plt2.xlim(low_lim, up_lim)
    plt2.tight_layout()
    if plot == True:
        tmp = str(filename) + str("/Dos") + str(extra) + str(".png")
        plt2.savefig(tmp)

    plt3 = get_publication_quality_plot(14, 10)
    # plt3.close()
    # plt3.clf()
    elt_dos = complete_dos.get_element_dos()
    at_symbs = set(run.atomic_symbols)
    for i, sym in enumerate(at_symbs):
        elt = elt_dos[Element(sym)]
        # color = cmap(float(i)/len(at_symbs))
        if i == 0:
            color = "b"
        elif i == 1:
            color = "g"
        elif i == 2:
            color = "r"
        elif i == 3:
            color = "c"
        elif i == 4:
            color = "m"
        elif i == 5:
            color = "y"
        elif i == 6:
            color = "k"
        else:
            color = "k"
            print("Use different color scheme")

        plt3.plot(en, elt.densities[Spin.up], color=color, linewidth=3, label=sym)
        plt3.plot(en, -elt.densities[Spin.down], color=color, linewidth=3)
        plt3.xlabel("Energies (eV)")
        plt3.ylabel("Elm. DOS (arb. unit.)")
        plt3.legend(prop={"size": 26})
        plt3.xlim(low_lim, up_lim)
        plt3.tight_layout()
    if plot == True:
        tmp = str(filename) + str("/EDos") + str(extra) + str(".png")
        plt3.savefig(tmp)

    return plt1, plt2, plt3


def plot_enc_convergence(
    directory="../vasp/examples/SiOptb88/", plot=False, filename="."
):
    """
    Plot convergence for plane-wave cut-off data
    Works only if jobs run through jarvis-tools framework
    
    Args:
       directory: parent directory for job run
    Returns:
           matplotlib object, converged cut-off value
    """

    x = []
    y = []
    for a in glob.glob(str(directory) + str("/*.json")):
        if "MAIN-RELAX" in a:
            main_inc = str(a.split(".json")[0]) + str("/") + str("INCAR")
            main_inc_obj = Incar.from_file(main_inc)
            convg_encut = float(main_inc_obj["ENCUT"])
        elif "ENCUT" in a:
            run = str(a.split(".json")[0]) + str("/") + str("vasprun.xml")
            contcar = Structure.from_file(
                (a.split(".json")[0]) + str("/") + str("CONTCAR")
            )
            vrun = Vasprun(run)
            infile = str(a.split(".json")[0]) + str("/") + str("INCAR")
            inc = Incar.from_file(infile)
            encut = inc["ENCUT"]
            en = float(vrun.final_energy)  # /float(contcar.composition.num_atoms)
            x.append(encut)
            y.append(en)
    order = np.argsort(x)
    xs = np.array(x)[order]
    ys = np.array(y)[order]
    plt = get_publication_quality_plot(14, 10)
    plt.ylabel("Energy (eV)")
    plt.plot(xs, ys, "s-", linewidth=2, markersize=10)
    plt.xlabel("Increment in ENCUT ")
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    plt.title(str("Converged at ") + str(int(convg_encut)) + str("eV"), fontsize=26)
    # filename=str('Encut.png')
    plt.tight_layout()
    if plot == True:
        plt.savefig(filename)
        plt.close()

    return plt, convg_encut


def plot_kp_convergence(
    directory="../vasp/examples/SiOptb88/", plot=False, filename="."
):
    """
    Plot convergence for k-points data
    Works only if jobs run through jarvis-tools framework
    
    Args:
       directory: parent directory for job run
    Returns:
           matplotlib object, converged k-points value
    """

    x = []
    y = []
    for a in glob.glob(str(directory) + str("/*.json")):
        if "MAIN-RELAX" in a:
            main_kp = str(a.split(".json")[0]) + str("/") + str("KPOINTS")
            main_kp_obj = Kpoints.from_file(main_kp)
            [k1, k2, k3] = main_kp_obj.kpts[0]
        elif "KPOINT" in a:
            k = int(float(str(a.split("-")[-1]).split(".json")[0]))
            run = str(a.split(".json")[0]) + str("/") + str("vasprun.xml")
            vrun = Vasprun(run)

            kpfile = str(a.split(".json")[0]) + str("/") + str("KPOINTS")
            contcar = Structure.from_file(
                (a.split(".json")[0]) + str("/") + str("CONTCAR")
            )
            kpoints = Kpoints.from_file(kpfile)
            [xx, yy, zz] = kpoints.kpts[0]
            en = float(vrun.final_energy)  # /float(contcar.composition.num_atoms)
            # en =float(vrun.final_energy)/float(contcar.composition.num_atoms)
            # print "ENERGYYY,AT",en,float(contcar.composition.num_atoms)
            x.append(k)
            y.append(en)

    order = np.argsort(x)
    xs = np.array(x)[order]
    ys = np.array(y)[order]
    len_xs = len(xs)
    xs1 = []
    ys1 = []
    target_ys = ys[-1]
    # print "target=",target_ys
    for i, el in enumerate(ys):
        if el <= (float(target_ys) + 0.002):
            xs1.append(xs[i])
            ys1.append(ys[i])
    # print "xs,ys=",xs,ys
    # print "xs1,ys1=",xs1,ys1
    left, bottom, width, height = [0.5, 0.5, 0.35, 0.35]
    plt = get_publication_quality_plot(14, 10)
    fig, ax1 = plt.subplots()
    plt.xlabel("Increment in K point", fontsize=20)
    plt.ylabel("Energy (eV)", fontsize=20)
    plt.title(
        str("Converged at ")
        + str(k1)
        + str("x")
        + str(k2)
        + str("x")
        + str(k3)
        + str(" ")
        + str("Automatic Mesh ")
        + str(max(x) - 25),
        fontsize=14,
    )
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ax2 = fig.add_axes([left, bottom, width, height])
    # ax2.plot(xs1,ys1, '.-',linewidth=2,markersize=10)
    ax1.plot(xs, ys, "s-", linewidth=2, markersize=10)
    el_list = sorted(list(contcar.symbol_set))
    search = "-".join([item for item in el_list])
    # ax1.xlabel('Increment in K point')
    # print "xs,ys"
    # print xs
    # print ys
    plt.plot(xs1, ys1, ".-", linewidth=2, markersize=10)
    plt.ylim([float(target_ys) + 0.002, float(target_ys) - 0.002])
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    filename = str("KDen.png")
    k_convg = str(xx) + str("x") + str(yy) + str("x") + str(zz)
    plt.tight_layout()
    if plot == True:
        plt.savefig(filename)
        plt.close()
    return plt, k_convg


def IP_optics(vrun="", filename=".", extra="", plot=False):
    """
    Optoelectronic properties derived from dielectric function 

    See http://www.wien2k.at/reg_user/textbooks/WIEN2k_lecture-notes_2013/optic_handout.pdf
    http://www.phys.ufl.edu/~tanner/Wooten-OpticalPropertiesOfSolids-2up.pdf

    Args:
       vrun: path to vasprun.xml
    Reurns:
          info: dictionary with optoelectronic information
          energy: energy bins generally 0 to 40 eV
          real_part: real part of dielectric function
          imag_part: imaginary part of dielectric function
          absorption: absorption coefficient in cm^-1
       
    TODO: shorten the code with multi-dimensional np array   
    """

    info = {}

    storedir = filename

    run = Vasprun(vrun, occu_tol=1e-2)
    erange = len(run.dielectric[0])
    en = []
    realx = []
    imagx = []
    absorpx = []
    refrx = []
    reflx = []
    eelsx = []
    extcx = []
    opt_conx = []
    realy = []
    imagy = []
    absorpy = []
    refry = []
    refly = []
    eelsy = []
    extcy = []
    opt_cony = []
    realz = []
    imagz = []
    absorpz = []
    refrz = []
    reflz = []
    eelsz = []
    extcz = []
    opt_conz = []
    H = 4.13566733 * (10 ** (-15))
    # c0=2.99792458
    c0 = 2.99792458 * (math.pow(10, 8))
    for i in range(0, erange - 1):
        en.append(run.dielectric[0][i])
        realx.append(run.dielectric[1][i][0])
        imagx.append(run.dielectric[2][i][0])
        ab_valx = 1.4142 * (
            (float(run.dielectric[0][i]) / float(H))
            * (
                float(
                    math.sqrt(
                        -run.dielectric[1][i][0]
                        + math.sqrt(
                            (run.dielectric[2][i][0]) * (run.dielectric[2][i][0])
                            + (run.dielectric[1][i][0]) * (run.dielectric[1][i][0])
                        )
                    )
                )
                / float(float(c0) * 100.0)
            )
        )
        absorpx.append(ab_valx)
        refr_valx = float(
            math.sqrt(
                (run.dielectric[1][i][0])
                + math.sqrt(
                    (run.dielectric[2][i][0]) * (run.dielectric[2][i][0])
                    + (run.dielectric[1][i][0]) * (run.dielectric[1][i][0])
                )
            )
        ) / float(1.4142)
        refrx.append(refr_valx)
        eels_valx = float(run.dielectric[2][i][0]) / float(
            (run.dielectric[2][i][0]) * (run.dielectric[2][i][0])
            + (run.dielectric[1][i][0]) * (run.dielectric[1][i][0])
        )
        eelsx.append(eels_valx)
        extc_valx = float(
            math.sqrt(
                -(run.dielectric[1][i][0])
                + math.sqrt(
                    (run.dielectric[2][i][0]) * (run.dielectric[2][i][0])
                    + (run.dielectric[1][i][0]) * (run.dielectric[1][i][0])
                )
            )
        ) / float(1.4142)
        extcx.append(extc_valx)
        refl_valx = ((refr_valx - 1) * (refr_valx - 1) + extc_valx * extc_valx) / (
            (refr_valx + 1) * (refr_valx + 1) + extc_valx * extc_valx
        )
        reflx.append(refl_valx)
        opt_valx = (
            float(float(run.dielectric[0][i]) / float(H))
            / float(4 * math.pi)
            * (run.dielectric[2][i][0])
        )
        opt_conx.append(opt_valx)

        realy.append(run.dielectric[1][i][1])
        imagy.append(run.dielectric[2][i][1])

        ab_valy = 1.4142 * (
            (float(run.dielectric[0][i]) / float(H))
            * (
                float(
                    math.sqrt(
                        -run.dielectric[1][i][1]
                        + math.sqrt(
                            (run.dielectric[2][i][1]) * (run.dielectric[2][i][1])
                            + (run.dielectric[1][i][1]) * (run.dielectric[1][i][1])
                        )
                    )
                )
                / float(float(c0) * 100.0)
            )
        )
        absorpy.append(ab_valy)
        refr_valy = float(
            math.sqrt(
                (run.dielectric[1][i][1])
                + math.sqrt(
                    (run.dielectric[2][i][1]) * (run.dielectric[2][i][1])
                    + (run.dielectric[1][i][1]) * (run.dielectric[1][i][1])
                )
            )
        ) / float(1.4142)
        refry.append(refr_valy)
        eels_valy = float(run.dielectric[2][i][1]) / float(
            (run.dielectric[2][i][1]) * (run.dielectric[2][i][1])
            + (run.dielectric[1][i][1]) * (run.dielectric[1][i][1])
        )
        eelsy.append(eels_valy)
        extc_valy = float(
            math.sqrt(
                -(run.dielectric[1][i][1])
                + math.sqrt(
                    (run.dielectric[2][i][1]) * (run.dielectric[2][i][1])
                    + (run.dielectric[1][i][1]) * (run.dielectric[1][i][1])
                )
            )
        ) / float(1.4142)
        extcy.append(extc_valy)
        refl_valy = ((refr_valy - 1) * (refr_valy - 1) + extc_valy * extc_valy) / (
            (refr_valy + 1) * (refr_valy + 1) + extc_valy * extc_valy
        )
        refly.append(refl_valy)
        opt_valy = (
            float(float(run.dielectric[0][i]) / float(H))
            / float(4 * math.pi)
            * (run.dielectric[2][i][1])
        )
        opt_cony.append(opt_valy)

        realz.append(run.dielectric[1][i][2])
        imagz.append(run.dielectric[2][i][2])

        ab_valz = 1.4142 * (
            (float(run.dielectric[0][i]) / float(H))
            * (
                float(
                    math.sqrt(
                        -run.dielectric[1][i][2]
                        + math.sqrt(
                            (run.dielectric[2][i][2]) * (run.dielectric[2][i][2])
                            + (run.dielectric[1][i][2]) * (run.dielectric[1][i][2])
                        )
                    )
                )
                / float(float(c0) * 100.0)
            )
        )
        absorpz.append(ab_valz)
        refr_valz = float(
            math.sqrt(
                (run.dielectric[1][i][2])
                + math.sqrt(
                    (run.dielectric[2][i][2]) * (run.dielectric[2][i][2])
                    + (run.dielectric[1][i][2]) * (run.dielectric[1][i][2])
                )
            )
        ) / float(1.4142)
        refrz.append(refr_valz)
        eels_valz = float(run.dielectric[2][i][2]) / float(
            (run.dielectric[2][i][2]) * (run.dielectric[2][i][2])
            + (run.dielectric[1][i][2]) * (run.dielectric[1][i][2])
        )
        eelsz.append(eels_valz)
        extc_valz = float(
            math.sqrt(
                -(run.dielectric[1][i][2])
                + math.sqrt(
                    (run.dielectric[2][i][2]) * (run.dielectric[2][i][2])
                    + (run.dielectric[1][i][2]) * (run.dielectric[1][i][2])
                )
            )
        ) / float(1.4142)
        extcz.append(extc_valz)
        refl_valz = ((refr_valz - 1) * (refr_valz - 1) + extc_valz * extc_valz) / (
            (refr_valz + 1) * (refr_valz + 1) + extc_valz * extc_valz
        )
        reflz.append(refl_valz)
        opt_valz = (
            float(float(run.dielectric[0][i]) / float(H))
            / float(4 * math.pi)
            * (run.dielectric[2][i][2])
        )
        opt_conz.append(opt_valz)

    if plot == True:
        # plt.close()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, realx, linewidth=2, label=r"$\epsilon_1x$")
        plt.plot(en, realy, linewidth=2, label=r"$\epsilon_1y$")
        plt.plot(en, realz, linewidth=2, label=r"$\epsilon_1z$")
        plt.legend(prop={"size": 26})
        plt.xlim([0, 30])
        plt.xlabel("Energy (eV)")
        plt.ylabel("Dielec. function (Real part)")
        filename = str(storedir) + str("/") + str("Real") + str(extra) + str(".png")
        plt.tight_layout()
        plt.savefig(filename)
        plt.tight_layout()
        plt.close()

        plt.clf()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, imagx, linewidth=2, label=r"$\epsilon_2x$")
        plt.plot(en, imagy, linewidth=2, label=r"$\epsilon_2y$")
        plt.plot(en, imagz, linewidth=2, label=r"$\epsilon_2z$")
        plt.xlabel("Energy (eV)")
        plt.xlim([0, 30])
        plt.ylabel("Dielec. function(Imag part)")
        filename = str(storedir) + str("/") + str("Imag") + str(extra) + str(".png")
        plt.tight_layout()
        plt.legend(prop={"size": 26})
        plt.savefig(filename)
        plt.close()

        plt.clf()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, refrx, linewidth=2, label=r"$n_x$")
        plt.plot(en, refry, linewidth=2, label=r"$n_y$")
        plt.plot(en, refrz, linewidth=2, label=r"$n_z$")
        plt.xlim([0, 30])
        plt.xlabel("Energy (eV)")
        plt.ylabel("Refractive index")
        filename = str(storedir) + str("/") + str("Refr") + str(extra) + str(".png")
        plt.tight_layout()
        plt.legend(prop={"size": 26})
        plt.savefig(filename)
        plt.close()

        plt.clf()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, extcx, linewidth=2, label=r"$k_x$")
        plt.plot(en, extcy, linewidth=2, label=r"$k_y$")
        plt.plot(en, extcz, linewidth=2, label=r"$k_z$")
        plt.xlim([0, 30])
        plt.xlabel("Energy (eV)")
        plt.ylabel("Extinction Coefficient")
        filename = str(storedir) + str("/") + str("Extc") + str(extra) + str(".png")
        plt.tight_layout()
        plt.legend(prop={"size": 26})
        plt.savefig(filename)
        plt.tight_layout()
        plt.close()

        plt.clf()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, absorpx, linewidth=2, label=r"$\alpha_x$")
        plt.plot(en, absorpy, linewidth=2, label=r"$\alpha_y$")
        plt.plot(en, absorpz, linewidth=2, label=r"$\alpha_z$")
        plt.xlim([0, 30])
        plt.xlabel("Energy (eV)")
        plt.ylabel("Absorption coefficient")
        filename = str(storedir) + str("/") + str("Absorp") + str(extra) + str(".png")
        plt.tight_layout()
        plt.legend(prop={"size": 26})
        plt.savefig(filename)
        plt.tight_layout()
        plt.close()

        plt.clf()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, eelsx, linewidth=2, label=r"$e_x$")
        plt.plot(en, eelsy, linewidth=2, label=r"$e_y$")
        plt.plot(en, eelsz, linewidth=2, label=r"$e_z$")
        plt.xlim([0, 30])
        plt.xlabel("Energy (eV)")
        plt.ylabel("Energy loss spectrum")
        filename = str(storedir) + str("/") + str("ELS") + str(extra) + str(".png")
        plt.tight_layout()
        plt.legend(prop={"size": 26})
        plt.savefig(filename)
        plt.tight_layout()
        plt.close()

        plt.clf()
        plt = get_publication_quality_plot(14, 10)
        plt.plot(en, opt_conx, linewidth=2, label=r"$\sigma_x$")
        plt.plot(en, opt_cony, linewidth=2, label=r"$\sigma_y$")
        plt.plot(en, opt_conz, linewidth=2, label=r"$\sigma_z$")
        plt.xlim([0, 30])
        plt.xlabel("Energy (eV)")
        plt.ylabel("Optical conductivity (Real part)")
        filename = str(storedir) + str("/") + str("Opt_con") + str(extra) + str(".png")
        plt.tight_layout()
        plt.legend(prop={"size": 26})
        plt.savefig(filename)
        plt.tight_layout()
        plt.close()

    info["absorption"] = [absorpx, absorpy, absorpz]
    info["energy"] = en
    info["real_part"] = [realx, realy, realz]
    info["imag_part"] = [imagx, imagy, imagz]
    info["refr_index"] = [refrx, refry, refrz]
    info["extinc"] = [extcx, extcy, extcz]
    info["reflectance"] = [reflx, refly, reflz]
    info["EELS"] = [eelsx, eelsy, eelsz]
    info["opt_conduct"] = [opt_conx, opt_cony, opt_conz]
    print(info["real_part"][0][0])
    return info


def elastic_props(outcar="", vacuum=False):
    """
    Obtain elastic tensor and calculate related properties

    Args:
        outcar: OUTCAR file path
        vacuum: whether the structure has vaccum such as 2D materials
        for vacuum structures bulk and shear mod. needs extra attenstion
        and elastic tensor are in Nm^-1 rather than GPa
    Returns:
          info: data for elastic tensor (in string and object representation), bulk, shear modulus, and phonon modes
    """
    if vacuum == True:
        contcar = Structure.from_file(str(outcar).replace("OUTCAR", "POSCAR"))
        ratio_c = 0.1 * float(
            abs(contcar.lattice.matrix[2][2])
        )  # *(10**9)*(10**-10) #N/m unit
    el_tens = "na"
    KV = "na"
    GV = "na"
    spin = "na"
    info = {}
    ratio_c = 1.0
    v = open(outcar, "r")
    lines = v.read().splitlines()
    for i, line in enumerate(lines):
        if "TOTAL ELASTIC MODULI (kBar)" in line:
            c11 = lines[i + 3].split()[1]
            c12 = lines[i + 3].split()[2]
            c13 = lines[i + 3].split()[3]
            c14 = lines[i + 3].split()[4]
            c15 = lines[i + 3].split()[5]
            c16 = lines[i + 3].split()[6]
            c21 = lines[i + 4].split()[1]
            c22 = lines[i + 4].split()[2]
            c23 = lines[i + 4].split()[3]
            c24 = lines[i + 4].split()[4]
            c25 = lines[i + 4].split()[5]
            c26 = lines[i + 4].split()[6]
            c31 = lines[i + 5].split()[1]
            c32 = lines[i + 5].split()[2]
            c33 = lines[i + 5].split()[3]
            c34 = lines[i + 5].split()[4]
            c35 = lines[i + 5].split()[5]
            c36 = lines[i + 5].split()[6]
            c41 = lines[i + 6].split()[1]
            c42 = lines[i + 6].split()[2]
            c43 = lines[i + 6].split()[3]
            c44 = lines[i + 6].split()[4]
            c45 = lines[i + 6].split()[5]
            c46 = lines[i + 6].split()[6]
            c51 = lines[i + 7].split()[1]
            c52 = lines[i + 7].split()[2]
            c53 = lines[i + 7].split()[3]
            c54 = lines[i + 7].split()[4]
            c55 = lines[i + 7].split()[5]
            c56 = lines[i + 7].split()[6]
            c61 = lines[i + 8].split()[1]
            c62 = lines[i + 8].split()[2]
            c63 = lines[i + 8].split()[3]
            c64 = lines[i + 8].split()[4]
            c65 = lines[i + 8].split()[5]
            c66 = lines[i + 8].split()[6]
            c11 = round(ratio_c * float(c11) / float(10), 1)
            c12 = round(ratio_c * float(c12) / float(10), 1)
            c13 = round(ratio_c * float(c13) / float(10), 1)
            c14 = round(ratio_c * float(c14) / float(10), 1)
            c15 = round(ratio_c * float(c15) / float(10), 1)
            c16 = round(ratio_c * float(c16) / float(10), 1)
            c21 = round(ratio_c * float(c21) / float(10), 1)
            c22 = round(ratio_c * float(c22) / float(10), 1)
            c23 = round(ratio_c * float(c23) / float(10), 1)
            c24 = round(ratio_c * float(c24) / float(10), 1)
            c25 = round(ratio_c * float(c25) / float(10), 1)
            c26 = round(ratio_c * float(c26) / float(10), 1)
            c31 = round(float(c31) / float(10), 1)
            c32 = round(float(c32) / float(10), 1)
            c33 = round(float(c33) / float(10), 1)
            c34 = round(float(c34) / float(10), 1)
            c35 = round(float(c35) / float(10), 1)
            c36 = round(float(c36) / float(10), 1)
            c41 = round(float(c41) / float(10), 1)
            c42 = round(float(c42) / float(10), 1)
            c43 = round(float(c43) / float(10), 1)
            c44 = round(float(c44) / float(10), 1)
            c45 = round(float(c45) / float(10), 1)
            c46 = round(float(c46) / float(10), 1)
            c51 = round(float(c51) / float(10), 1)
            c52 = round(float(c52) / float(10), 1)
            c53 = round(float(c53) / float(10), 1)
            c54 = round(float(c54) / float(10), 1)
            c55 = round(float(c55) / float(10), 1)
            c56 = round(float(c56) / float(10), 1)
            c61 = round(float(c61) / float(10), 1)
            c62 = round(float(c62) / float(10), 1)
            c63 = round(float(c63) / float(10), 1)
            c64 = round(float(c64) / float(10), 1)
            c65 = round(float(c65) / float(10), 1)
            c66 = round(float(c66) / float(10), 1)
            KV = float((c11 + c22 + c33) + 2 * (c12 + c23 + c31)) / float(9)
            GV = float(
                (c11 + c22 + c33) - (c12 + c23 + c31) + 3 * (c44 + c55 + c66)
            ) / float(15)
            KV = round(KV, 3)
            GV = round(GV, 3)
            break
    v.close()

    # Convenient string representation for storage

    el_tens = (
        str(c11)
        + str(",")
        + str(c12)
        + str(",")
        + str(c13)
        + str(",")
        + str(c14)
        + str(",")
        + str(c15)
        + str(",")
        + str(c16)
        + str(",")
        + str(c21)
        + str(",")
        + str(c22)
        + str(",")
        + str(c23)
        + str(",")
        + str(c24)
        + str(",")
        + str(c25)
        + str(",")
        + str(c26)
        + str(",")
        + str(c31)
        + str(",")
        + str(c32)
        + str(",")
        + str(c33)
        + str(",")
        + str(c34)
        + str(",")
        + str(c35)
        + str(",")
        + str(c36)
        + str(",")
        + str(c41)
        + str(",")
        + str(c42)
        + str(",")
        + str(c43)
        + str(",")
        + str(c44)
        + str(",")
        + str(c45)
        + str(",")
        + str(c46)
        + str(",")
        + str(c51)
        + str(",")
        + str(c52)
        + str(",")
        + str(c53)
        + str(",")
        + str(c54)
        + str(",")
        + str(c55)
        + str(",")
        + str(c56)
        + str(",")
        + str(c61)
        + str(",")
        + str(c62)
        + str(",")
        + str(c63)
        + str(",")
        + str(c64)
        + str(",")
        + str(c65)
        + str(",")
        + str(c66)
    )

    cij = np.empty((6, 6), dtype=float)
    elast = np.array(el_tens.split(","), dtype="float")
    count = 0
    for ii in range(6):
        for jj in range(6):
            cij[ii][jj] = elast[count]
    el_tens2 = ElasticTensor.from_voigt(cij)

    modes = []
    try:
        for i in lines:
            if "cm-1" in i and "meV" in i:
                mod = float(i.split()[7])
                if mod not in modes:
                    modes.append(float(mod))
    except:
        pass

    info["el_tens_str"] = el_tens
    info["el_tens_obj"] = el_tens2
    info["KV"] = KV
    info["GV"] = GV
    info["modes"] = modes

    return info


def magnetic_moment(oszicar=""):
    """
    Obtain total orbital magnetic moment
    Magnetic moment info is found in both OSZICAR and OUTCAR
    We prefer OSZICAR values

    Args:
        oszicar: path to OSZICAR
    Returns:
          magmom: orbital magnetic moment in Bohr-magneton
    """
    magmom = Oszicar(oszicar).ionic_steps[-1]["mag"]
    # magmom_out=out.total_mag #for outcar
    return magmom


def get_area(strt=""):
    """
    Get area of a Structure object along z-direction

    Args: 
        strt: Structure object
    Returns:
           area
    """

    m = strt.lattice.matrix
    area = np.linalg.norm(np.cross(m[0], m[1]))
    return area


def get_spacegroup(strt=""):
    """
    Get spacegroup of a Structure pbject
    
    Args:
        strt: Structure object
    Returns:
           num: spacegroup number
           symb: spacegroup symbol
    """
    finder = SpacegroupAnalyzer(strt)
    num = finder.get_space_group_number()
    symb = finder.get_space_group_symbol()
    return num, symb


# vrun='../vasp/examples/SiOptb88/MAIN-OPTICS-bulk@mp_149/vasprun.xml'
# out='../vasp/examples/SiOptb88/MAIN-ELASTIC-bulk@mp_149/OUTCAR'
# osz='../vasp/examples/SiOptb88/MAIN-RELAX-bulk@mp_149/OSZICAR'
# m=magnetic_moment(osz)
# print (m)
# info=elastic_props(out)
# print (info['KV'])
# op=IP_optics(vrun,plot=True)
##op['real_part'][0][0]
# print (op)
# plot_kp_convergence()
