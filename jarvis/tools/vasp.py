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
