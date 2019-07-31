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


def ip_optics(vrun="", filename=".", extra="", plot=False):
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


if __name__ == "__main__":
    run = path = str(
        os.path.join(
            os.path.dirname(__file__),
            "../vasp/examples/SiOptb88/MAIN-MBJ-bulk@mp_149/vasprun.xml",
        )
    )
    ep = ip_optics(vrun=run)
    print('val',ep["absorption"][0][0])
