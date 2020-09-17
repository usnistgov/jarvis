"""Module to make XML for JARVIS-REST-API schema."""
import os
import glob
from jarvis.core.atoms import Atoms
import yaml
from jarvis.core.atoms import get_supercell_dims
from jarvis.analysis.thermodynamics.energetics import get_twod_defect_energy
from jarvis.io.vasp.outputs import Oszicar, Vasprun, Outcar
from matplotlib.pyplot import imread
from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.analysis.topological.spillage import Spillage
from jarvis.db.jsonutils import loadjson
from jarvis.analysis.phonon.ir import ir_intensity
from jarvis.analysis.structure.neighbors import NeighborsAnalysis
from jarvis.analysis.solarefficiency.solar import SolarEfficiency
from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
from jarvis.analysis.thermodynamics.energetics import form_enp
from jarvis.tasks.boltztrap.run import run_boltztrap
from jarvis.tasks.phonopy.run import run_phonopy
# from jarvis.io.phonopy.outputs import bandstructure_plot
from jarvis.ai.pkgs.utils import get_ml_data
from jarvis.analysis.diffraction.xrd import XRD
import numpy as np
from jarvis.analysis.elastic.tensor import ElasticTensor
from jarvis.io.boltztrap.outputs import BoltzTrapOutput
from jarvis.io.vasp.outputs import parse_raman_dat
from jarvis.core.utils import array_to_string
from jarvis.core.utils import stringdict_to_xml
from jarvis.ai.descriptors.cfid import CFID
from jarvis.io.wannier.outputs import WannierHam


cfid_x, cfid_y, cfid_jids = get_ml_data()


def get_cfid_descriptors(jid="JVASP-1002", atoms="", make_cfid=True):
    """Get CFID pre-computed descriptors for a JID."""
    for i, j in zip(cfid_x, cfid_jids):
        if j == jid:
            return i
    if make_cfid:
        print("Calculating CFID descriptors.")
        cfid = CFID(atoms)
        return cfid.get_comp_descp().tolist()


class VaspToApiXmlSchema(object):
    """Module to convert VASP data to XML schema for API."""

    def __init__(
        self,
        folder="",
        meta_data={
            "id_file": "JARVIS-ID",
            "ehull": "",
            "data_source": "JARVIS-DFT-VASP",
            "material_type": "bulk",
            "exfoliation_energy": "",
        },
    ):
        """Initialize class."""
        self.folder = folder
        if "1L" in folder:
            meta_data["material_type"] = "SingleLayer"
        if "2L" in folder:
            meta_data["material_type"] = "BiLayer"
        if "3L" in folder:
            meta_data["material_type"] = "TriLayer"
        self.meta_data = meta_data

    def ebandstruct(self, vrun="", kp=""):
        """Get electronic bandstucre related data."""
        line = ""
        try:
            f = open(kp, "r")
            lines = f.read().splitlines()
            f.close()

            vrun = Vasprun(vrun)
            band_indir_gap = ""
            band_dir_gap = ""
            try:
                band_indir_gap = str(round(vrun.get_indir_gap[0], 2))
                band_dir_gap = str(round(vrun.get_dir_gap, 2))
            except Exception:
                print("Cannot get bandgap.", kp)

            line += "<band_indir_gap>" + band_indir_gap + "</band_indir_gap>"
            line += "<band_dir_gap>" + band_dir_gap + "</band_dir_gap>"

            fermi_velocities = ""
            try:
                fermi_velocities = ",".join(map(str, vrun.fermi_velocities[0]))
            except Exception:
                print("Cannot get the Fermi-velocity.", kp)
                pass
            line += (
                '<fermi_velocities>"'
                + fermi_velocities
                + '"</fermi_velocities>'
            )

            kp_labels = []
            kp_labels_points = []
            for ii, i in enumerate(lines):
                if ii > 2:
                    tmp = i.split()
                    if len(tmp) == 5:
                        tmp = str(tmp[4])
                        if len(kp_labels) == 0:
                            kp_labels.append(tmp)
                            kp_labels_points.append(ii - 3)
                        elif tmp != kp_labels[-1]:
                            kp_labels.append(tmp)
                            kp_labels_points.append(ii - 3)
            tmp = 0.0
            zero_efermi = True
            if zero_efermi:
                tmp = float(vrun.efermi)
            spin = 0
            up_bands_y = []
            up_bands_x = []
            for i, ii in enumerate(vrun.eigenvalues[spin][:, :, 0].T - tmp):
                y = ",".join(map(str, ii))
                x = ",".join(map(str, range(0, len(ii))))
                up_bands_y.append(y)
                up_bands_x.append(x)
            spin = 1
            down_bands_y = []
            down_bands_x = []
            for i, ii in enumerate(vrun.eigenvalues[spin][:, :, 0].T - tmp):
                y = ",".join(map(str, ii))
                x = ",".join(map(str, range(0, len(ii))))
                down_bands_y.append(y)
                down_bands_x.append(x)
            line += (
                '<kp_labels_points>"'
                + ",".join(map(str, kp_labels_points))
                + '"</kp_labels_points>'
            )
            line += (
                '<kp_labels>"'
                + ",".join(map(str, kp_labels))
                + '"</kp_labels>'
            )
            line += (
                '<spin_up_bands_x>"'
                + ";".join(map(str, up_bands_x))
                + '"</spin_up_bands_x>'
            )
            line += (
                '<spin_up_bands_y>"'
                + ";".join(map(str, up_bands_y))
                + '"</spin_up_bands_y>'
            )
            line += (
                '<spin_down_bands_x>"'
                + ";".join(map(str, down_bands_x))
                + '"</spin_down_bands_x>'
            )
            line += (
                '<spin_down_bands_y>"'
                + ";".join(map(str, down_bands_y))
                + '"</spin_down_bands_y>'
            )
        except Exception:
            print("Could not get bandstructure info,vrun,kp", vrun, kp)
            pass

        return line

    def encut_kp(self):
        """Get ENCUT, cut-off convergence data."""
        folder = self.folder
        os.chdir(folder)
        info = {}

        try:
            encut_files = []
            encut_values = []
            encut_based_energies = []
            for i in glob.glob("*.json"):
                if "ENCUT" in i:
                    encut_values.append(
                        int(str(i.split("-")[-1]).split(".json")[0])
                    )
                    encut_files.append(i)
                    erun = os.path.join(
                        folder, i.split(".json")[0], "vasprun.xml"
                    )
                    energy = float(Vasprun(erun).final_energy)
                    encut_based_energies.append(energy)

            kplength_files = []
            kp_values = []
            kp_based_energies = []
            for i in glob.glob("*.json"):
                if "KPOINT" in i:
                    kp_values.append(
                        int(str(i.split("-")[-1]).split(".json")[0])
                    )
                    kplength_files.append(i)
                    krun = os.path.join(
                        folder, i.split(".json")[0], "vasprun.xml"
                    )
                    energy = float(Vasprun(krun).final_energy)
                    kp_based_energies.append(energy)
            encut_values = np.array(encut_values)
            encut_based_energies = np.array(encut_based_energies)
            order = np.argsort(encut_values)
            encut_values = encut_values[order].tolist()
            encut_based_energies = encut_based_energies[order].tolist()
            kp_values = np.array(kp_values)
            kp_based_energies = np.array(kp_based_energies)
            order = np.argsort(kp_values)
            kp_values = kp_values[order].tolist()
            kp_based_energies = kp_based_energies[order].tolist()

            info = {}
            info["converged_encut"] = str(sorted(list(set(encut_values)))[-6])
            info["converged_kpoint_length"] = str(
                sorted(list(set(kp_values)))[-6]
            )
            info["kp_values"] = "'" + ",".join(map(str, kp_values)) + "'"
            info["kp_based_energies"] = (
                "'" + ",".join(map(str, kp_based_energies)) + "'"
            )
            info["encut_values"] = "'" + ",".join(map(str, encut_values)) + "'"
            info["encut_based_energies"] = (
                "'" + ",".join(map(str, encut_based_energies)) + "'"
            )
            self.encut_kp_info = info
        except Exception:
            print("Check the convergence folders are named properly")
            pass
        os.chdir(folder)

        return info

    def electronic_dos_info(self, vrun):
        """Get electronic dos information with , _ and ; seperators."""
        line = ""
        try:
            total_dos = vrun.total_dos
            # energies = np.array(total_dos[0])
            spdf_dos = vrun.get_spdf_dos()
            atom_dos = vrun.get_atom_resolved_dos()

            line = ""
            for i, ii in enumerate(total_dos):
                if i == 0:
                    line = (
                        line
                        + "<edos_energies>'"
                        + ",".join(map(str, np.array(ii) - vrun.efermi))
                        + "'</edos_energies>"
                    )
                elif i == 1:
                    line += (
                        "<total_edos_up>'"
                        + ",".join(map(str, ii))
                        + "'</total_edos_up>"
                    )
                elif i == 2:
                    line += (
                        "<total_edos_down>'"
                        + ",".join(map(str, ii))
                        + "'</total_edos_down>"
                    )

            line += "<spdf_dos>"
            for i, j in spdf_dos.items():
                line += (
                    "<"
                    + str(i)
                    + ">"
                    # + ">'"
                    + "'"
                    + ",".join(map(str, j))
                    + "'"
                    + "</"
                    # + "'</"
                    + str(i)
                    + ">"
                )

            line += "</spdf_dos>"
            line += "<elemental_dos>"
            spin_up_info = atom_dos["spin_up_info"]
            spin_down_info = atom_dos["spin_down_info"]
            line += '<spin_up_info>"'
            for i, j in spin_up_info.items():
                line += str(i) + str("_") + ",".join(map(str, j)) + ";"
            line += '"</spin_up_info>'

            line += '<spin_down_info>"'
            # line += '<spin_down_info>"'
            for i, j in spin_down_info.items():
                line += str(i) + str("_") + ",".join(map(str, j)) + ";"
            line += '"</spin_down_info>'
            # line += '"</spin_down_info>'
            line += "</elemental_dos>"

            fermi_velocities = ""
            try:
                fermi_velocities = ",".join(map(str, vrun.fermi_velocities[0]))
            except Exception:
                print("Cannot get the dos Fermi-velocity.")
                pass
            line += (
                '<fermi_velocities>"'
                + fermi_velocities
                + '"</fermi_velocities>'
            )

        except Exception:
            print("Cannot get DOS data", vrun)
        return line

    def get_xrd(self, atoms):
        """Get XRD pattern."""
        line = ""
        try:
            two_thetas, d_hkls, intensities = XRD().simulate(atoms=atoms)
            line += (
                '<two_thetas>"'
                + array_to_string(two_thetas)
                + '"</two_thetas>'
            )
            line += '<d_hkls>"' + array_to_string(d_hkls) + '"</d_hkls>'
            line += (
                '<intensities>"' + array_to_string(d_hkls) + '"</intensities>'
            )
        except Exception:
            print("Cannot gt XRD pattern.")
            pass
        return line

    def electronic_band_struct(self, vasprun="", kpoints_file_path=""):
        """Get electronic_band_struct."""
        line = self.ebandstruct(vrun=vasprun, kp=kpoints_file_path)
        return line

    def main_band(self):
        """Get bandstructure from MAIN-BAND-*."""
        folder = self.folder
        os.chdir(folder)
        data = ""
        info = {}
        for i in glob.glob("MAIN-BAND*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                kp = os.getcwd() + "/" + folder + "/KPOINTS"
                data = self.electronic_band_struct(
                    vasprun=vrun, kpoints_file_path=kp
                )
            except Exception:
                print("Cannot get ebands info", folder)
                pass
        os.chdir(folder)
        info["main_bands_info"] = data
        return info

    def main_hse06_band(self):
        """Get bandstructure from MAIN-HSE-*."""
        folder = self.folder
        os.chdir(folder)
        data = ""
        info = {}
        for i in glob.glob("MAIN-HSE*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                kp = os.getcwd() + "/" + folder + "/KPOINTS"
                data = self.electronic_band_struct(
                    vasprun=vrun, kpoints_file_path=kp
                )
            except Exception:
                print("Cannot get hseband info", folder)
                pass
        os.chdir(folder)
        info["main_hse_bands_info"] = data
        return info

    def main_pbe0_band(self):
        """Get bandstructure from MAIN-PBE0-*."""
        folder = self.folder
        os.chdir(folder)
        data = ""
        info = {}
        for i in glob.glob("MAIN-PBE0*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                kp = os.getcwd() + "/" + folder + "/KPOINTS"
                data = self.electronic_band_struct(
                    vasprun=vrun, kpoints_file_path=kp
                )
            except Exception:
                print("Cannot get ebands info", folder)
                pass
        os.chdir(folder)
        info["main_pbe0_bands_info"] = data
        return info

    def get_wannier_lines(self, data={}, energy_tol=4):
        """Get data for wannier- quality check."""
        # Dense mesh
        line = ""
        mesh_vrun_x = []
        mesh_vrun_y = []
        line += (
            "<maxdiff_mesh>"
            + str(data["info_mesh"]["maxdiff"])
            + "</maxdiff_mesh>"
        )
        line += (
            "<maxdiff_bz>" + str(data["info_bz"]["maxdiff"]) + "</maxdiff_bz>"
        )
        for i, ii in enumerate(np.array(data["info_mesh"]["eigs_vrun"]).T):
            y = ",".join(map(str, ii))
            x = ",".join(map(str, range(0, len(ii))))
            mesh_vrun_y.append(y)
            mesh_vrun_x.append(x)
        line += (
            '<mesh_vrun_x>"'
            + ";".join(map(str, mesh_vrun_x))
            + '"</mesh_vrun_x>'
        )
        line += (
            '<mesh_vrun_y>"'
            + ";".join(map(str, mesh_vrun_y))
            + '"</mesh_vrun_y>'
        )

        min_arr = []
        erange = [-energy_tol, energy_tol]
        dd = {}
        kpp = data["info_mesh"]["eigs_vrun"]
        vasp = np.array(data["info_mesh"]["eigs_vrun"])
        wann = np.array(data["info_mesh"]["eigs_wan"])
        for k in range(len(kpp)):
            for n in wann[k]:
                diff_arr = []
                if n > erange[0] and n < erange[1]:
                    for v in vasp[k]:
                        diff = abs(n - v)
                        diff_arr.append(diff)
                if diff_arr != []:
                    tmp = np.min(diff_arr)
                    dd.setdefault(n, tmp)
                    min_arr.append(tmp)

        line += (
            '<mesh_difference_values>"'
            + str(",".join(map(str, dd.values())))
            + '"</mesh_difference_values>'
        )
        line += (
            '<mesh_difference_keys>"'
            + str(",".join(map(str, dd.keys())))
            + '"</mesh_difference_keys>'
        )

        mesh_wann_x = []
        mesh_wann_y = []
        for i, ii in enumerate(np.array(data["info_mesh"]["eigs_wan"]).T):
            y = ",".join(map(str, ii))
            x = ",".join(map(str, range(0, len(ii))))
            mesh_wann_y.append(y)
            mesh_wann_x.append(x)
        line += (
            '<mesh_wann_x>"'
            + ";".join(map(str, mesh_wann_x))
            + '"</mesh_wann_x>'
        )
        line += (
            '<mesh_wann_y>"'
            + ";".join(map(str, mesh_wann_y))
            + '"</mesh_wann_y>'
        )

        # High symmetry BZ k-points

        bz_vrun_x = []
        bz_vrun_y = []
        for i, ii in enumerate(np.array(data["info_bz"]["eigs_vrun"]).T):
            y = ",".join(map(str, ii))
            x = ",".join(map(str, range(0, len(ii))))
            bz_vrun_y.append(y)
            bz_vrun_x.append(x)
        line += (
            '<bz_vrun_x>"' + ";".join(map(str, bz_vrun_x)) + '"</bz_vrun_x>'
        )
        line += (
            '<bz_vrun_y>"' + ";".join(map(str, bz_vrun_y)) + '"</bz_vrun_y>'
        )

        min_arr = []
        erange = [-energy_tol, energy_tol]
        dd = {}
        kpp = data["info_bz"]["eigs_vrun"]
        vasp = np.array(data["info_bz"]["eigs_vrun"])
        wann = np.array(data["info_bz"]["eigs_wan"])
        for k in range(len(kpp)):
            for n in wann[k]:
                diff_arr = []
                if n > erange[0] and n < erange[1]:
                    for v in vasp[k]:
                        diff = abs(n - v)
                        diff_arr.append(diff)
                if diff_arr != []:
                    tmp = np.min(diff_arr)
                    dd.setdefault(n, tmp)
                    min_arr.append(tmp)

        line += (
            '<bz_difference_values>"'
            + str(",".join(map(str, dd.values())))
            + '"</bz_difference_values>'
        )
        line += (
            '<bz_difference_keys>"'
            + str(",".join(map(str, dd.keys())))
            + '"</bz_difference_keys>'
        )
        return line

    def wannier_comparison_plot(self, energy_tol=4):
        """Compare Wannier and DFT data on mesh and BZ."""
        folder = self.folder
        os.chdir(folder)
        line = ""
        info = {}

        for i in glob.glob("MAIN-WANN*"):
            try:
                folder = i
                wann_dat = os.getcwd() + "/" + folder + "/wannier90_hr.dat"
                if os.path.exists(wann_dat):
                    wann_json = wann_dat = (
                        os.getcwd() + "/" + folder + "/wannier_comparison.json"
                    )
                    if os.path.exists(wann_json):
                        print("Getting data from preexisiting wannier json.")
                        data = loadjson(wann_json)
                        line += self.get_wannier_lines(data=data, energy_tol=4)
                    else:
                        print("Running jarvis-wannier solver.")
                        wann_ham = WannierHam(filename=wann_dat)
                        bands_vrun = wann_dat.replace(
                            "wannier90_hr.dat", "vasprun.xml"
                        ).replace("WANN", "SOCSCFBAND")
                        mesh_vrun = wann_dat.replace(
                            "wannier90_hr.dat", "vasprun.xml"
                        ).replace("WANN", "SOC")
                        bz = wann_ham.compare_dft_wann(
                            vasprun_path=bands_vrun, plot=False
                        )
                        mesh = wann_ham.compare_dft_wann(
                            vasprun_path=mesh_vrun, plot=False
                        )
                        data = {}
                        data["info_bz"] = bz
                        data["info_mesh"] = mesh
                        line += self.get_wannier_lines(data=data, energy_tol=4)

            except Exception:
                print("Cannot get wannier data.", folder)
                pass
        info["wannier_band_comparison"] = line
        return info

    def vacancy_formation_optb88vdw(self):
        """Get vacancy formation energy data."""
        folder = self.folder
        os.chdir(folder)

        line = ""
        info = {}

        try:
            vac_path = os.path.join(folder, "MAIN-VACANCY")
            if os.path.exists(vac_path):
                for i in glob.glob(vac_path + "/JVASP*"):
                    vac_vrun_path = os.path.join(
                        folder, "MAIN-VACANCY", i, "vasprun.xml"
                    )
                    outcar = Outcar(
                        os.path.join(folder, "MAIN-VACANCY", i, "OUTCAR")
                    )
                    if outcar.converged:
                        vac_vrun = Vasprun(vac_vrun_path)
                        tmp = i.split("/")[-1].split("_")
                        jid = tmp[0]
                        atom = tmp[2]
                        Ef = get_twod_defect_energy(
                            vrun=vac_vrun, jid=jid, atom=atom
                        )
                        line += (
                            jid + "_" + atom + "_" + str(round(Ef, 3)) + ","
                        )
        except Exception:
            print("Cannot get vacancy info.", folder)
            pass
        info["vacancy_formation_energy"] = line
        os.chdir(folder)
        return info

    def loptics_optoelectronics(self, vrun=""):
        """Get optoelctronic data."""
        line = ""

        try:
            lvrun = Vasprun(vrun)
            reals, imags = lvrun.dielectric_loptics
            energies = reals[:, 0]
            line = ""
            line += (
                "<energies>'" + ",".join(map(str, energies)) + "'</energies>"
            )
            line += '<loptics_dielectric_constant>"'
            line += str(reals[:, 1][0]) + ","
            line += str(reals[:, 2][0]) + ","
            line += str(reals[:, 3][0]) + ","
            line += str(reals[:, 4][0]) + ","
            line += str(reals[:, 5][0]) + ","
            line += str(reals[:, 6][0])
            line += '"</loptics_dielectric_constant>'
            lopt = ""
            try:
                max_linopt_eps = max(
                    reals[:, 1][0],
                    reals[:, 2][0],
                    reals[:, 3][0],
                    reals[:, 4][0],
                    reals[:, 5][0],
                    reals[:, 6][0],
                )
                lopt = str(round(max_linopt_eps, 2))
            except Exception:
                pass

            line += "<max_linopt_eps>" + lopt + "</max_linopt_eps>"
            fermi_velocities = ""
            try:
                fermi_velocities = ",".join(
                    map(str, lvrun.fermi_velocities[0])
                )
            except Exception:
                print("Cannot get the lepsilon Fermi-velocity.", vrun)
                pass
            line += (
                '<fermi_velocities>"'
                + fermi_velocities
                + '"</fermi_velocities>'
            )
            for i in np.arange(1, reals.shape[1]):
                line += (
                    "<real_"
                    + str(i)
                    + ">'"
                    + ",".join(map(str, reals[:, i]))
                    + "'</real_"
                    + str(i)
                    + ">"
                )
                line += (
                    "<imag_"
                    + str(i)
                    + ">'"
                    + ",".join(map(str, imags[:, i]))
                    + "'</imag_"
                    + str(i)
                    + ">"
                )
        except Exception:
            print("Cannot get loptics data", vrun)
            pass
        eff_slme = ""
        eff_sq = ""
        dirgap = ""
        indirgap = ""
        try:
            dirgap = round(lvrun.get_dir_gap, 3)
            indirgap = round(lvrun.get_indir_gap[0], 3)
            en, abz = lvrun.avg_absorption_coefficient
            abz = abz * 100
            eff_slme = SolarEfficiency().slme(
                en, abz, indirgap, indirgap, plot_current_voltage=False
            )
            # print("SLME", 100 * eff)
            eff_sq = SolarEfficiency().calculate_SQ(indirgap)
            eff_slme = round(100 * eff_slme, 2)
            eff_sq = round(100 * eff_sq, 2)
        except Exception:
            print("Cannot get solar data.")
            pass
        line += "<solar_slme>" + str(eff_slme) + "</solar_slme>"
        line += "<solar_sq>" + str(eff_sq) + "</solar_sq>"
        line += "<opto_dir_gap>" + str(dirgap) + "</opto_dir_gap>"
        line += "<opto_indir_gap>" + str(indirgap) + "</opto_indir_gap>"

        return line

    def main_optics_semilocal(self):
        """Get optoelctronic data from MAIN-OPTICS*."""
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-OPTICS*.json"):

            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                data = self.loptics_optoelectronics(vrun)
            except Exception:
                print("Cannot get semilocal optics data", folder)
                pass
        os.chdir(folder)
        info["main_optics_info"] = data
        self.main_optics_semilocal_info = data
        return info

    def main_optics_mbj(self):
        """Get optoelctronic data from MAIN-MBJ*."""
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-MBJ*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                data = self.loptics_optoelectronics(vrun)
            except Exception:
                print("Cannot get mbj optics data", folder)
                pass
        os.chdir(folder)
        info["main_optics_mbj_info"] = data
        self.main_optics_mbj_info = data
        return info

    def spillage(
        self,
        soc_wf="MAIN-SOCSCFBAND-bulk@JVASP-1067_mp-541837/WAVECAR",
        nonsoc_wf="MAIN-MAGSCFBAND-bulk@JVASP-1067_mp-541837/WAVECAR",
        spl_json_existing="/rk2/knc6/DB/SPILLAGE/all_spillage.json",
        jid="",
    ):
        """Get spillage data."""
        line = ""

        try:
            sp = {}
            print("soc_wf", soc_wf)
            tmp_jid = jid
            print("tmp_jid", tmp_jid)
            if tmp_jid != "":
                spl_dat = loadjson(spl_json_existing)
                for i in spl_dat:
                    if i["jid"] == tmp_jid:
                        sp = i
                        sp["kpoints"] = list(
                            [[0, 0, 0] for i in range(len(sp["spillage_k"]))]
                        )
            if sp == {}:
                # As spillage calculation may take a long time
                sp = Spillage(
                    wf_noso=nonsoc_wf, wf_so=soc_wf
                ).overlap_so_spinpol()
            soc_bands = self.electronic_band_struct(
                vasprun=soc_wf.replace("WAVECAR", "vasprun.xml"),
                kpoints_file_path=(soc_wf.replace("WAVECAR", "KPOINTS")),
            )
            line += "<soc_bands>" + soc_bands + "</soc_bands>"
            nonsoc_bands = self.electronic_band_struct(
                vasprun=nonsoc_wf.replace("WAVECAR", "vasprun.xml"),
                kpoints_file_path=(nonsoc_wf.replace("WAVECAR", "KPOINTS")),
            )
            line += "<nonsoc_bands>" + nonsoc_bands + "</nonsoc_bands>"
            spillage_k = ",".join(map(str, sp["spillage_k"]))
            line += "<spillage_k>'" + spillage_k + "'</spillage_k>"
            spillage_kpoints = ",".join(
                [";".join(map(str, i)) for i in sp["kpoints"]]
            )
            line += (
                "<spillage_kpoints>'"
                + spillage_kpoints
                + "'</spillage_kpoints>"
            )
            max_spillage = round(max(sp["spillage_k"]), 3)
            line += "<max_spillage>" + str(max_spillage) + "</max_spillage>"
        except Exception:
            print("Cannot get spillage info", soc_wf, nonsoc_wf)
            pass
        return line

    def main_soc_spillage(self, jid=""):
        """Get spillage data from folder."""
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-SOCSCFBAND*.json"):
            try:
                folder = i.split(".json")[0]
                soc_wf = os.getcwd() + "/" + folder + "/WAVECAR"
                nonsoc_wf = (
                    os.getcwd()
                    + "/"
                    + folder.replace("SOC", "MAG")
                    + "/WAVECAR"
                )
                data = self.spillage(
                    soc_wf=soc_wf, nonsoc_wf=nonsoc_wf, jid=jid
                )
            except Exception:
                print("Cannot get spillage info", folder)
                pass
        os.chdir(folder)
        info["main_spillage_info"] = data
        self.main_soc_spillage_info = info
        return info

    def basic_info(
        self,
        include_dos_info=True,
        include_neighbor_info=True,
        use_cfid_data=True,
    ):
        """Get data from MAIN-RELAX folder."""
        main_folder = self.folder
        os.chdir(main_folder)
        info = {}
        info["source_folder"] = self.folder
        for i, j in self.meta_data.items():
            info[i] = j
        id_file = self.meta_data["id_file"]
        for i in glob.glob("MAIN-RELAX*.json"):
            folder = i.split(".json")[0]
            main_vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
            # vrun = Vasprun(main_vrun)
            # atoms = vrun.all_structures[-1]
            try:
                vrun = Vasprun(main_vrun)
                atoms = vrun.all_structures[-1]
                info["XRD"] = self.get_xrd(atoms)
                f = open(os.path.join(folder, "..", id_file), "r")
                lines = f.read().splitlines()
                f.close()
                id = lines[0]
                info["id"] = id
                formula = atoms.composition.reduced_formula
                info["formula"] = formula
                elements = ",".join(atoms.uniq_species)
                info["elements"] = elements
                info["number_uniq_species"] = len(atoms.uniq_species)
                method = ""
                try:
                    f = open(os.path.join(folder, "..", "FUNCTIONAL"), "r")
                    lines = f.read().splitlines()
                    f.close()
                    method = lines[0]
                except Exception:
                    print('Cannot find "FUNCTIONAL"')
                    pass
                if method == "PBEBO":
                    method = "OptB88vdW"
                if method == "PBEOR":
                    method = "OptPBEvdW"
                if method == "PBEMK":
                    method = "OptB86bvdW"
                info["method"] = method
                final_energy = vrun.final_energy
                num_atoms = atoms.num_atoms
                fen = ""
                if method == "OptB88vdW":
                    fen = form_enp(atoms=atoms, total_energy=final_energy)
                info["formation_energy"] = fen
                rel_en = final_energy / num_atoms
                info["relaxed_energy"] = round(rel_en, 3)
                spg = Spacegroup3D(atoms)
                prim_atoms = spg.primitive_atoms
                conv_atoms = spg.conventional_standard_structure
                spg_numb = spg.space_group_number
                info["spacegroup_number"] = spg_numb
                info["point_group_symbol"] = spg.point_group_symbol
                # spg_symb = spg.space_group_symbol
                info["spg_space_group_symbol"] = spg.space_group_symbol
                crys_system = spg.crystal_system
                info["crys_system"] = crys_system
                conv_params = conv_atoms.lattice.parameters
                a_conv = round(conv_params[0], 2)
                b_conv = round(conv_params[1], 2)
                c_conv = round(conv_params[2], 2)
                info["a_conv"] = a_conv
                info["b_conv"] = b_conv
                info["c_conv"] = c_conv
                alpha_conv = round(conv_params[3], 2)
                beta_conv = round(conv_params[4], 2)
                gamma_conv = round(conv_params[5], 2)
                info["alpha_conv"] = alpha_conv
                info["beta_conv"] = beta_conv
                info["gamma_conv"] = gamma_conv
                prim_params = prim_atoms.lattice.parameters
                a_prim = round(prim_params[0], 2)
                b_prim = round(prim_params[1], 2)
                c_prim = round(prim_params[2], 2)
                info["a_prim"] = a_prim
                info["b_prim"] = b_prim
                info["c_prim"] = c_prim
                alpha_prim = round(prim_params[3], 2)
                beta_prim = round(prim_params[4], 2)
                gamma_prim = round(prim_params[5], 2)
                info["alpha_prim"] = alpha_prim
                info["beta_prim"] = beta_prim
                info["gamma_prim"] = gamma_prim
                prim_natoms = prim_atoms.num_atoms
                conv_natoms = conv_atoms.num_atoms
                info["prim_natoms"] = prim_natoms
                info["conv_natoms"] = conv_natoms
                density = round(atoms.density, 3)
                volume = round(atoms.volume, 3)
                info["density"] = density

                info["volume"] = volume
                packing_fr = ""
                try:
                    packing_fr = round(atoms.packing_fraction, 3)
                except Exception:
                    print("Cannot calculate packing_fr, check atomic_radii")
                    pass
                info["packing_fr"] = packing_fr

                dim = get_supercell_dims(conv_atoms)
                # dim=[1,1,1]
                xyz = (
                    conv_atoms.make_supercell_matrix(dim)
                    .center_around_origin()
                    .get_xyz_string
                )
                info["xyz"] = '"' + str(xyz).replace("\n", "\\n") + '"'
                info["poscar_conv"] = (
                    '"'
                    + str(
                        conv_atoms.make_supercell_matrix(dim).get_string()
                    ).replace("\n", "\\n")
                    + '"'
                )

                info["contcar"] = (
                    '"' + str(atoms.get_string()).replace("\n", "\\n") + '"'
                )

                # info["xyz"] = '"' + str(xyz).replace("\n", "\\n") + '"'
                rdf_bins = np.arange(0.1, 10.2, 0.1)
                ang1_bins = np.arange(1, 181.0, 1)
                ang2_bins = np.arange(1, 181.0, 1)
                dhd_bins = np.arange(1, 181.0, 1)
                rdf_hist = []
                ang1_hist = []
                ang2_hist = []
                dhd_hist = []
                info["cfid_descs"] = ""
                # Use pre-computed CFID dataset with chemo-structural features
                if use_cfid_data:
                    try:
                        cfid_descs = get_cfid_descriptors(jid=id, atoms=atoms)
                        if cfid_descs is not None:
                            include_neighbor_info = False
                            rdf_hist = cfid_descs[820:920]
                            ang1_hist = cfid_descs[920:1099]
                            ang2_hist = cfid_descs[1099:1278]
                            dhd_hist = cfid_descs[1278:1457]
                            info["cfid_descs"] = (
                                "'" + ",".join(map(str, cfid_descs)) + "'"
                            )
                    except Exception:
                        print("Cannot obtain cfid dataset", id)
                        pass
                if include_neighbor_info:
                    nbr = NeighborsAnalysis(atoms)
                    rdf_bins, rdf_hist, nn = nbr.get_rdf()
                    nbr = NeighborsAnalysis(atoms, max_cut=10.0)
                    ang1_hist, ang1_bins = nbr.ang_dist_first()
                    ang2_hist, ang2_bins = nbr.ang_dist_second()
                    dhd_hist, dhd_bins = nbr.get_ddf()
                info["rdf_bins"] = "'" + ",".join(map(str, rdf_bins)) + "'"
                info["rdf_hist"] = "'" + ",".join(map(str, rdf_hist)) + "'"

                info["ang1_bins"] = "'" + ",".join(map(str, ang1_bins)) + "'"
                info["ang1_hist"] = "'" + ",".join(map(str, ang1_hist)) + "'"
                info["ang2_bins"] = "'" + ",".join(map(str, ang2_bins)) + "'"
                info["ang2_hist"] = "'" + ",".join(map(str, ang2_hist)) + "'"
                info["dihedral_bins"] = (
                    "'" + ",".join(map(str, dhd_bins)) + "'"
                )
                info["dihedral_hist"] = (
                    "'" + ",".join(map(str, dhd_hist)) + "'"
                )

                scf_indir_gap = vrun.get_indir_gap[0]
                scf_dir_gap = vrun.get_dir_gap
                info["scf_indir_gap"] = round(scf_indir_gap, 3)
                info["scf_dir_gap"] = round(scf_dir_gap, 3)
                try:
                    oszicar = Oszicar(os.path.join(folder, "OSZICAR"))
                    magmom = oszicar.magnetic_moment
                    info["magmom"] = round(float(magmom), 3)
                except Exception:
                    print("Check oszicar")
                    pass
                if include_dos_info:
                    main_relax_dos = self.electronic_dos_info(vrun)
                    info["main_relax_dos"] = main_relax_dos
            except Exception:
                pass
        os.chdir(main_folder)
        self.basic_info_dict = info
        return info

    def dfpt_related(self, vrun="", out=""):
        """Get DFPT data."""
        vrun = Vasprun(vrun)
        data = vrun.dfpt_data
        out = Outcar(out)
        line = ""
        line += (
            "<epsilon>'"
            + ",".join(
                [
                    ";".join(map(str, data["epsilon"]["epsilon"][:, i]))
                    for i in range(0, 3)
                ]
            )
            + "'</epsilon>"
        )
        line += (
            "<epsilon_rpa>"
            + ",".join(map(str, data["epsilon"]["epsilon_rpa"].flatten()))
            + "</epsilon_rpa>"
        )
        line += (
            "<epsilon_ion>"
            + ",".join(map(str, data["epsilon"]["epsilon_ion"].flatten()))
            + "</epsilon_ion>"
        )
        max_eps = np.max(np.abs(np.array(data["epsilon"]["epsilon"])))
        line += "<max_eps>" + str(max_eps) + "</max_eps>"
        born_charges = data["born_charges"]
        atoms = vrun.all_structures[-1]
        spgg = Spacegroup3D(atoms)
        wycs = spgg._dataset["wyckoffs"]
        # spg = str(spgg._dataset["number"])
        natoms = atoms.num_atoms
        combs = []
        line += "<born_effective_charge>"
        for k in range(natoms):
            comb = str(atoms.elements[k]) + "," + str(wycs[k])
            if comb not in combs:
                p = (
                    str(atoms.elements[k])
                    + ","
                    + str(wycs[k])
                    + ","
                    + ",".join(map(str, born_charges[k].flatten()))
                    + ";"
                )
                line += p
                combs.append(comb)
        line += "</born_effective_charge>"
        # print (line)
        # vrun_eigs = data["phonon_eigenvalues"]
        pza = out.piezoelectric_tensor[1]
        max_pza = np.max(np.abs(pza))
        line += (
            "<max_piezo_stress_coeff>"
            + str(max_pza)
            + "</max_piezo_stress_coeff>"
        )
        pz = ",".join([";".join(map(str, pza[:, i])) for i in range(0, 6)])
        line += (
            '<dfpt_piezoelectric_tensor>"'
            + pz
            + '"</dfpt_piezoelectric_tensor>'
        )

        phonon_eigenvalues = out.phonon_eigenvalues
        phonon_eigenvectors = data["phonon_eigenvectors"]
        masses = data["masses"]
        born_charges = data["born_charges"]
        x, y = ir_intensity(
            phonon_eigenvectors=phonon_eigenvectors,
            phonon_eigenvalues=phonon_eigenvalues,
            masses=masses,
            born_charges=born_charges,
        )

        line += "<ir_intensity>'"
        line += ",".join(map(str, x)) + ";" + ",".join(map(str, y))
        line += "'</ir_intensity>"

        xx, yy = ir_intensity(
            smoothen=False,
            phonon_eigenvectors=phonon_eigenvectors,
            phonon_eigenvalues=phonon_eigenvalues,
            masses=masses,
            born_charges=born_charges,
        )
        xx = np.array(xx)
        yy = np.array(yy)
        max_ir_mode = ""
        min_ir_mode = ""
        try:
            tol = 1e-5
            non_zero_modes = yy > tol
            max_ir_mode = max(yy[non_zero_modes])
            min_ir_mode = min(yy[non_zero_modes])
        except Exception:
            print("Cannot find max min IR mode.", out)
        line += "<max_ir_mode>" + str(round(max_ir_mode, 2)) + "</max_ir_mode>"
        line += "<min_ir_mode>" + str(round(min_ir_mode, 2)) + "</min_ir_mode>"

        return line

    def main_lepsilon(self):
        """Get DFPT data from folder."""
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""
        for i in glob.glob("MAIN-LEPSILON*.json"):
            try:
                folder = i.split(".json")[0]
                vrun = os.getcwd() + "/" + folder + "/vasprun.xml"
                out = os.getcwd() + "/" + folder + "/OUTCAR"
                data = self.dfpt_related(out=out, vrun=vrun)
            except Exception:
                print("Cannot find LEPSION data", folder)
                pass
        os.chdir(folder)
        info["main_lepsilon_info"] = data
        self.main_lepsilon_info = info
        return info

    def elastic_props(
        self, outcar="MAIN-ELASTIC-bulk@mp-149/OUTCAR",
    ):
        """Get elastic property data."""
        out = Outcar(outcar)
        try:
            out_phonon_eigenvalues = out.phonon_eigenvalues
        except Exception:
            print("Cannot get Outcar phonons.")
            pass
        vacuum = False
        atoms = None
        unit_system = "GPa"
        if "Layer" in self.meta_data["material_type"]:
            vacuum = True
            atoms = Atoms.from_poscar(outcar.replace("OUTCAR", "POSCAR"))
            unit_system = "Nm^-1"
        cij = np.array(out.elastic_props(atoms=atoms, vacuum=vacuum)["cij"])
        d = ElasticTensor(cij).to_dict()
        if "Layer" in self.meta_data["material_type"]:
            d = {}
        line = ""
        line += (
            '<cij>"'
            + ";".join([",".join(map(str, cij[:, i])) for i in range(0, 6)])
            + '"</cij>'
        )
        # print("line", line)
        for i, j in d.items():
            if i != "raw_et_tensor":
                line += (
                    "<" + str(i) + ">" + str(round(j, 2)) + "</" + str(i) + ">"
                )
        totdos = outcar.replace("OUTCAR", "total_dos.dat")
        line += '<unit_system>"' + unit_system + '"</unit_system>'
        cwd = str(os.getcwd())
        if not os.path.isfile(totdos):
            run_phonopy(outcar.split("/OUTCAR")[0])
        os.chdir(cwd)
        if os.path.isfile(totdos):

            mesh_yaml = outcar.replace("OUTCAR", "mesh.yaml")
            # band_yaml = outcar.replace("OUTCAR", "band.yaml")
            with open(mesh_yaml, "r") as f:
                doc = yaml.load(f)
            nmodes = doc["phonon"][0]["band"]
            ph_modes = []

            for p in nmodes:
                ph_modes.append(p["frequency"])
            ph_modes = sorted(set(ph_modes))
            f = open(totdos, "r")
            freq = []
            pdos = []
            for lines in f.readlines():
                if not str(lines.split()[0]).startswith("#"):
                    #   print (lines)
                    # else:
                    freq.append(float(lines.split()[0]))
                    pdos.append(float(lines.split()[1]))
            line += (
                "<phonon_modes>"
                + ",".join(map(str, [round(i, 2) for i in ph_modes]))
                + "</phonon_modes>"
            )
            line += (
                "<min_fd_phonon_mode>"
                + str(round(min(out_phonon_eigenvalues), 1))
                + "</min_fd_phonon_mode>"
            )
            line + '<outcar_phonons>"' + ",".join(
                map(str, out_phonon_eigenvalues)
            ) + '"</outcar_phonons>'
            line += (
                "<phonon_dos_frequencies>'"
                + ",".join(map(str, freq))
                + "'</phonon_dos_frequencies>"
            )
            line += (
                "<phonon_dos_intensity>'"
                + ",".join(map(str, pdos))
                + "'</phonon_dos_intensity>"
            )

            # Comment until 4 MB text size error
            # frequencies, distances, labels,label_points = bandstructure_plot(
            #    band_yaml
            # )
            # tmp = ""
            # for i in range(np.array(frequencies).shape[1]):
            #    tmp += ",".join(map(str, np.array(frequencies)[:, i])) + ";"
            # line += (
            #    "<phonon_bandstructure_distances>'"
            #    + ",".join(map(str, distances))
            #    + "'</phonon_bandstructure_distances>"
            # )
            # line += (
            #    "<phonon_bandstructure_frequencies>'"
            #    + tmp
            #    + "'</phonon_bandstructure_frequencies>"
            # )
            # line += (
            #    "<phonon_bandstructure_labels>'"
            #    + ",".join(map(str, labels))
            #    + "'</phonon_bandstructure_labels>"
            # )
            # line += (
            #    "<phonon_bandstructure_label_points>'"
            #    + ",".join(map(str, label_points))
            #    + "'</phonon_bandstructure_label_points>"
            # )
        return line

    def main_elastic(self):
        """Get elastic property data from folder."""
        folder = self.folder
        os.chdir(folder)
        info = {}
        data = ""

        for i in glob.glob("MAIN-ELAST*.json"):
            try:
                folder = i.split(".json")[0]
                out = os.getcwd() + "/" + folder + "/OUTCAR"
                data = self.elastic_props(out)
            except Exception:
                print("Cannot get ELAST folder", folder)
                pass
        os.chdir(folder)
        info["main_elastic_info"] = data
        self.main_elastic_info = info
        return info

    def raman_data(self):
        """Get Raman intensity."""
        folder = self.folder
        os.chdir(folder)
        info = {}
        line = ""

        for i in glob.glob("RAMANDIR*/vasp_raman.dat"):
            try:
                raman_file = os.path.join(os.getcwd(), i)
                ram = open(raman_file, "r")
                lines = ram.read().splitlines()
                ram.close()
                an = raman_file.replace("vasp_raman.dat", "POSCAR")
                pos = Atoms.from_poscar(an)
                if len(lines) == 3 * pos.num_atoms + 1:
                    data = parse_raman_dat(raman_file)
                    line += (
                        '<frequencies>"'
                        + array_to_string(data["freqs"])
                        + '"</frequencies>'
                    )
                    line += (
                        '<activity>"'
                        + array_to_string(data["activity"])
                        + '"</activity>'
                    )
                    line += (
                        '<alpha>"'
                        + array_to_string(data["alpha"])
                        + '"</alpha>'
                    )
                    line += (
                        '<beta2>"'
                        + array_to_string(data["beta2"])
                        + '"</beta2>'
                    )
                    line += (
                        '<indices>"'
                        + array_to_string(data["beta2"])
                        + '"</indices>'
                    )
                    line += (
                        "<max_raman_mode>"
                        + str(round(max(data["freqs"]), 2))
                        + "</max_raman_mode>"
                    )
            except Exception:
                print("Cannot process Raman data.")
                pass
        os.chdir(folder)
        info["raman_dat"] = line
        return info

    def effective_mass_data(self,):
        """Get effective-mass data."""
        folder = self.folder
        info = {}
        line = ""
        os.chdir(folder)

        try:
            trap = os.path.join(folder, "trap_info.json")
            data = loadjson(trap)
            avg = data["avg_mass"]
            ndat = avg["n"]["300"][0]["data"]
            pdat = avg["p"]["300"][0]["data"]
            line += (
                '<electron_mass_300K>"'
                + str(round(ndat[0][0], 2))
                + ","
                + str(round(ndat[1][1], 2))
                + ","
                + str(round(ndat[2][2], 2))
                + '"</electron_mass_300K>'
            )
            line += (
                '<hole_mass_300K>"'
                + str(round(pdat[0][0], 2))
                + ","
                + str(round(pdat[1][1], 2))
                + ","
                + str(round(pdat[2][2], 2))
                + '"</hole_mass_300K>'
            )
            line += (
                "<max_electron_mass_300K>"
                + str(round(max(ndat[0][0], ndat[1][1], ndat[2][2]), 2))
                + "</max_electron_mass_300K>"
            )
            line += (
                "<max_hole_mass_300K>"
                + str(round(max(pdat[0][0], pdat[1][1], pdat[2][2]), 2))
                + "</max_hole_mass_300K>"
            )

        except Exception:
            print("Cannot get trap info.")
            pass
        os.chdir(folder)
        info["effective_mass"] = line
        return info

    def boltztrap_data(
        self,
        path="MAIN-RELAX-bulk@mp-149/boltztrap",
        temperature=600,
        doping=1e20,
    ):
        """Get transport data."""
        info = {}
        cwd = str(os.getcwd())
        if not os.path.exists(path):
            run_boltztrap(path.split("/boltztrap")[0])
        os.chdir(cwd)
        all_data = BoltzTrapOutput(path).to_dict()
        # info["all_data"] = all_data
        small_p = all_data["condtens_fixdoping"]["p"][temperature]
        small_n = all_data["condtens_fixdoping"]["n"][temperature]
        for i, j in small_p.items():
            if j["N_cm3"] == doping:
                tmp = j
                pseeb = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(
                            tmp["seeb"].reshape(3, 3) * 1e6
                        )
                    ]
                )
                pcond = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(tmp["cond"].reshape(3, 3))
                        / 1e14
                    ]
                )
                ppf = pseeb ** 2 * pcond / 1e6
                pkappa = np.linalg.eigvals(tmp["kappa"].reshape(3, 3))

                info["pseeb"] = [round(k, 2) for k in pseeb.tolist()]
                info["pcond"] = [round(k, 2) for k in pcond.tolist()]
                info["ppf"] = [round(k, 2) for k in ppf.tolist()]
                info["pkappa"] = [round(k, 2) for k in pkappa.tolist()]

        for i, j in small_n.items():
            if j["N_cm3"] == -1 * doping:
                tmp = j
                nseeb = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(
                            tmp["seeb"].reshape(3, 3) * 1e6
                        )
                    ]
                )
                ncond = np.array(
                    [
                        np.real(r)
                        for r in np.linalg.eigvals(tmp["cond"].reshape(3, 3))
                        / 1e14
                    ]
                )
                npf = nseeb ** 2 * ncond / 1e6
                nkappa = np.linalg.eigvals(tmp["kappa"].reshape(3, 3))

                info["nseeb"] = np.array(
                    [round(np.real(r), 2) for r in nseeb.tolist()]
                )
                info["ncond"] = np.array(
                    [round(np.real(r), 2) for r in ncond.tolist()]
                )
                info["npf"] = [round(k, 2) for k in npf.tolist()]
                info["nkappa"] = [round(k, 2) for k in nkappa.tolist()]
        line = ""
        for i, j in info.items():
            if isinstance(j, list):
                tmp = ";".join(map(str, j))
            else:
                tmp = j
            line += (
                "<"
                + str(i)
                + ">'"
                + ",".join(map(str, j))
                + "'</"
                + str(i)
                + ">"
            )
        # print (line)
        # print (all_data)
        return line

    def main_boltz_data(self):
        """Get transport data from folder."""
        folder = self.folder
        info = {}
        data = ""
        os.chdir(folder)
        for i in glob.glob("MAIN-RELAX*.json"):
            try:
                folder = i.split(".json")[0]
                path = os.getcwd() + "/" + folder + "/boltztrap"
                data = self.boltztrap_data(path=path)
            except Exception:
                print("Cannot find boltztrap folder", folder)
                pass
        os.chdir(folder)
        info["boltztrap_info"] = data
        self.boltztrap_info_dict = info
        return info

    def image_to_string(
        self, img_path="2DSTM/PNG_JSON3/JVASP-60776_mp-19795_pos.jpg",
    ):
        """Transform image to string."""
        # 2D array only
        fig = imread(img_path)
        line = (
            str(fig.shape[0])
            + "_"
            + str(fig.shape[1])
            + "_"
            + ",".join(map(str, fig[:, :, 0].flatten()))
        )
        # Comment when 4 MBlimit is over
        line = ""
        return line

    def main_stm_neg(self):
        """Get neg. bias constant height STM image."""
        folder = self.folder
        line = ""
        for i in glob.glob("MAIN-STM-NEG*.json"):
            try:
                folder = i.split(".json")[0]
                pchg = os.getcwd() + "/" + folder + "/PARCHG"
                TH_STM = TersoffHamannSTM(chg_name=pchg)
                filename = os.getcwd() + "/" + folder + "/testh.jpg"
                t1 = TH_STM.constant_height(filename=filename)
                zcut = t1["zcut"]
                print(zcut)
                # line += "<zcut>" + str(zcut) + "</zcut>"
                line += self.image_to_string(img_path=filename)
                cmd = "rm -rf " + filename
                if os.path.isfile(filename):
                    os.system(cmd)
            except Exception:
                print("Couldnt make negative bias STM image", folder)
                pass
        return line

    def main_stm_pos(self):
        """Get pos. bias constant height STM image."""
        folder = self.folder
        line = ""
        for i in glob.glob("MAIN-STM-POS*.json"):
            try:
                folder = i.split(".json")[0]
                pchg = os.getcwd() + "/" + folder + "/PARCHG"
                TH_STM = TersoffHamannSTM(chg_name=pchg)
                filename = os.getcwd() + "/" + folder + "/testh.jpg"
                t1 = TH_STM.constant_height(filename=filename)
                zcut = t1["zcut"]
                print(zcut)
                # line += "<zcut>" + str(zcut) + "</zcut>"
                line += self.image_to_string(img_path=filename)
                cmd = "rm -rf " + filename
                if os.path.isfile(filename):
                    os.system(cmd)
            except Exception:
                print("Couldnt make positive bias STM image", folder)
                pass
        return line

    def efg_tensor(self):
        """Get electric field gradient tensor."""
        main_folder = self.folder
        info = {}
        efg_dat = ""
        max_efg = ""
        os.chdir(main_folder)
        for i in glob.glob("MAIN-LEFG*.json"):
            try:
                folder = i.split(".json")[0]
                out = os.getcwd() + "/" + folder + "/OUTCAR"
                atoms = Atoms.from_poscar(out.replace("OUTCAR", "POSCAR"))
                out = Outcar(out)
                spgg = Spacegroup3D(atoms)
                wycs = spgg._dataset["wyckoffs"]
                natoms = atoms.num_atoms
                combs = []
                efg = out.efg_raw_tensor
                efg_dat = ""
                max_efg = np.max(np.abs(np.array(efg)))
                for k in range(natoms):
                    comb = str(atoms.elements[k]) + "," + str(wycs[k])
                    if comb not in combs:
                        efg_dat = (
                            efg_dat
                            + comb
                            + ","
                            + ",".join(map(str, efg[k].flatten()))
                            + ";"
                        )
                        combs.append(comb)
            except Exception:
                print("Cannot get LEFG folder", folder)
                pass

        info["efg_raw_tensor"] = efg_dat
        info["max_efg"] = max_efg
        info["max_efg_eta"] = ""
        # self.efg_raw_tensor_info = info
        os.chdir(main_folder)
        return info

    def write_xml(self, filename="temp.xml"):
        """Get overall XML file."""
        with open(filename, "w") as f:
            line = '<?xml version="1.0" encoding="UTF-8"?>\n'
            line += '<?xml-stylesheet type="text/xsl" '
            line += 'href="jarvisdft.xsl"?>\n<basic_info>'
            f.write(str(line))
            line = (
                "<convergence_info>"
                + stringdict_to_xml(self.encut_kp())
                + "</convergence_info>"
                + "\n"
            )
            f.write(str(line))
            line = stringdict_to_xml(self.wannier_comparison_plot()) + "\n"
            f.write(str(line))

            line = stringdict_to_xml(self.vacancy_formation_optb88vdw()) + "\n"
            f.write(str(line))
            line = stringdict_to_xml(self.raman_data()) + "\n"
            f.write(str(line))
            b_info = self.basic_info()
            jid = b_info["id"]
            line = (
                "<main_relax_info>"
                + stringdict_to_xml(b_info)
                + "</main_relax_info>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_band>"
                + stringdict_to_xml(self.main_band(), enforce_string=False)
                + "</main_band>"
            )
            f.write(str(line))
            line = (
                "<main_hse06_band>"
                + stringdict_to_xml(
                    self.main_hse06_band(), enforce_string=False
                )
                + "</main_hse06_band>"
            )
            f.write(str(line))

            line = stringdict_to_xml(self.effective_mass_data()) + "\n"
            f.write(str(line))

            line = (
                "<main_pbe0_band>"
                + stringdict_to_xml(
                    self.main_pbe0_band(), enforce_string=False
                )
                + "</main_pbe0_band>"
            )
            f.write(str(line))
            line = (
                "<main_optics_semilocal>"
                + stringdict_to_xml(self.main_optics_semilocal())
                + "</main_optics_semilocal>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_optics_mbj>"
                + stringdict_to_xml(self.main_optics_mbj())
                + "</main_optics_mbj>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_elastic>"
                + stringdict_to_xml(self.main_elastic())
                + "</main_elastic>"
                + "\n"
            )
            f.write(str(line))
            line = (
                "<main_boltz>"
                + stringdict_to_xml(self.main_boltz_data())
                + "</main_boltz>"
                + "\n"
            )
            f.write(str(line))
            line = stringdict_to_xml(self.main_lepsilon()) + "\n"
            f.write(str(line))
            line = stringdict_to_xml(self.main_soc_spillage(jid=jid)) + "\n"
            f.write(str(line))
            line = (
                stringdict_to_xml(self.efg_tensor(), enforce_string=True)
                + "\n"
            )
            f.write(str(line))

            line = (
                "<main_stm_neg>" + str(self.main_stm_neg()) + "</main_stm_neg>"
            )
            f.write(str(line))
            line = (
                "<main_stm_pos>" + str(self.main_stm_pos()) + "</main_stm_pos>"
            )
            f.write(str(line))
            f.write("</basic_info>\n")
            """
            line = "<stm_image>" + image_to_string() + "</stm_image>"+ "\n"
            f.write(str(line))
            line ='</basic_info>'
            f.write(str(line))
            """


"""
if __name__ == "__main__":
    folder = "/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO"
    filename = "JVASP-1002.xml"
    VaspToApiXmlSchema(folder=folder).write_xml(filename=filename)

    folder = "/rk2/knc6/JARVIS-DFT/TE-bulk/mp-541837_bulk_PBEBO"
    filename = "JVASP-1067.xml"
    VaspToApiXmlSchema(folder=folder).write_xml(filename=filename)

    folder = "/rk2/knc6/JARVIS-DFT/2D-1L/POSCAR-mp-2815-1L.vasp_PBEBO"
    filename = "JVASP-664.xml"
    VaspToApiXmlSchema(folder=folder).write_xml(filename=filename)

    # directories = ["/rk2/knc6/JARVIS-DFT/Elements-bulkk/mp-149_bulk_PBEBO"]
    # filenames = ["JVASP-1002.xml"]
"""
